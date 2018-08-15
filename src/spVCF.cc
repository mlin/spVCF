#include "spVCF.h"
#include <vector>
#include <sstream>
#include <algorithm>
#include <cstdint>
#include <stdexcept>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <climits>
#include <assert.h>
#include "strlcpy.h"
#include "htslib/kseq.h"
#include "htslib/tbx.h"

using namespace std;

namespace spVCF {

// split s on delim & return strlen(s). s is damaged by side-effect
template<typename Out>
size_t split(char *s, char delim, Out result, uint64_t maxsplit = ULLONG_MAX) {
    string delims("\n");
    delims[0] = delim;
    char *cursor = s;
    char *token = strsep(&cursor, delims.c_str());
    char *last_token = token;
    uint64_t i = 0;
    while (token) {
        *(result++) = last_token = token;
        if (++i < maxsplit) {
            token = strsep(&cursor, delims.c_str());
        } else {
            if (cursor) {
                *result = last_token = cursor;
            }
            break;
        }
    }
    return last_token ? (last_token + strlen(last_token) - s) : 0;
}

template<typename Out>
size_t split(string& s, char delim, Out result, uint64_t maxsplit = ULLONG_MAX) {
    return split(&s[0], delim, result, maxsplit);
}

// because std::ostringstream is too slow :(
class OStringStream {
public:
    OStringStream(size_t initial_capacity)
        : buf_size_(initial_capacity)
        , cursor_(0)
        {
            buf_ = make_unique<char[]>(buf_size_+1);
            buf_[0] = 0;
        }
    OStringStream() : OStringStream(64) {}
    OStringStream(const OStringStream&) = delete;

    inline void Add(char c) {
        if (remaining() == 0) {
            grow();
            assert(remaining() > 0);
        }

        buf_[cursor_++] = c;
        buf_[cursor_] = 0;
    }
    inline OStringStream& operator<<(char c) {
        Add(c);
        return *this;
    }

    void Add(const char* s) {
        size_t rem = remaining();
        size_t len = strlcpy(&buf_[cursor_], s, rem+1);
        if (len <= rem) {
            cursor_ += len;
            assert(buf_[cursor_] == 0);
        } else {
            cursor_ += rem;
            assert(buf_[cursor_] == 0);
            grow(len - rem);
            return Add(s + rem);
        }
    }
    inline OStringStream& operator<<(const char* s) {
        Add(s);
        return *this;
    }
    inline OStringStream& operator<<(const std::string& s) {
        return *this << s.c_str();
    }

    inline const char* Get() const {
        assert(buf_[cursor_] == 0);
        return &buf_[0];
    }

    inline size_t Size() const {
        assert(buf_[cursor_] == 0);
        return cursor_;
    }

    void Clear() {
        buf_[0] = cursor_ = 0;
    }

private:
    inline size_t remaining() const {
        assert(buf_[cursor_] == 0);
        assert(cursor_ <= buf_size_);
        return buf_size_ - cursor_;
    }

    size_t grow(size_t hint = 0) {
        buf_size_ = max(2*buf_size_, hint+buf_size_);
        auto buf = make_unique<char[]>(buf_size_+1);
        memcpy(&buf[0], &buf_[0], cursor_);
        buf[cursor_] = 0;
        swap(buf_,buf);
    }

    // actual size of buf_ should always be buf_size_+1 to accommodate NUL
    // invariants: cursor_ <= buf_size_ && buf_[cursor_] == 0
    unique_ptr<char[]> buf_;
    size_t buf_size_, cursor_;
};

// Base class for encoder/decoder with common state & error-handling
class TranscoderBase : public Transcoder {
public:
    TranscoderBase() = default;
    TranscoderBase(const TranscoderBase&) = delete;

    transcode_stats Stats() override { return stats_; }

protected:
    void fail(const string& msg) {
        ostringstream ss;
        ss << "spvcf: " << msg << " (line " << line_number_ << ")";
        throw runtime_error(ss.str());
    }

    // state to be updated by derived classes
    uint64_t line_number_ = 0;
    transcode_stats stats_;
};

class EncoderImpl : public TranscoderBase {
public:
    EncoderImpl(uint64_t checkpoint_period, bool sparse, bool squeeze)
        : checkpoint_period_(checkpoint_period), sparse_(sparse), squeeze_(squeeze)
        {}
    EncoderImpl(const EncoderImpl&) = delete;
    const char* ProcessLine(char* input_line) override;

private:
    void Squeeze(const vector<char*>& line);

    uint64_t checkpoint_period_ = 0;
    bool sparse_ = true;
    bool squeeze_ = false;

    vector<string> dense_entries_; // main state memory
    string chrom_;
    uint64_t since_checkpoint_ = 0, checkpoint_pos_ = 0;

    OStringStream buffer_;
    vector<string> roundDP_table_;
};

const char* EncoderImpl::ProcessLine(char* input_line) {
    ++line_number_;
    // Pass through header lines
    if (*input_line == 0 || *input_line == '#') {
        return input_line;
    }
    ++stats_.lines;

    // Split the tab-separated line
    vector<char*> tokens;
    tokens.reserve(dense_entries_.size() + 9);
    size_t linesz = split(input_line, '\t', back_inserter(tokens));
    if (tokens.size() < 10) {
        fail("Invalid: fewer than 10 columns");
    }
    uint64_t N = tokens.size() - 9;

    if (dense_entries_.empty()){ // First line: allocate the dense entries
        dense_entries_.resize(N);
        stats_.N = N;
    } else if (dense_entries_.size() != N) { // Subsequent line -- check expected # columns
        for (int i = 9; i < tokens.size(); i++) {
            const string t = tokens[i];
            if (!t.empty() && t[0] == '"') {
                fail("Input seems to be sparse-encoded already");
            }
        }
        fail("Inconsistent number of samples");
    }

    if (squeeze_) {
        Squeeze(tokens);
    }

    buffer_.Clear();

    // Pass through first nine columns
    buffer_ << tokens[0];
    for (int i = 1; i < 9; i++) {
        buffer_ << '\t';
        if (i != 7 || !sparse_) {
            buffer_ << tokens[i];
        } else {
            // prepend spVCF_checkpointPOS to the INFO column, conveying the POS of the
            // last checkpoint (full dense row), which is useful for random
            // access & partial decoding of the file.
            //
            // TODO: create an appropriate header line for spVCF_checkpointPOS?
            string INFO = tokens[7];
            string newINFO;
            if (!INFO.empty() && INFO != ".") {
                newINFO = INFO;
                newINFO = ";" + newINFO;
            }
            newINFO = "spVCF_checkpointPOS=" + to_string(checkpoint_pos_) + newINFO;
            buffer_ << newINFO;
        }
    }

    if (!sparse_) {
        for (int i = 9; i < tokens.size(); i++) {
            buffer_ << '\t' << tokens[i];
        }
        return buffer_.Get();
    }

    uint64_t quote_run = 0; // current run-length of quotes across the row
    uint64_t sparse_cells = 0;
    // Iterate over the columns, compare each entry with the last entry
    // recorded densely.
    for (uint64_t s = 0; s < N; s++) {
        string& m = dense_entries_[s];
        const char* t = tokens[s+9];
        if (*t == '"') {
            fail("Input seems to be sparse-encoded already");
        }
        if (m.empty() || strcmp(m.c_str(), t) != 0) {
            // Entry doesn't match the last one recorded densely for this
            // column. Output any accumulated run of quotes in the current row,
            // then this new entry, and update the state appropriately.
            if (quote_run) {
                buffer_ << "\t\"";
                if (quote_run > 1) {
                    buffer_ << to_string(quote_run);
                }
                quote_run = 0;
                ++sparse_cells;
            }
            buffer_ << '\t' << t;
            ++sparse_cells;
            m = t;
        } else {
            // Entry matches; add to the current run of quotes
            quote_run++;
        }
    }
    // Output final run of quotes
    if (quote_run) {
        buffer_ << "\t\"";
        if (quote_run > 1) {
            buffer_ << to_string(quote_run);
        }
        ++sparse_cells;
    }

    // CHECKPOINT -- return a densely-encode row -- if we've switched to a new
    // chromosome OR we've hit the specified period OR we've passed half the
    // period and this line is mostly dense anyway.
    if (chrom_ != tokens[0] ||
        (checkpoint_period_ > 0 &&
         (since_checkpoint_ >= checkpoint_period_ ||
         (since_checkpoint_*2 >= checkpoint_period_ && sparse_cells*2 >= N)))) {
        buffer_.Clear();
        for (int t = 0; t < tokens.size(); t++) {
            if (t > 0) {
                buffer_ << '\t';
            }
            buffer_ << tokens[t];
            if (t >= 9) {
                dense_entries_[t-9] = tokens[t];
            }
            assert(tokens.size() == stats_.N+9);
        }
        since_checkpoint_ = 0;
        errno = 0;
        uint64_t POS = strtoull(tokens[1], nullptr, 10);
        if (errno) {
            fail("Couldn't parse POS");
        }
        if (chrom_ == tokens[0] && POS < checkpoint_pos_) {
            fail("input VCF not sorted (detected decreasing POS)");
        }
        checkpoint_pos_ = POS;
        chrom_ = tokens[0];
        ++stats_.checkpoints;
        return buffer_.Get();
    }
    ++since_checkpoint_;

    stats_.sparse_cells += sparse_cells;
    auto sparse_pct = 100*sparse_cells/N;
    if (sparse_pct <= 25) {
        ++stats_.sparse75_lines;
    }
    if (sparse_pct <= 10) {
        ++stats_.sparse90_lines;
    }
    if (sparse_pct <= 1) {
        ++stats_.sparse99_lines;
    }

    return buffer_.Get();
}

// Truncate cells to GT:DP, and round DP down to a power of two, if
//   - AD is present and indicates zero read depth for alternate alleles;
//   - VR is present and zero
// All cells (and the FORMAT specification) are reordered to begin with
// GT:DP, followed by any remaining fields.
//
// Each element of line is modified in-place.
void EncoderImpl::Squeeze(const vector<char*>& line) {
    if (roundDP_table_.empty()) {
        // precompute a lookup table of DP values, rounded down to a power of
        // two and stringified
        roundDP_table_.push_back("0");
        for (uint64_t DP = 1; DP < 10000; DP++) {
            uint64_t rDP = uint64_t(pow(2, floor(log2(DP))));
            assert(rDP <= DP && DP < rDP*2);
            roundDP_table_.push_back(to_string(rDP));
        }
    }

    // parse the FORMAT field
    vector<string> format;
    size_t formatsz = split(line[8], ':', back_inserter(format));

    // locate fields of interest
    if (format[0] != "GT") {
        fail("cells don't start with genotype (GT)");
    }
    int iDP = -1;
    auto pDP = find(format.begin(), format.end(), "DP");
    if (pDP != format.end()) {
        iDP = pDP - format.begin();
        assert (iDP > 0 && iDP < format.size());
    }
    int iAD = -1;
    auto pAD = find(format.begin(), format.end(), "AD");
    if (pAD != format.end()) {
        iAD = pAD - format.begin();
        assert (iAD > 0 && iAD < format.size());
    }
    int iVR = -1;
    auto pVR = find(format.begin(), format.end(), "VR");
    if (pVR != format.end()) {
        iVR = pVR - format.begin();
        assert (iVR > 0 && iVR < format.size());
    }

    // compute the new field order and update FORMAT
    vector<size_t> permutation;
    permutation.push_back(0);
    if (iDP >= 1) {
        permutation.push_back(iDP);
    }
    for (size_t i = 1; i < format.size(); i++) {
        if (i != iDP) {
            permutation.push_back(i);
        }
    }
    OStringStream new_format;
    new_format << "GT";
    for (const auto i : permutation) {
        if (i > 0) {
            new_format << ":" << format[i];
        }
    }
    assert(new_format.Size() <= formatsz);
    strcpy(line[8], new_format.Get());

    // reusable buffers to save allocations in upcoming inner loop
    OStringStream new_cell;
    vector<char*> entries;

    // proceed through all cells
    for (int s = 9; s < line.size(); s++) {
        entries.clear();
        // parse individual entries
        size_t cellsz = split(line[s], ':', back_inserter(entries));
        if (entries.empty()) {
            fail("empty cell");
        }

        // decide if conditions exist to truncate this cell to GT:DP
        bool truncate = false;
        if (iAD > 0 && entries.size() > iAD) {
            // does AD have any non-zero values after the first value?
            char *c = strchr(entries[iAD], ',');
            if (c) {
                for (; (*c == '0' || *c == ','); c++);
                if (*c == 0) {
                    truncate = true;
                }
            }
        }
        if (iVR > 0 && entries.size() >= iVR) {
            // is VR zero?
            if (strcmp(entries[iVR], "0") == 0) {
                truncate = true;
            }
        }

        // construct revised cell, beginning with GT:DP, then any remaining fields
        new_cell.Clear();
        new_cell << entries[0];
        if (iDP > 0) {
            assert(permutation[1] == iDP);
            if (entries.size() > iDP) {
                if (truncate) {
                    // round down the DP value
                    errno = 0;
                    uint64_t DP = strtoull(entries[iDP], nullptr, 10);
                    if (errno) {
                        fail("Couldn't parse DP");
                    }
                    new_cell << ':';
                    if (DP < roundDP_table_.size()) {
                        new_cell << roundDP_table_[DP];
                    } else {
                        uint64_t rDP = pow(2, floor(log2(DP)));
                        assert(rDP <= DP && DP < rDP*2);
                        new_cell << to_string(rDP);
                    }
                } else {
                    new_cell << ':' << entries[iDP];
                }
            } else {
                new_cell << ":.";
            }
        }
        if (truncate) {
            ++stats_.squeezed_cells;
        } else {
            // copy over remaining fields
            for (size_t i = 2; i < permutation.size(); i++) {
                new_cell << ':';
                if (entries.size() > permutation[i]) {
                    new_cell << entries[permutation[i]];
                } else {
                    new_cell << '.';
                }
            }
        }

        // write revised cell back into line.
        assert(new_cell.Size() <= cellsz);
        strcpy(line[s], new_cell.Get());
    }
}

unique_ptr<Transcoder> NewEncoder(uint64_t checkpoint_period, bool sparse, bool squeeze) {
    return make_unique<EncoderImpl>(checkpoint_period, sparse, squeeze);
}

class DecoderImpl : public TranscoderBase {
public:
    DecoderImpl() = default;
    DecoderImpl(const DecoderImpl&) = delete;
    const char* ProcessLine(char* input_line) override;

private:
    vector<string> dense_entries_;
    OStringStream buffer_;
};

const char* DecoderImpl::ProcessLine(char *input_line) {
    ++line_number_;
    // Pass through header lines
    if (*input_line == 0 || *input_line == '#') {
        return input_line;
    }
    ++stats_.lines;

    // Split the tab-separated line
    vector<char*> tokens;
    tokens.reserve(dense_entries_.size());
    split(input_line, '\t', back_inserter(tokens));
    if (tokens.size() < 10) {
        fail("Invalid project VCF: fewer than 10 columns");
    }

    // Figure out number of dense columns, the number of columns on the first line
    uint64_t N = dense_entries_.empty() ? (tokens.size() - 9) : dense_entries_.size();
    if (dense_entries_.empty()) {
        dense_entries_.resize(N);
        stats_.N = N;
    }
    assert(dense_entries_.size() == N);

    // Pass through first nine columns
    buffer_.Clear();
    buffer_ << tokens[0];
    for (int i = 1; i < 9; i++) {
        buffer_ << '\t';
        if (i != 7) {
            buffer_ << tokens[i];
        } else {
            // Strip the spVCF_checkpointPOS INFO field if present
            string INFO = tokens[7];
            if (INFO.substr(0, 20) == "spVCF_checkpointPOS=") {
                auto p = INFO.find(';');
                if (p != string::npos) {
                    buffer_ << INFO.substr(p+1);
                } else {
                    buffer_ << '.';
                }
            } else {
                buffer_ << INFO;
            }
        }
    }

    // Iterate over the sparse columns
    uint64_t sparse_cells = (tokens.size()-9), dense_cursor = 0;
    for (uint64_t sparse_cursor = 0; sparse_cursor < sparse_cells; sparse_cursor++) {
        const char* t = tokens[sparse_cursor+9];
        if (*t == 0) {
            fail("empty cell");
        } else if (*t != '"') {
            // Dense entry - remember it and copy it to the output
            // TODO: Perhaps fill QC fields with missing values (.) if they were squeezed out.
            //       The VCF spec does however say "Trailing fields can be dropped"
            if (dense_cursor >= N) {
                fail("Greater-than-expected number of columns implied by sparse encoding");
            }
            dense_entries_[dense_cursor++] = t;
            buffer_ << '\t' << t;
        } else {
            // Sparse entry - determine the run length
            uint64_t r = 1;
            if (t[1]) { // strlen(t) > 1
                errno = 0;
                auto s = strtoull(t+1, nullptr, 10);
                if (errno) {
                    fail("Undecodable sparse cell");
                }
                r = s;
            }
            // Output the implied run of entries from the remembered state
            if (dense_cursor + r > N) {
                ostringstream msg;
                msg << "Greater-than-expected number of columns implied by sparse encoding"
                    << " (expected N=" << N << ")";
                fail(msg.str());
            }
            for (uint64_t p = 0; p < r; p++) {
                if (dense_entries_[dense_cursor].empty()) {
                    fail("Missing preceding dense cells");
                }
                buffer_ << '\t' << dense_entries_[dense_cursor++];
            }
            assert(dense_cursor <= N);
        }
    }
    if (dense_cursor != N) {
        ostringstream msg;
        msg << "Unexpected number of columns implied by sparse encoding"
            << " (expected N=" << N << ", got " << dense_cursor << ")";
        fail(msg.str());
    }

    stats_.sparse_cells += sparse_cells;
    auto sparse_pct = 100*sparse_cells/N;
    if (sparse_pct <= 25) {
        ++stats_.sparse75_lines;
    }
    if (sparse_pct <= 10) {
        ++stats_.sparse90_lines;
    }
    if (sparse_pct <= 1) {
        ++stats_.sparse99_lines;
    }

    return buffer_.Get();
}

unique_ptr<Transcoder> NewDecoder() {
    return make_unique<DecoderImpl>();
}

class TabixIterator {
    htsFile *fp_;
    tbx_t *tbx_;

    hts_itr_t *it_;
    kstring_t str_ = {0,0,0};
    bool valid_ = false;

    TabixIterator(htsFile *fp, tbx_t *tbx, hts_itr_t *it) {
        fp_ = fp;
        tbx_ = tbx;
        it_ = it;
        assert(tbx_ && fp_ && it_);
    }

public:
    ~TabixIterator() {
        tbx_itr_destroy(it_);
        if (str_.s) {
            free(str_.s);
        }
    }

    static unique_ptr<TabixIterator> Open(htsFile *fp, tbx_t *tbx, const char* region) {
        if (!tbx) {
            return nullptr;
        }
        hts_itr_t *it = tbx_itr_querys(tbx, region);
        if (!it) {
            return nullptr;
        }
        auto ans = unique_ptr<TabixIterator>(new TabixIterator(fp, tbx, it));
        ans->Next();
        return ans;
    }

    bool Valid() const {
        return valid_ && str_.s;
    }

    const char* Line() const {
        if (Valid()) {
            assert(ks_len(&str_) == strlen(str_.s));
            return str_.s;
        }
        return nullptr;
    }

    void Next() {
        valid_ = (tbx_itr_next(fp_, tbx_, it_, &str_) >= 0);
        assert(!valid_ || str_.s);
    }
};

void TabixSlice(const std::string& spvcf_gz, std::vector<std::string> regions, std::ostream& out) {
    // Open the file
    auto fp = shared_ptr<htsFile>(hts_open(spvcf_gz.c_str(), "r"),
                                  [] (htsFile* f) { if (f && hts_close(f)) throw runtime_error("hts_close"); });
    if (!fp) {
        throw runtime_error("Failed to open " + spvcf_gz);
    }

    // Open the index file
    auto tbx = shared_ptr<tbx_t>(tbx_index_load(spvcf_gz.c_str()),
                                 [] (tbx_t *t) { if (t) tbx_destroy(t); });
    if (!tbx) {
        throw runtime_error("Falied to open .tbi/.csi index of " + spvcf_gz);
    }

    // Copy the header lines
    kstring_t str = {0,0,0};
    while (hts_getline(fp.get(), KS_SEP_LINE, &str) >= 0) {
        if (!str.l || str.s[0]!= tbx->conf.meta_char) {
            break;
        }
        out << str.s << '\n';
    }

    for (const auto& region : regions) {
        // Parse the region as either 'chrom' or 'chrom:lo-hi'
        string region_chrom, region_hi;
        uint64_t region_lo = ULLONG_MAX;
        auto c = region.find(':');
        if (c == string::npos) {
            region_chrom = region;
        } else {
            region_chrom = region.substr(0, c);
            auto d = region.find('-', c);
            if (c == 0 || d == string::npos || d <= c+1 || d >= region.size()-1) {
                throw runtime_error("invalid region " + region);
            }
            errno = 0;
            region_lo = strtoull(region.substr(c+1, d-c-1).c_str(), nullptr, 10);
            if (errno) {
                throw runtime_error("invalid region lo " + region);
            }
            region_hi = region.substr(d+1);
            assert(region_hi.size());
        }
        assert(region_chrom.size());

        // Read the first line in this region
        auto itr = TabixIterator::Open(fp.get(), tbx.get(), region.c_str());
        if (!itr || !itr->Valid()) {
            continue;
        }

        // extract INFO spVCF_checkpointPOS=ck.
        vector<char*> tokens;
        string linecpy(itr->Line());
        split(linecpy, '\t', back_inserter(tokens), 9);
        if (tokens.size() < 10) {
            throw runtime_error("read line with fewer than 10 columns");
        }
        string INFO = tokens[7];
        if (INFO.substr(0, 20) != "spVCF_checkpointPOS=") {
            // This first line happens to be a checkpoint, so we can just copy
            // all the encoded spVCF. This occurs when slicing a whole
            // chromosome since the first line for each chromosome is a
            // checkpoint, but it can also happen fortuitously in the middle.
            for (; itr->Valid(); itr->Next()) {
                out << itr->Line() << '\n';
            }
            continue;
        }
        if (region_lo == ULLONG_MAX) {
            throw runtime_error("First line for chromosome was not a checkpoint: " + region);
        }
        errno = 0;
        uint64_t ck = strtoull(INFO.substr(20, INFO.find(';')).c_str(), nullptr, 10);
        if (errno || ck >= region_lo) {
            throw runtime_error("invalid spVCF_checkpointPOS field");
        }

        // Reopen iterator on chrom:ck-hi
        assert(region_hi.size());
        string ck_region = region_chrom + ":" + to_string(ck) + "-" + region_hi;
        itr = TabixIterator::Open(fp.get(), tbx.get(), ck_region.c_str());
        if (!itr || !itr->Valid()) {
            throw runtime_error("couldn't open checkpoint region " + ck_region + " before " + region);
        }

        // Find the first checkpoint in this expanded region (it's not
        // guaranteed to be the very first result in all cases)
        uint64_t linepos = ULLONG_MAX;
        while (true) {
            linecpy = itr->Line();
            tokens.clear();
            split(linecpy, '\t', back_inserter(tokens), 9);
            if (tokens.size() < 10) {
                throw runtime_error("read line with fewer than 10 columns");
            }
            errno = 0;
            linepos = strtoull(tokens[1], nullptr, 10);
            if (errno) {
                throw runtime_error("invalid POS " + string(tokens[1]) + " while looking for checkpoint in " + ck_region);
            }
            if (string(tokens[7]).substr(0, 20) != "spVCF_checkpointPOS=") {
                break;
            }
            itr->Next();
            // We're expecting to find the checkpoint before region_lo
            if (!itr->Valid() || linepos >= region_lo) {
                throw runtime_error("couldn't find checkpoint in " + ck_region + " before " + region);
            }
        }

        // From the checkpoint, run spVCF decoder until we see a line with POS >= region_lo
        // TODO: the correct condition may be END >= region_lo. Needs careful testing.
        auto decoder = NewDecoder();
        while (true) {
            linecpy = itr->Line();

            string decoded_line = decoder->ProcessLine(&linecpy[0]);
            linecpy = decoded_line;
            tokens.clear();
            split(linecpy, '\t', back_inserter(tokens), 9);
            if (tokens.size() < 10) {
                throw runtime_error("read line with fewer than 10 columns");
            }

            errno = 0;
            linepos = strtoull(tokens[1], nullptr, 10);
            if (errno) {
                throw runtime_error("invalid POS " + string(tokens[1]) + " while looking for checkpoint in " + ck_region);
            }

            itr->Next();

            if (linepos >= region_lo) {
                // output this row as a new checkpoint
                out << decoded_line << '\n';
                break;
            }

            if (!itr->Valid()) {
                throw runtime_error("Couldn't resume from checkpoint " + ck_region + " for " + region);
            }
        }

        // Copy subsequent lines up until the next checkpoint, setting
        // spVCF_checkpointPOS=linepos.
        for (; itr->Valid(); itr->Next()) {
            linecpy = itr->Line();
            tokens.clear();
            split(linecpy, '\t', back_inserter(tokens), 9);
            if (tokens.size() < 10) {
                throw runtime_error("read line with fewer than 10 columns");
            }
            INFO = tokens[7];
            if (INFO.substr(0, 20) != "spVCF_checkpointPOS=") {
                break;
            }

            string newINFO = "spVCF_checkpointPOS=" + to_string(linepos);
            auto sc = INFO.find(';');
            if (sc != string::npos) {
                newINFO = newINFO + ";" + INFO.substr(sc+1);
            }

            for (int i = 0; i < tokens.size(); i++) {
                if (i) {
                    out << '\t';
                }
                if (i != 7) {
                    out << tokens[i];
                } else {
                    out << newINFO;
                }
            }
            out << '\n';
        }

        // Copy remaining lines unmodified
        for (; itr->Valid(); itr->Next()) {
            out << itr->Line() << '\n';
        }
    }
}

}

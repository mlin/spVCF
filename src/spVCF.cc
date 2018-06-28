#include "spVCF.h"
#include <vector>
#include <sstream>
#include <algorithm>
#include <cstdint>
#include <stdexcept>
#include <stdlib.h>
#include <assert.h>

using namespace std;

namespace spVCF {

// https://stackoverflow.com/a/236803
template<typename Out>
void split(const std::string &s, char delim, Out result) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        *(result++) = move(item);
    }
}

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
    EncoderImpl(uint64_t checkpoint_period, bool squeeze) : checkpoint_period_(checkpoint_period), squeeze_(squeeze) {}
    EncoderImpl(const EncoderImpl&) = delete;
    string ProcessLine(const string& input_line) override;

private:
    void Squeeze(vector<string>& line);

    uint64_t checkpoint_period_ = 0, since_checkpoint_ = 0;
    bool squeeze_ = false;
    vector<string> dense_entries_;
};

string EncoderImpl::ProcessLine(const string& input_line) {
    ++line_number_;
    // Pass through header lines
    if (input_line.empty() || input_line[0] == '#') {
        return input_line;
    }
    ++stats_.lines;

    // Split the tab-separated line
    vector<string> tokens;
    tokens.reserve(dense_entries_.size());
    split(input_line, '\t', back_inserter(tokens));
    if (tokens.size() < 10) {
        fail("Invalid: fewer than 10 columns");
    }
    uint64_t N = tokens.size() - 9;

    if (dense_entries_.empty()){ // First line: allocate the dense entries
        dense_entries_.resize(N);
        stats_.N = N;
    } else if (dense_entries_.size() != N) { // Subsequent line -- check expected # columns
        for (int i = 9; i < tokens.size(); i++) {
            const string& t = tokens[i];
            if (!t.empty() && t[0] == '"') {
                fail("Input seems to be sparse-encoded already");
            }
        }
        fail("Inconsistent number of samples");
    }

    if (squeeze_) {
        Squeeze(tokens);
    }

    // Pass through first nine columns
    ostringstream output_line;
    output_line << tokens[0];
    for (int i = 1; i < 9; i++) {
        output_line << '\t' << tokens[i];
    }

    uint64_t quote_run = 0; // current run-length of quotes across the row
    uint64_t sparse_cells = 0;
    // Iterate over the columns, compare each entry with the last entry
    // recorded densely.
    for (uint64_t s = 0; s < N; s++) {
        string& m = dense_entries_[s];
        const string& t = tokens[s+9];
        if (!t.empty() && t[0] == '"') {
            fail("Input seems to be sparse-encoded already");
        }
        if (m.empty() || m.size() != t.size() || m != t) {
            // Entry doesn't match the last one recorded densely for this
            // column. Output any accumulated run of quotes in the current row,
            // then this new entry, and update the state appropriately.
            if (quote_run) {
                output_line << '\t' << '"';
                if (quote_run > 1) {
                    output_line << quote_run;
                }
                quote_run = 0;
                ++sparse_cells;
            }
            output_line << '\t' << t;
            ++sparse_cells;
            m = t;
        } else {
            // Entry matches; add to the current run of quotes
            quote_run++;
        }
    }
    // Output final run of quotes
    if (quote_run) {
        output_line << '\t' << '"';
        if (quote_run > 1) {
            output_line << quote_run;
        }
        ++sparse_cells;
    }

    // CHECKPOINT -- return a densely-encoded row -- if we've hit the specified
    // period OR if we've passed half the period and this line is mostly dense
    // anyway
    if (checkpoint_period_ > 0 &&
        (since_checkpoint_ >= checkpoint_period_ ||
         (since_checkpoint_*2 >= checkpoint_period_ && sparse_cells*2 >= N))) {
        ostringstream cp;
        for (int t = 0; t < tokens.size(); t++) {
            if (t > 0) {
                cp << '\t';
            }
            cp << tokens[t];
            if (t >= 9) {
                dense_entries_[t-9] = tokens[t];
            }
            assert(tokens.size() == stats_.N+9);
        }
        since_checkpoint_ = 0;
        ++stats_.checkpoints;
        return cp.str();
    }
    ++since_checkpoint_;

    stats_.sparse_cells += sparse_cells;
    auto sparse_pct = 100*sparse_cells/N;
    if (sparse_pct <= 10) {
        ++stats_.sparse90_lines;
    }
    if (sparse_pct <= 1) {
        ++stats_.sparse99_lines;
    }

    return output_line.str();
}

// Discard QC measures (FORMAT fields) for entries which have no ALT allele
// hard-called
void EncoderImpl::Squeeze(vector<string>& line) {
    const auto& spec = line[8];
    if ((spec.size() >= 3 && spec.substr(0,3) != "GT:") && spec != "GT") {
        fail("cells don't start with genotype (GT)");
    }

    vector<string> fields, gt;
    for (int s = 9; s < line.size(); s++) {
        auto& cell = line[s];
        size_t p = 0;
        for (; p < cell.size() && cell[p] != ':'; p++) {
            switch (cell[p]) {
                case '0':
                case '.':
                case '/':
                case '|':
                    break;
                default:
                    // some ALT allele called here -- this cell is not to be
                    // touched.
                    p = cell.size();
            }
        }

        if (p < cell.size()) {
            // truncate cell to just the GT
            assert(cell[p] == ':');
            cell = cell.substr(0, p);
            ++stats_.squeezed_cells;
        }
    }
}

unique_ptr<Transcoder> NewEncoder(uint64_t checkpoint_period, bool squeeze) {
    return make_unique<EncoderImpl>(checkpoint_period, squeeze);
}

class DecoderImpl : public TranscoderBase {
public:
    DecoderImpl() = default;
    DecoderImpl(const DecoderImpl&) = delete;
    string ProcessLine(const string& input_line) override;

private:
    vector<string> dense_entries_;
};

string DecoderImpl::ProcessLine(const string& input_line) {
    ++line_number_;
    // Pass through header lines
    if (input_line.empty() || input_line[0] == '#') {
        return input_line;
    }
    ++stats_.lines;

    // Split the tab-separated line
    vector<string> tokens;
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
    ostringstream output_line;
    output_line << tokens[0];
    for (int i = 1; i < 9; i++) {
        output_line << '\t' << tokens[i];
    }

    // Iterate over the sparse columns
    uint64_t sparse_cells = (tokens.size()-9), dense_cursor = 0;
    for (uint64_t sparse_cursor = 0; sparse_cursor < sparse_cells; sparse_cursor++) {
        const string& t = tokens[sparse_cursor+9];
        if (t.empty()) {
            fail("empty cell");
        }
        if (t[0] != '"') {
            // Dense entry - remember it and copy it to the output
            // TODO: Perhaps fill QC fields with missing values (.) if they were squeezed out.
            //       The VCF spec does however say "Trailing fields can be dropped"
            if (dense_cursor >= N) {
                fail("Greater-than-expected number of columns implied by sparse encoding");
            }
            dense_entries_[dense_cursor++] = t;
            output_line << '\t' << t;
        } else {
            // Sparse entry - determine the run length
            uint64_t r = 1;
            if (t.size() > 1) {
                errno = 0;
                auto s = strtoul(t.substr(1).c_str(), nullptr, 10);
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
                output_line << '\t' << dense_entries_[dense_cursor++];
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
    if (sparse_pct <= 10) {
        ++stats_.sparse90_lines;
    }
    if (sparse_pct <= 1) {
        ++stats_.sparse99_lines;
    }

    return output_line.str();
}

unique_ptr<Transcoder> NewDecoder() {
    return make_unique<DecoderImpl>();
}

}

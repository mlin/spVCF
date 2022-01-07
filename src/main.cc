#include <iomanip>
#include <locale>
#include <fstream>
#include <getopt.h>
#include <unistd.h>
#include <assert.h>
#include <deque>
#include <thread>
#include <future>
#include <mutex>
#include "spVCF.h"

using namespace std;

enum class CodecMode {
    encode,
    squeeze_only,
    decode
};

void check_input_format(CodecMode mode, const string& first_line) {
    if (first_line.size() >= 2 && uint8_t(first_line[0]) == 0x1F && uint8_t(first_line[1]) == 0x8B) {
        throw runtime_error("input appears gzipped; decompress or pipe through `gzip -dc` for use with this tool");
    }
    const string vcf_startswith = (mode == CodecMode::decode) ? "##fileformat=spVCF" : "##fileformat=VCF";
    if (first_line.size() < vcf_startswith.size() ||
        first_line.substr(0, vcf_startswith.size()) != vcf_startswith) {
        cerr << "[WARN] input doesn't begin with "
             << vcf_startswith
             << "; this tool expects uncompressed VCF/spVCF format"
             << endl;
    }
}

// Run encoder in a multithreaded way by buffering batches of input lines and
// spawning a thread to work on each batch. Below, main_codec has a simpler
// single-threaded default way to run the codec.
spVCF::transcode_stats multithreaded_encode(CodecMode mode, uint64_t checkpoint_period, bool squeeze, double roundDP_base,
                                            size_t thread_count, istream& input_stream, ostream& output_stream) {
    assert(mode != CodecMode::decode);

    mutex mu;
    // mu protects the following two variables:
    deque<future<pair<spVCF::transcode_stats,shared_ptr<vector<string>>>>> output_batches;
    bool input_complete = false;

    // spawn "sink" task to await output blocks in order and write them to output_stream
    future<spVCF::transcode_stats> sink = async(launch::async, [&]() {
        spVCF::transcode_stats ans;
        while (true) {
            future<pair<spVCF::transcode_stats,shared_ptr<vector<string>>>> batch;
            {
                lock_guard<mutex> lock(mu);
                if (output_batches.empty()) {
                    if (input_complete) {
                        break;
                    } else {
                        using namespace chrono_literals;
                        this_thread::sleep_for(100us);
                        continue;
                    }
                }
                batch = move(output_batches.front());
                output_batches.pop_front();
            }
            auto rslt = move(batch.get());
            for (auto& line : *rslt.second) {
                output_stream << line << '\n';
                if (!output_stream.good()) {
                    throw runtime_error("I/O error");
                }
            }
            ans += rslt.first;
        }
        return ans;
    });

    // worker task to process a batch of input lines into output_batches
    auto worker = [&](shared_ptr<vector<string>> input_batch) {
        unique_ptr<spVCF::Transcoder> tc = spVCF::NewEncoder(checkpoint_period, (mode == CodecMode::encode), squeeze, roundDP_base);
        auto output_lines = make_shared<vector<string>>();
        for (auto& input_line : *input_batch) {
            output_lines->push_back(tc->ProcessLine(&input_line[0]));
        }
        return make_pair(tc->Stats(), output_lines);
    };

    // In driver thread, read batches of lines from input_stream, and spawn a worker
    // for each batch. If the number of outstanding tasks hits thread_count, wait
    // for some of them to drain.
    auto input_batch = make_shared<vector<string>>();
    size_t input_batch_size = 0;
    string input_line;
    if (getline(input_stream, input_line)) {
        check_input_format(mode, input_line);
        do {
            if (!input_line.empty() && input_line[0] != '#') {
                // TODO: it would be nice to cut off the batch at the end of each
                // chromosome, to guarantee identical checkpoint positions between
                // single- and multi-threaded encoder
                ++input_batch_size;
            }
            size_t reserve = input_line.size() * 5 / 4;
            input_batch->push_back(move(input_line));
            input_line.reserve(reserve);
            assert(input_line.empty());
            if (input_batch_size >= checkpoint_period) {
                while (true) {
                    lock_guard<mutex> lock(mu);
                    if (output_batches.size() >= thread_count) {
                        using namespace chrono_literals;
                        this_thread::sleep_for(100us);
                        continue;
                    }
                    output_batches.push_back(async(launch::async, worker, input_batch));
                    break;
                }
                input_batch = make_shared<vector<string>>();
                input_batch_size = 0;
            }
        } while (getline(input_stream, input_line));
    }
    {
        lock_guard<mutex> lock(mu);
        if (!input_batch->empty()) {
            output_batches.push_back(async(launch::async, worker, input_batch));
        }
        input_complete = true;
    }
    if (!input_stream.eof() || input_stream.bad()) {
        throw runtime_error("I/O error");
    }

    return sink.get();
}

void help_codec(CodecMode mode) {
    switch (mode) {
        case CodecMode::encode:
            cout << "spvcf encode: Encode Project VCF to Sparse Project VCF" << endl;
            cout << GIT_REVISION << "    " << __TIMESTAMP__ << endl << endl
                 << "spvcf encode [options] [in.vcf|-]" << endl
                 << "Reads VCF text from standard input if filename is empty or -" << endl << endl
                 << "Options:" << endl
                 << "  -o,--output out.spvcf  Write to out.spvcf instead of standard output" << endl
                 << "  -n,--no-squeeze        Disable lossy QC squeezing transformation (lossless run-encoding only)" << endl
                 << "  -r,--resolution        Resolution parameter r for DP rounding, rDP=floor(r^floor(log_r(DP)))" << endl
                 << "                           (default: 2.0; to increase resolution set 1.0<r<2.0)" << endl
                 << "  -p,--period P          Ensure checkpoints (full dense rows) at this period or less (default: 1000)" << endl
                 << "  -t,--threads N         Use multithreaded encoder with this number of worker threads" << endl
                 << "  -q,--quiet             Suppress statistics printed to standard error" << endl
                 << "  -h,--help              Show this help message" << endl << endl;
            break;
        case CodecMode::squeeze_only:
            cout << "spvcf squeeze: Squeeze Project VCF (without run-encoding)" << endl
                 << GIT_REVISION << "    " << __TIMESTAMP__ << endl << endl;
            cout << "spvcf squeeze [options] [in.vcf|-]" << endl
                 << "Reads VCF text from standard input if filename is empty or -" << endl << endl
                 << "Options:" << endl
                 << "  -o,--output out.vcf    Write to out.vcf instead of standard output" << endl
                 << "  -r,--resolution        Resolution parameter r for DP rounding, rDP=floor(r^floor(log_r(DP)))" << endl
                 << "                           (default: 2.0; to increase resolution set 1.0<r<2.0)" << endl
                 << "  -t,--threads N         Use multithreaded encoder with this many worker threads" << endl
                 << "  -q,--quiet             Suppress statistics printed to standard error" << endl
                 << "  -h,--help              Show this help message" << endl << endl;
            cout << "Squeezing is a lossy transformation to selectively reduce entropy in pVCF QC values." << endl
                 << "Truncate cells to GT:DP, with DP rounded down to a power of two, if: " << endl
                 << "- AD is present and indicates zero read depth for alternate alleles; OR" << endl
                 << "- VR is present and zero" << endl
                 << "May reorder fields within all cells." << endl << endl;
            break;
        case CodecMode::decode:
            cout << "spvcf decode: decode Sparse Project VCF to Project VCF" << endl;
            cout << GIT_REVISION << "    " << __TIMESTAMP__ << endl << endl
                 << "spvcf decode [options] [in.spvcf|-]" << endl
                 << "Reads spVCF text from standard input if filename is empty or -" << endl << endl
                 << "Options:" << endl
                 << "  -o,--output out.vcf  Write to out.vcf instead of standard output" << endl
                 << "  -q,--quiet           Suppress statistics printed to standard error" << endl
                 << "  -h,--help            Show this help message" << endl << endl;
            break;
    }
}

int main_codec(int argc, char *argv[], CodecMode mode) {
    bool squeeze = true;
    bool quiet = false;
    string output_filename;
    uint64_t checkpoint_period = 1000;
    size_t thread_count = 1;
    double roundDP_base = 2.0;

    static struct option long_options[] = {
        {"help", no_argument, 0, 'h'},
        {"no-squeeze", no_argument, 0, 'n'},
        {"period", required_argument, 0, 'p'},
        {"resolution", required_argument, 0, 'r'},
        {"threads", required_argument, 0, 't'},
        {"quiet", no_argument, 0, 'q'},
        {"output", required_argument, 0, 'o'},
        {0, 0, 0, 0}
    };

    int c;
    while (-1 != (c = getopt_long(argc, argv, "hnp:r:qo:t:",
                                  long_options, nullptr))) {
        switch (c) {
            case 'h':
                help_codec(mode);
                return 0;
            case 'n':
                if (mode != CodecMode::encode) {
                    help_codec(mode);
                    return -1;
                }
                squeeze = false;
                break;
            case 'p':
                if (mode == CodecMode::decode) {
                    help_codec(mode);
                    return -1;
                }
                errno=0;
                checkpoint_period = strtoull(optarg, nullptr, 10);
                if (errno) {
                    cerr << "spvcf: couldn't parse --period" << endl;
                    return -1;
                }
                break;
            case 'r':
                if (mode == CodecMode::decode) {
                    help_codec(mode);
                    return -1;
                }
                errno=0;
                roundDP_base = strtod(optarg, nullptr);
                if (errno || roundDP_base <= 1.0) {
                    cerr << "spvcf: invalid --resolution" << endl;
                    return -1;
                }
                break;
            case 't':
                if (mode == CodecMode::decode) {
                    help_codec(mode);
                    return -1;
                }
                errno=0;
                thread_count = strtoull(optarg, nullptr, 10);
                if (errno) {
                    cerr << "spvcf: couldn't parse --threads" << endl;
                    return -1;
                }
                break;
            case 'q':
                quiet = true;
                break;
            case 'o':
                output_filename = string(optarg);
                if (output_filename.empty()) {
                    help_codec(mode);
                    return -1;
                }
                break;
            default: 
                help_codec(mode);
                return -1;
        }
    }

    string input_filename;
    if (optind == argc-1) {
        input_filename = string(argv[optind]);
    } else if (optind != argc) {
        help_codec(mode);
        return -1;
    }

    // Set up input & output streams
    std::ios_base::sync_with_stdio(false);
    istream* input_stream = &cin;
    cin.tie(nullptr);
    unique_ptr<ifstream> input_box;
    if (!input_filename.empty() && input_filename != "-") {
        input_box = make_unique<ifstream>(input_filename);
        if (!input_box->good()) {
            throw runtime_error("Failed to open input file");
        }
        input_stream = input_box.get();
    } else if (isatty(STDIN_FILENO)) {
        help_codec(mode);
        return -1;
    }

    ostream* output_stream = &cout;
    unique_ptr<ofstream> output_box;
    if (!output_filename.empty()) {
        output_box = make_unique<ofstream>(output_filename);
        if (output_box->bad()) {
            throw runtime_error("Failed to open output file");
        }
        output_stream = output_box.get();
    }
    output_stream->setf(ios_base::unitbuf);

    // Encode or decode
    spVCF::transcode_stats stats;
    if (thread_count <= 1) {
        unique_ptr<spVCF::Transcoder> tc;
        if (mode == CodecMode::decode) {
            tc = spVCF::NewDecoder();
        } else {
            tc = spVCF::NewEncoder(checkpoint_period, (mode == CodecMode::encode), squeeze, roundDP_base);
        }
        string input_line;
        if (getline(*input_stream, input_line)) {
            check_input_format(mode, input_line);
            do {
                *output_stream << tc->ProcessLine(&input_line[0]);
                *output_stream << '\n';
                if (input_stream->fail() || input_stream->bad() || !output_stream->good()) {
                    throw runtime_error("I/O error");
                }
            } while (getline(*input_stream, input_line));
        }
        if (!input_stream->eof() || input_stream->bad()) {
            throw runtime_error("I/O error");
        }
        stats = tc->Stats();
    } else {
        assert(mode != CodecMode::decode);
        stats = multithreaded_encode(mode, checkpoint_period, squeeze, roundDP_base, thread_count,
                                     *input_stream, *output_stream);
    }

    // Close up
    if (output_box) {
        output_box->close();
        if (output_box->fail()) {
            throw runtime_error("Failed to close output file");
        }
    }

    // Output stats
    if (!quiet) {
        cerr.imbue(locale(""));
        cerr << "N = " << fixed << stats.N << endl;
        cerr << "dense cells = " << fixed << stats.N*stats.lines << endl;
        if (squeeze) {
            cerr << "squeezed cells = " << fixed << stats.squeezed_cells << endl;
        }
        if (mode != CodecMode::squeeze_only) {
            cerr << "sparse cells = " << fixed << stats.sparse_cells << endl;
            cerr << "lines (non-header) = " << fixed << stats.lines << endl;
            cerr << "lines (75% sparse) = " << fixed << stats.sparse75_lines << endl;
            cerr << "lines (90% sparse) = " << fixed << stats.sparse90_lines << endl;
            cerr << "lines (99% sparse) = " << fixed << stats.sparse99_lines << endl;
        }
        if (mode == CodecMode::encode) {
            cerr << "checkpoints = " << fixed << stats.checkpoints << endl;
        }
    }

    return 0;
}

void help_tabix() {
    cout << "spvcf tabix: use a .tbi index to slice a spVCF bgzip file by genomic range" << endl;
    cout << GIT_REVISION << "    " << __TIMESTAMP__ << endl << endl
         << "spvcf tabix [options] in.spvcf.gz chr1:1000-2000 [chr2 ...]" << endl
         << "Requires tabix index present e.g. in.spvcf.gz.tbi. Includes all header lines." << endl << endl
         << "Options:" << endl
         << "  -o,--output out.spvcf  Write to out.spvcf instead of standard output" << endl
         << "  -h,--help              Show this help message" << endl << endl;
}

int main_tabix(int argc, char *argv[]) {
    string output_filename;

    static struct option long_options[] = {
        {"help", no_argument, 0, 'h'},
        {"output", required_argument, 0, 'o'},
        {0, 0, 0, 0}
    };

    int c;
    while (-1 != (c = getopt_long(argc, argv, "ho:",
                                  long_options, nullptr))) {
        switch (c) {
            case 'h':
                help_tabix();
                return 0;
            case 'o':
                output_filename = string(optarg);
                if (output_filename.empty()) {
                    help_tabix();
                    return -1;
                }
                break;
            default:
                help_tabix();
                return -1;
        }
    }

    if (optind+2 > argc) {
        help_tabix();
        return -1;
    }

    string input_filename = argv[optind++];
    vector<string> regions;
    while (optind < argc) {
        regions.push_back(argv[optind++]);
    }

    ostream* output_stream = &cout;
    unique_ptr<ofstream> output_box;
    if (!output_filename.empty()) {
        output_box = make_unique<ofstream>(output_filename);
        if (output_box->bad()) {
            throw runtime_error("Failed to open output file");
        }
        output_stream = output_box.get();
    }

    spVCF::TabixSlice(input_filename, regions, *output_stream);
    return 0;
}

void help() {
    cout << "spvcf: Sparse Project VCF tool" << endl;
    cout << GIT_REVISION << "    " << __TIMESTAMP__ << endl << endl
         << "subcommands:" << endl
         << "  encode   encode Project VCF to spVCF" << endl
         << "  squeeze  squeeze Project VCF" << endl
         << "  decode   decode spVCF to Project VCF" << endl
         << "  tabix    use a .tbi index to slice a spVCF bgzip file by genomic range" << endl
         << "  help     show this help message" << endl << endl;
}

int main(int argc, char *argv[]) {
    if (argc <= 1) {
        help();
        return -1;
    }

    string subcommand = argv[1];
    if (subcommand == "help" || subcommand == "-h" || subcommand == "--help") {
        help();
        return 0;
    }

    optind = 2;
    if (subcommand == "encode") {
        return main_codec(argc, argv, CodecMode::encode);
    } else if (subcommand == "squeeze") {
        return main_codec(argc, argv, CodecMode::squeeze_only);
    } else if (subcommand == "decode") {
        return main_codec(argc, argv, CodecMode::decode);
    } else if (subcommand == "tabix") {
        return main_tabix(argc, argv);
    }

    help();
    return -1;
}

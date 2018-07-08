#include <iostream>
#include <iomanip>
#include <locale>
#include <fstream>
#include <getopt.h>
#include <unistd.h>
#include "spVCF.h"

using namespace std;

void usage() {
    cout << "spvcf: Encode Project VCF to Sparse Project VCF or vice-versa. " << endl;
    cout << "       " << GIT_REVISION << "    " << __TIMESTAMP__ << endl << endl
         << "spvcf [options] [in.vcf|-]" << endl
         << "Reads VCF text from standard input if filename is empty or -" << endl << endl
         << "Options:" << endl
         << "  -o,--output out.spvcf  Write to out.spvcf instead of standard output" << endl
         << "  -d,--decode            Decode from the sparse format instead of encoding to" << endl
         << "  -p,--period P          Ensure checkpoints (full dense rows) at this period or less (default: 1000)" << endl
         << "  -q,--quiet             Suppress statistics printed to standard error" << endl
         << "  -h,--help              Show this usage message" << endl << endl
         << "Lossy transformation to increase compression: " << endl
         << "  -S,--squeeze           Truncate cells to GT:DP, with DP rounded down to a power of two, if: " << endl
         << "                         - AD is present and indicates zero read depth for alternate alleles; OR" << endl
         << "                         - VR is present and zero" << endl
         << "                         Reorders fields within all cells."<< endl << endl;
}

int main(int argc, char *argv[]) {
    // Parse arguments
    bool squeeze = false;
    bool decode = false;
    bool quiet = false;
    string output_filename;
    uint64_t checkpoint_period = 1000;

    static struct option long_options[] = {
        {"help", no_argument, 0, 'h'},
        {"decode", no_argument, 0, 'd'},
        {"squeeze", no_argument, 0, 'S'},
        {"period", required_argument, 0, 'p'},
        {"quiet", no_argument, 0, 'q'},
        {"output", required_argument, 0, 'o'},
        {0, 0, 0, 0}
    };

    int c;
    while (-1 != (c = getopt_long(argc, argv, "hSdp:qo:",
                                  long_options, nullptr))) {
        switch (c) {
            case 'h':
                usage();
                return 0;
            case 'S':
                squeeze = true;
                break;
            case 'd':
                decode = true;
                break;
            case 'p':
                errno=0;
                checkpoint_period = strtoull(optarg, nullptr, 10);
                if (errno) {
                    cerr << "spvcf: couldn't parse --period" << endl;
                    return -1;
                }
                break;
            case 'q':
                quiet = true;
                break;
            case 'o':
                output_filename = string(optarg);
                if (output_filename.empty()) {
                    usage();
                    return -1;
                }
                break;
            default: 
                usage();
                return -1;
        }
    }

    if (squeeze && decode) {
        cerr << "spvcf: --squeeze and --decode should not be specified together" << endl;
        return -1;
    }

    string input_filename;
    if (optind == argc-1) {
        input_filename = string(argv[optind]);
    } else if (optind != argc) {
        usage();
        return -1;
    }

    // Set up input & output streams
    istream* input_stream = &cin;
    cin.tie(nullptr);
    unique_ptr<ifstream> input_box;
    if (!input_filename.empty()) {
        input_box = make_unique<ifstream>(input_filename);
        if (!input_box->good()) {
            throw runtime_error("Failed to open input file");
        }
        input_stream = input_box.get();
    } else if (isatty(STDIN_FILENO)) {
        usage();
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

    // Encode or decode
    unique_ptr<spVCF::Transcoder> tc;
    if (decode) {
        tc = spVCF::NewDecoder();
    } else {
        tc = spVCF::NewEncoder(checkpoint_period, squeeze);
    }
    string input_line;
    for (; getline(*input_stream, input_line); ) {
        *output_stream << tc->ProcessLine(&input_line[0]) << '\n';
        if (input_stream->fail() || input_stream->bad() || !output_stream->good()) {
            throw runtime_error("I/O error");
        }
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
        auto stats = tc->Stats();
        cerr.imbue(locale(""));
        cerr << "N = " << fixed << stats.N << endl;
        cerr << "dense cells = " << fixed << stats.N*stats.lines << endl;
        if (squeeze) {
            cerr << "squeezed cells = " << fixed << stats.squeezed_cells << endl;
        }
        cerr << "sparse cells = " << fixed << stats.sparse_cells << endl;
        cerr << "lines (non-header) = " << fixed << stats.lines << endl;
        cerr << "lines (75% sparse) = " << fixed << stats.sparse75_lines << endl;
        cerr << "lines (90% sparse) = " << fixed << stats.sparse90_lines << endl;
        cerr << "lines (99% sparse) = " << fixed << stats.sparse99_lines << endl;
        if (!decode) {
            cerr << "checkpoints = " << fixed << stats.checkpoints << endl;
        }
    }

    return 0;
}

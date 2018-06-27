#include <iostream>
#include <iomanip>
#include <locale>
#include <fstream>
#include <getopt.h>
#include <unistd.h>
#include "spVCF.h"

using namespace std;

void usage() {
    cout << "spvcf [options] [in.vcf|-]" << endl;
    cout << "Read from standard input if input filename is empty or -" << endl;
    cout << "Options:" << endl;
    cout << "  -o,--output out.spvcf  Write to out.spvcf instead of standard output" << endl;
    cout << "  -S,--squeeze           Discard QC measures from cells with no ALT allele called" << endl;
    cout << "  -d,--decode            Decode from the sparse format instead of encoding to" << endl;
    cout << "  -q,--quiet             Suppress statistics printed to standard error" << endl;
    cout << "  -h,--help              Show this usage message" << endl;
    cout << "source revision: " << GIT_REVISION << endl;
}

int main(int argc, char *argv[]) {
    // Parse arguments
    bool squeeze = false;
    bool decode = false;
    bool quiet = false;
    string output_filename;

    static struct option long_options[] = {
        {"help", no_argument, 0, 'h'},
        {"squeeze", no_argument, 0, 'S'},
        {"decode", no_argument, 0, 'd'},
        {"quiet", no_argument, 0, 'q'},
        {"output", required_argument, 0, 'o'},
        {0, 0, 0, 0}
    };

    int c;
    while (-1 != (c = getopt_long(argc, argv, "hSdqo:",
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
    unique_ptr<spVCF::Transcoder> tc = decode ? spVCF::NewDecoder() : spVCF::NewEncoder(squeeze);
    for (string input_line; getline(*input_stream, input_line); ) {
        *output_stream << tc->ProcessLine(input_line) << endl;
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
        cerr << "lines (90% sparse) = " << fixed << stats.sparse90_lines << endl;
        cerr << "lines (99% sparse) = " << fixed << stats.sparse99_lines << endl;
    }

    return 0;
}

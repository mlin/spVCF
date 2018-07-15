#include <iomanip>
#include <locale>
#include <fstream>
#include <getopt.h>
#include <unistd.h>
#include "spVCF.h"

using namespace std;

void help_codec(bool decode) {
    if (!decode) {
        cout << "spvcf encode: Encode Project VCF to Sparse Project VCF" << endl;
        cout << GIT_REVISION << "    " << __TIMESTAMP__ << endl << endl
            << "spvcf encode [options] [in.vcf|-]" << endl
            << "Reads VCF text from standard input if filename is empty or -" << endl << endl
            << "Options:" << endl
            << "  -o,--output out.spvcf  Write to out.spvcf instead of standard output" << endl
            << "  -p,--period P          Ensure checkpoints (full dense rows) at this period or less (default: 1000)" << endl
            << "  -q,--quiet             Suppress statistics printed to standard error" << endl
            << "  -h,--help              Show this help message" << endl << endl
            << "Lossy transformation to increase compression: " << endl
            << "  -S,--squeeze           Truncate cells to GT:DP, with DP rounded down to a power of two, if: " << endl
            << "                         - AD is present and indicates zero read depth for alternate alleles; OR" << endl
            << "                         - VR is present and zero" << endl
            << "                         Reorders fields within all cells."<< endl << endl;
    } else {
        cout << "spvcf decode: decode Sparse Project VCF to Project VCF" << endl;
        cout << GIT_REVISION << "    " << __TIMESTAMP__ << endl << endl
             << "spvcf decode [options] [in.spvcf|-]" << endl
             << "Reads spVCF text from standard input if filename is empty or -" << endl << endl
             << "Options:" << endl
             << "  -o,--output out.vcf  Write to out.vcf instead of standard output" << endl
             << "  -q,--quiet           Suppress statistics printed to standard error" << endl
             << "  -h,--help            Show this help message" << endl << endl;
    }
}

int main_codec(int argc, char *argv[], bool decode) {
    bool squeeze = false;
    bool quiet = false;
    string output_filename;
    uint64_t checkpoint_period = 1000;

    static struct option long_options[] = {
        {"help", no_argument, 0, 'h'},
        {"squeeze", no_argument, 0, 'S'},
        {"period", required_argument, 0, 'p'},
        {"quiet", no_argument, 0, 'q'},
        {"output", required_argument, 0, 'o'},
        {0, 0, 0, 0}
    };

    int c;
    while (-1 != (c = getopt_long(argc, argv, "hSp:qo:",
                                  long_options, nullptr))) {
        switch (c) {
            case 'h':
                help_codec(decode);
                return 0;
            case 'S':
                if (decode) {
                    help_codec(decode);
                    return -1;
                }
                squeeze = true;
                break;
            case 'p':
                if (decode) {
                    help_codec(decode);
                    return -1;
                }
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
                    help_codec(decode);
                    return -1;
                }
                break;
            default: 
                help_codec(decode);
                return -1;
        }
    }

    string input_filename;
    if (optind == argc-1) {
        input_filename = string(argv[optind]);
    } else if (optind != argc) {
        help_codec(decode);
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
        help_codec(decode);
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
    unique_ptr<spVCF::Transcoder> tc;
    if (decode) {
        tc = spVCF::NewDecoder();
    } else {
        tc = spVCF::NewEncoder(checkpoint_period, squeeze);
    }
    string input_line;
    for (; getline(*input_stream, input_line); ) {
        *output_stream << tc->ProcessLine(&input_line[0]);
        *output_stream << '\n';
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
         << "  encode  encode Project VCF to spVCF" << endl
         << "  decode  decode spVCF to Project VCF" << endl
         << "  tabix   use a .tbi index to slice a spVCF bgzip file by genomic range" << endl
         << "  help    show this help message" << endl << endl;
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
        return main_codec(argc, argv, false);
    } else if (subcommand == "decode") {
        return main_codec(argc, argv, true);
    } else if (subcommand == "tabix") {
        return main_tabix(argc, argv);
    }

    help();
    return -1;
}

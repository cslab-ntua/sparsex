#include "bench.h"
#include <cstring>
#include <getopt.h>

static const char *program_name;

static struct option long_options[] = {
    {"directory",   required_argument,  0, 'd'},
    {"file",        required_argument,  0, 'f'},
    {"library",     required_argument,  0, 'l'},
    {"help",        no_argument,        0, 'h'}
};

void PrintUsage(std::ostream &os)
{
    os << "Usage:  " << program_name 
       << " -f <mmf_file> [-l library]\n"
       << "\t" << program_name
       << " -d <directory> [-l library]\n"
       << "\t-f    Run SpMV kernel on file.\n"
       << "\t-d    Run SpMV kernel on all files in a directory.\n"
       << "\t-l    Use one of the available libraries (MKL, pOSKI, LIBCSX).\n"
       << "\t-h    Print this help message and exit.\n";
}

int main(int argc, char **argv)
{
    char c;
    int option_index = 0;
    char *directory = 0;
    char *file = 0;
    char *library = 0;

    program_name = argv[0];
    while ((c = getopt_long(argc, argv, "d:f:l:h", long_options,
                            &option_index)) != -1) {
        switch (c) {
        case 'd':
            directory = optarg;
            break;
        case 'f':
            file = optarg;
            break;
        case 'l':
            library = optarg;
            break;
        case 'h':
            PrintUsage(std::cerr);
            exit(0);
        default:
            PrintUsage(std::cerr);
            exit(1);
        }
    }

    if (optind < argc || argc < 2) {
        PrintUsage(std::cerr);
        exit(1);
    }

    if (directory)
        Bench_Directory(directory, library, NULL, 128, 5);
    else if (file)
        Bench_Matrix(file, library, NULL, 128, 5);

    return 0;
}

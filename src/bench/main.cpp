/*
 * Copyright (C) 2014, Computing Systems Laboratory (CSLab), NTUA.
 * Copyright (C) 2014, Athena Elafrou
 * All rights reserved.
 *
 * This file is distributed under the BSD License. See LICENSE.txt for details.
 */

/**
 * \file main.cpp
 * \brief Benchmarking utility
 *
 * \author Computing Systems Laboratory (CSLab), NTUA
 * \date 2011&ndash;2014
 * \copyright This file is distributed under the BSD License. See LICENSE.txt
 * for details.
 */

#include "Bench.hpp"
#include <cstring>
#include <getopt.h>

static const char *program_name;

extern unsigned int OUTER_LOOPS;
extern unsigned long LOOPS;
extern unsigned int NR_THREADS;

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

static unsigned int GetOptionOuterLoops()
{
    const char *loops_env = getenv("OUTER_LOOPS");
    unsigned int ret = 5;
    
    if (loops_env) {
        ret = atoi(loops_env);
        if (ret < 0)
            ret = 0;
    }
    
    return ret;
}

static unsigned long GetOptionIterations()
{
    const char *loops_env = getenv("LOOPS");
    unsigned long ret = 128;
    
    if (loops_env) {
        ret = atoi(loops_env);
        if (ret < 0)
            ret = 0;
    }
    
    return ret;
}

static unsigned int GetOptionNrThreads()
{
    const char *threads_env = getenv("NR_THREADS");
    unsigned int ret = 1;
    
    if (threads_env) {
        ret = atoi(threads_env);
        if (ret < 0)
            ret = 0;
    }
    
    return ret;
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

    NR_THREADS = GetOptionNrThreads();
    OUTER_LOOPS = GetOptionOuterLoops();
    LOOPS = GetOptionIterations();

    if (directory)
        Bench_Directory(directory, library, NULL);
    else if (file)
        Bench_Matrix(file, library, NULL);

    return 0;
}

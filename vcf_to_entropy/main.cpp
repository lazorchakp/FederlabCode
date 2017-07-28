/*
    Peter Lazorchak
    lazorchakp@jhu.edu
    5/25/2017
*/

#include <fstream>
#include <iostream>

#include <unistd.h>
#include <string>

// #include <getopt.h>

#include "functions.h"

int main(int argc, char *argv[]) {
    if (argc < 2) {
        usage(argv[0]);
        return 1;
    }

    bool gl = true;
    std::string outFileNameStr("undefined");
    char opt = 0;

    while ((opt = getopt(argc, argv, "o:c")) != -1) {
        switch (opt) {
        case 'o':
            outFileNameStr = optarg;
            break;
        case 'c':
            gl = false;
            break;
        case '?':
        default:
            usage(argv[0]);
            return 1;
        }
    }

    std::ifstream infile;
    infile.open(argv[argc - 1]);
    // check if the file opened properly
    if (!infile) {
        std::cerr << "Error: unable to open " << argv[argc - 1] << std::endl;
        return 1;
    }

    // check if the file has the correct vcf header
    if (!checkHeader(infile)) {
        std::cerr << "Error: missing vcf header in " << argv[argc - 1]
            << std::endl;
        infile.close();
        return 1;
    }

    // name the output file if its name was not specified by command line opt
    std::string inFileNameStr(argv[argc - 1]);
    if (outFileNameStr.compare("undefined") == 0) {
        // find_last_of removes the file extension
        outFileNameStr = inFileNameStr.substr(0, inFileNameStr.find_last_of('.'))
            + ".entropy.dat";
    }

    std::cout << "Converting " << argv[argc - 1] << "..." << std::endl;

    // generate a new file with the entropy program's expected format
    if (!convert(infile, gl, (char *)outFileNameStr.c_str())) {
        std::cerr << "Error: failed to convert " << argv[argc - 1]
            << std::endl;
        infile.close();
        return 1;
    }

    std::cout << "Conversion successful" << std::endl;
    infile.close();
}

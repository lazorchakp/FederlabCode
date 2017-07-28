/*
    Peter Lazorchak
    lazorchakp@jhu.edu
    6/5/2017
*/

#include <iostream>
#include <fstream>
#include <map>
#include <utility>
#include <string>

#include <unistd.h>

#include "func.h"

using std::cerr;
using std::endl;
using std::map;
using std::pair;
using std::string;

int main(int argc, char *argv[]) {
    if (argc < 2) {
        usage(argv[0]);
        return 1;
    }

    dataContainer data;
    char opt = 0;

    while ((opt = getopt(argc, argv, "i:g:o:a:")) != -1) {
        switch (opt) {
        case 'i':
            data.inFileName = optarg;
            break;
        case 'g':
            data.groupsFileName = optarg;
            break;
        case 'o':
            data.outFileName = optarg;
            break;
        case 'a':
            data.analyticsName = optarg;
            break;
        case '?':
        default:
            usage(argv[0]);
            return 1;
        }
    }

    if (data.inFileName.compare("undefined") == 0) {
        cerr << "Error: name of vcf file not specified (use -i)" << endl;
        return 1;
    }

    if (data.groupsFileName.compare("undefined") == 0) {
        cerr <<  "Error: name of groups file not specified (use -g)" << endl;
        return 1;
    }

    // name the output file if its name was not specified by command line opt
    if (data.outFileName.compare("undefined") == 0) {
        // find_last_of removes the file extension
        data.outFileName = data.inFileName.substr(0,
            data.inFileName.find_last_of('.')) + ".ld_filtered.vcf";
    }

    // load and filter the data
    cerr << "Loading snp groups" << endl;
    if (!loadSnpGroups(&data)) {
        return 1;
    }
    cerr << "Filtering snps" << endl;
    if (!filter(&data)) {
        return 1;
    }
    cerr << "Filtration successful" << endl;

    return 0;
}

/*
    Peter Lazorchak
    lazorchakp@jhu.edu
    5/25/2017
*/

#include <fstream>
#include <string>
#include <sstream>
#include <stdlib.h>
#include <iostream>

#include "functions.h"

using std::string;

// see functions.h for function descriptions

void usage(char *cmd) {
    std::cerr << "\nUsage: " << cmd << " [OPTIONS] <infile.vcf>\n"
        << std::endl;
    std::cerr <<
        "Summary of options:\n\n"
        "-o <outfile_name> [default: <infile>.entropy.dat]\n"
        "\tThe name of the output file\n\n"
        "-c [default: genotype likelihood format]\n"
        "\tOutput file will be in count format\n"
    << std::endl;
}

bool checkHeader(std::ifstream &infile) {
    string line;
    if (getline(infile, line)) {
        string cmp("##fileformat=VCFv"); // excludes version number
        line = line.substr(0, cmp.length());
        // if the beginning of the first line exactly matches cmp
        if (line.compare(cmp) == 0) {
            return true;
        }
    }
    return false;
}

bool convert(std::ifstream &infile, bool gl, char *outFileName) {
    std::ofstream outfile;
    outfile.open(outFileName);
    // check if there was a problem opening the output file
    if (!outfile) {
        std::cout << "Failed to open " << outFileName << std::endl;
        return false;
    }
    // ignore every line beginning with "##"
    string line("##");
    int dataBegin = 0;
    while (line.substr(0, 2).compare("##") == 0) {
        dataBegin = infile.tellg();
        getline(infile, line);
    }
    // dataBegin is now at the beginning of the line with a single '#'

    if (!gl) {
        return convertCount(infile, outfile);
    } else {
        // need to back up a line before calling getIndividuals
        infile.seekg(dataBegin);
        outfile << getIndividuals(infile, false, outfile) << ' ';
        outfile << getLoci(infile) << std::endl;
        // avoid re-reading the beginning lines of the file
        infile.clear();
        infile.seekg(dataBegin);
        getIndividuals(infile, true, outfile);
        outfile << std::endl;
        return convertGl(infile, outfile);
    }
}

bool convertCount(std::ifstream &infile, std::ofstream &outfile) {
    // ignore every line beginning with '#'
    string line("#");
    while (line[0] == '#') {
        getline(infile, line);
    }

    while (!line.empty()) {
        std::stringstream stream(line);
        string word;
        outfile << "locus ";
        stream >> word;
        outfile << word << ' '; // CHROM
        stream >> word;
        outfile << word << ' '; // POS
        // the "error" component would belong here
        outfile << std::endl;

        // skip ID, REF, ALT, QUAL, FILTER, INFO, and FORMAT
        for (int i = 0; i < 7; ++i) {
            stream >> word;
        }

        // run once for each individual
        while (stream >> word) {
            // format of word is #/#:<count1>,<count2>:# etc
            // TODO regex would be better here
            // finds the first colon index, then adds one to get to count1
            size_t firstNum = word.find(':') + 1;
            if (word[firstNum] == '.') {
                outfile << ". ." << std::endl;
                continue;
            }
            size_t secondColon = word.find(':', firstNum);
            // finds the comma index, then adds one to get to count2
            size_t secondNum = word.find(',', firstNum) + 1;
            // converts substrings of count1 and count2 to longs, and appends
            // data containing only '.' will be left as '.'
            outfile << word.substr(firstNum, secondNum - 1 - firstNum);
            outfile << ' ';
            // old version: strtol(word.substr(...).c_str(), NULL, 10)
            // the above line would convert '.' to 0
            outfile << word.substr(secondNum, secondColon - secondNum);
            outfile << std::endl;
        }

        getline(infile, line);
    }

    outfile.close();
    return true;
}

int getIndividuals(std::ifstream &infile, bool write, std::ofstream &outfile) {
    // shouldn't be necessary, but just in case
    string line("##");
    while (line.substr(0, 2).compare("##") == 0) {
        getline(infile, line);
    }
    // line must be "#CHROM POS ... FORMAT <ind1> <ind2> ..."
    // skip over column headers (there are 9 of these)
    std::stringstream stream(line);
    string word;
    for (int i = 0; i < 9; ++i) {
        stream >> word;
    }
    // the next read will be the first individual
    // count the number of reads made on the rest of the line
    int nind = 0;
    if (write) {
        if (stream >> word) {
            outfile << word; // prevents there from being an extra space
            ++nind;
        }
        while (stream >> word) {
            outfile << ' ' << word;
            ++nind;
        }
    } else {
        while (stream >> word) {
            ++nind;
        }
    }
    return nind;
}

int getLoci(std::ifstream &infile) {
    // ignore every line beginning with '#'
    string line("#");
    while (line[0] == '#') {
        getline(infile, line);
    }

    // start at 1 because the first line is currently stored
    int nloci = 1;
    while (getline(infile, line)) {
        nloci++;
    }
    return nloci;
}

bool convertGl(std::ifstream &infile, std::ofstream &outfile) {
    // ignore every line beginning with '#'
    string line("#");
    while (line[0] == '#') {
        getline(infile, line);
    }

    while (!line.empty()) {
        std::stringstream stream(line);
        string word;
        stream >> word;
        outfile << word; // add the CHROM data to outfile
        stream >> word;
        outfile << ':' << word << ' '; // add the POS data to the outfile

        // skip ID, REF, ALT, QUAL, FILTER, INFO, and FORMAT
        for (int i = 0; i < 7; ++i) {
            stream >> word;
        }

        while (stream >> word) {
            size_t firstNum = word.find_last_of(':') + 1;
            if (word[firstNum] == '.') {
                // data must be missing
                outfile << " . . .";
            } else {
                size_t secondNum = word.find(',', firstNum) + 1;
                size_t thirdNum = word.find(',', secondNum) + 1;
                outfile << ' ';
                outfile << word.substr(firstNum, secondNum - 1 - firstNum);
                outfile << ' ';
                outfile << word.substr(secondNum, thirdNum - 1 - secondNum);
                outfile << ' ';
                outfile << word.substr(thirdNum);
            }
        }
        outfile << std::endl;
        getline(infile, line);
    }
    outfile.close();
    return true;
}

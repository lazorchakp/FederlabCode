/*
    Peter Lazorchak
    lazorchakp@jhu.edu
    6/5/2017
*/

// see func.h for function descriptions

#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <utility>
#include <string>
#include <vector>
#include <algorithm>

#include "func.h"

using std::string;
using std::map;

std::ofstream analytics;

void usage(char *cmd) {
    std::cerr << "\nUsage: " << cmd
        << " -i <infile.vcf> -g <groupsfile> [OPTIONS]\n" << std::endl;
    std::cerr <<
"Summary of options:\n\n"
"-i <infile_name> [REQUIRED]\n"
"\tThe name of the .vcf file to be filtered (this file will not be modified)\n\n"
"-g <groupsfile_name> [REQUIRED]\n"
"\tThe name of the file containing LD group info, each line formatted as\n"
"\t<groupID> <CHROM>_<POS>\n\n"
"\tnote: <CHROM>_<POS> must exactly match CHROM and POS in the input .vcf file\n\n"
"-o <outfile_name> [default: <infile>.ld_filtered.vcf]\n"
"\tThe name of the filterd output file to be generated\n\n"
"-a <analytics_name> [default: do not output analytics]\n"
"\tOutput snp analytics as a csv file with the specified name\n"
    << std::endl;
}

bool loadSnpGroups(dataContainer *data) {
    std::ifstream groupsFile(data->groupsFileName);
    if (!groupsFile) {
        std::cerr << "Error: unable to open " << data->groupsFileName
            << std::endl;
        return false;
    }

    string line;
    std::map<string, int> groupIDs; // contains groupIDs and their int codes
    string group;
    string snp;
    int groupCode = 0; // the first group will have int value of 0
    while(getline(groupsFile, line)) {
        std::stringstream stream(line);
        stream >> group;
        if (groupIDs.count(group) == 0) {
            groupIDs[group] = groupCode;
            groupCode++;
        }
        stream >> snp;
        data->snpGroups[snp] = groupIDs[group];
    }

    groupsFile.close();
    return true;
}

bool filter(dataContainer *data) {
    // open input and output files and ensure there were no problems
    std::ifstream infile(data->inFileName);
    if (!infile) {
        std::cerr << "Error: unable to open " << data->inFileName << std::endl;
        return false;
    }
    std::ofstream outfile(data->outFileName);
    if (!outfile) {
        std::cerr << "Error: unable to open " << data->outFileName
            << std::endl;
        infile.close();
        return false;
    }
    if (data->analyticsName.compare("undefined") != 0) {
        analytics.open(data->analyticsName);
        if (!analytics) {
            std::cerr << "Error: unable to open " << data->analyticsName
                << std::endl;
            infile.close();
            outfile.close();
            return false;
        }
        analytics << "snp,group,median,missing,score" << std::endl;
    }

    string line;
    getline(infile, line);
    // copy every line starting with '#' in the input file to the output file
    while (line[0] == '#') {
        outfile << line << std::endl;
        getline(infile, line);
    }
    // line currently stores the first line of real data
    string snp;
    string pos;
    while (!line.empty()) {
        // get the snp in <CHROM>_<POS> format
        std::stringstream stream(line);
        stream >> snp;
        snp.push_back('_');
        stream >> pos;
        snp += pos;
        // if the snp is in an LD group, decide if it is the best candidate
        if (data->snpGroups.count(snp) == 1) {
            updateBestGroupScore(data, snp, stream, line);
        } else {
            outfile << line << std::endl;
        }
        getline(infile, line);
    }

    writeBestSnps(data, outfile);

    infile.close();
    outfile.close();
    if (analytics.is_open()) {
        analytics.close();
    }
    return true;
}

void updateBestGroupScore(dataContainer *data, string& snp,
        std::stringstream& stream, string& line) {
    string word;
    int missing = 0;
    double median = 0.0;
    std::map<int, int> hist; // a histogram of for incoming read values
    int validInds = 0;
    // stream has already passed over CHROM and POS
    // skip over ID, REF, ALT, QUAL, FILTER, INFO, & FORMAT
    for (int i = 0; i < 7; i++) {
        stream >> word;
    }
    while (stream >> word) {
        // handle missing data
        if (word[0] == '.') {
            ++missing;
        } else {
            // format is n/n:n,n:N:n:n,n,n where n=number and N=desired number
            int colonBefore = word.find_first_of(':', word.find(':') + 1);
            int colonAfter = word.find_first_of(':', colonBefore);
            int currReads = atoi(word.substr(colonBefore + 1,
                colonAfter - colonBefore - 1).c_str());
            ++hist[currReads];
            ++validInds;
        }
    }
    // find the median
    median = findMedian(hist, validInds);
    
    // conditionally add the snp to bestGroupScore
    double score = median / missing;
    int groupCode = data->snpGroups[snp];
    std::pair<double, string> *best = &(data->bestGroupScore[groupCode]);
    if (best->first < score) {
        best->first = score;
        best->second = line;
    }
    if (analytics.is_open()) {
        analytics << snp << ',' << groupCode << ',' << median << ',' << missing
            << ',' << score << std::endl;
    }
}

double findMedian(map<int, int>& hist, int size) {
    double median = 0.0;
    auto forward = hist.begin();
    int count = size / 2;
    // this loop finds the value immediately before the median for odd data
    // sets, or the first of the two values to be averaged for even sets
    while (count > 0) {
        count -= forward->second;
        median = forward->first;
        ++forward;
    }
    // if the next value is different, update the median
    // if count is less than 0, we must have been in the middle of a bin anyway
    if (count == 0) {
        if (size % 2 == 0) {
            median = (median + forward->first) / 2.0;
        } else {
            median = forward->first;
        }
    }
    return median;
}

void writeBestSnps(dataContainer *data, std::ofstream& outfile) {
    for (auto itr = data->bestGroupScore.begin();
            itr != data->bestGroupScore.end(); ++itr) {
        outfile << itr->second.second << std::endl;
    }
}

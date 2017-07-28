/*
    Peter Lazorchak
    lazorchakp@jhu.edu
    6/5/2017
*/

#ifndef FUNC_H
#define FUNC_H

#include <map>
#include <utility>
#include <string>
#include <fstream>

using std::string;
using std::pair;

/*
    Stores the essential data and variables to simplify passing them between
    functions.
*/
struct dataContainer {
    // string is CHROM_POS, int is LD group
    std::map<string, int> snpGroups;

    // first int is LD group; pair contains the highest score in the group and
    // the full line in the vcf file associated with this score
    std::map<int, pair<double, string>>  bestGroupScore;

    string inFileName;
    string groupsFileName;
    string outFileName;
    string analyticsName;
    bool analyze;

    dataContainer() : inFileName("undefined"), groupsFileName("undefined"),
        outFileName("undefined"), analyticsName("undefined") {}
};

/*
    Displays the correct usage of the program.
*/
void usage(char *cmd);

/*
    Loads snps and their corresponding group into the snpGroups map.
    Returns true on success.
*/
bool loadSnpGroups(dataContainer *data);

/*
    Wrapper for filtering the data. Ultimately writes out a filtered vcf file.
    Returns true on success.
*/
bool filter(dataContainer *data);

/*
    Scores a snp for data completeness. If the snp has a higher score than the
    current high score for its group, it becomes the new best snp in the group.
*/
void updateBestGroupScore(dataContainer *data, string& snp,
        std::stringstream& stream, string& line);

/*
    Finds the median of histogram hist on a set of data with size elements.
*/
double findMedian(std::map<int, int>& hist, int size);

/*
    Writes the best snps for each LD group to the bottom of outfile in standard
    vcf line format.
*/
void writeBestSnps(dataContainer *data, std::ofstream& outfile);

#endif
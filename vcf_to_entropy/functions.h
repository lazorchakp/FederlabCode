/*
    Peter Lazorchak
    lazorchakp@jhu.edu
    5/25/2017
*/

#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <fstream>

/*
    Displays the correct usage of the program.
*/
void usage(char *cmd);

/*
    Checks if a file has the proper vcf header. Returns true for valid header.
*/
bool checkHeader(std::ifstream &infile);

/*
    Converts data from the vcf format to the format expected by the entropy
    program. Writes to a file with the specified name. Returns true on success.
    This calls either convertCount() or convertGl().
*/
bool convert(std::ifstream &infile, bool gl, char *outFileName);

/*
    Converts vcf data to data in the count format.
*/
bool convertCount(std::ifstream &infile, std::ofstream &outfile);

/*
    Gets the number of individuals in the vcf file. After running, the next
    line from getLine() should be the first locus line.
    If write is true, will write the individual names to outfile.
*/
int getIndividuals(std::ifstream &infile, bool write, std::ofstream &outfile);

/*
    Gets the number of loci in the vcf file by counting lines not beginning
    with '#'.
*/
int getLoci(std::ifstream &infile);

/*
    Converts vcf data to data in the genotype likelihood format.
*/
bool convertGl(std::ifstream &infile, std::ofstream &outfile);

#endif
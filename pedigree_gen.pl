#!/usr/bin/perl

# Peter Lazorchak
# lazorchakp@gmail.com
# 6/16/17

# Highly specialized script for converting individual names map to a pedigree
# (in the LINKAGE format) to be used with lep-MAP3

# Modify the regexps for differently formatted input files
# General format:
# <ind_name><separator><crossinfo>
# example:  MID27.1_s_6_sequence    cross14-mother

use strict;
use warnings;
use Getopt::Long;
use constant {MALE => 1, FEMALE => 2, NA => 0};

my ($infile, $outfile);

GetOptions(
    'i=s'   =>  \$infile,
    'o=s'   =>  \$outfile
);

die "Must specify infile (-i) and outfile (-o). See comments for info.\n"
    if (!defined($infile) || !defined($outfile));

open(IN, $infile) or die "Cannot open $infile\n";

my (@inds, %mothers, %fathers);
my $nind = 0;

while (<IN>) {
    my @line = split;
    $line[1] =~ /([a-z]*\d+)[\.-]?(.+)/i;
    my $family = $1;
    my $id = $2;
    # store family, ind_name, father, mother, and sex
    $inds[$nind] = [$family, $line[0], NA, NA, NA];
    if ($id =~ /(mother|female)/) {
        $mothers{$family} = $line[0];
        $inds[$nind][4] = FEMALE;
    } elsif ($id =~ /(father|male)/) {
        $fathers{$family} = $line[0];
        $inds[$nind][4] = MALE;
    }
    ++$nind;
}

close(IN);

open(OUT, ">$outfile") or die "Cannot open $outfile\n";

foreach my $param (0..4) { # for each row in the output file (inner array)
    print OUT ".\t.";
    foreach my $k (0..$nind - 1) { # for each individual (outer array)
        if (!($inds[$k][4] eq NA)) { # if this is a parent
            print OUT "\t$inds[$k][$param]";
        } else {
            if ($inds[$k][2] eq NA) { # if we don't have parent data yet
                $inds[$k][2] = $fathers{$inds[$k][0]};
                $inds[$k][3] = $mothers{$inds[$k][0]};
            }
            print OUT "\t$inds[$k][$param]";
        }
    }
    print OUT "\n";
}

# we don't know any phenotypes, so print all 0s
print OUT ".\t.";
foreach my $k (0..$nind - 1) {
    print OUT "\t".NA;
}
print OUT "\n";

close(OUT);

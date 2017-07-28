#!/usr/bin/perl

# Peter Lazorchak
# lazorchakp@gmail.com
# 6/14/17

# Creates a populations file for use with CLUMPAK
# This would be easy, but the Grant apple & haw individuals have names that
# begin with MID and don't contain any site/race info.
# This program requires a Grant map file (grantfile) with the general format
# <IND_NAME>  <Site><Host>

# This program has regexps highly specific to the two input files. It is
# recommended that these be verified/modified for different data sets.

use strict;
use warnings;
use Getopt::Long;
use File::Spec;
use File::Basename;
use feature "switch";

my ($infile, $grantfile, $outfile);
GetOptions(
    'i=s'   =>  \$infile,
    'g=s'   =>  \$grantfile,
    'o=s'   =>  \$outfile,
);

if (!defined($infile) || !defined($grantfile) || !defined($outfile)) {
    die "Must specify names file (-i), grantfile (-g), & outfile (-o)\n";
}

my %midhash;
open(MID, $grantfile) or die ("Can't open grantfile");
while (<MID>) {
    my @line = split;
    $line[1] =~ /^([A-Z][a-z]+)([A-Z][a-z]+$)/;
    my $host = lc $2;
    $midhash{$line[0]} = "$1"."_$host";
}
close(MID);

open(NAMES, $infile) or die ("Can't open names file");
open(OUT, ">$outfile") or die ("Can't open out file");
while (<NAMES>) {
    if (/(haw|app)(le)?\.(urbana|dowagiac|BrazosBend|SouthBend)/i) {
        my $host = lc $1;
        my $site = ucfirst $3;
        print OUT "$site"."_$host\n";
    } elsif (/Fenn\.(haw|app)/i) {
        my $host = lc $1;
        print OUT "Fenn_$host\n";
    } elsif (/^cf\d*/i) {
        my $site = "unknown";
        for ($_) {
            when (/^cf\d*\.a2?\./i) {$site = "SFA"} # A. or A2.
            when (/^cf\d*\.g(2|ai)\./i) {$site = "GAI"} # G2. or GAI.
            when (/^cf\d*\.kn\./i) {$site = "KN"} # KN.
            when (/^cf\d*\.rfl\./i) {$site = "RFL"} # RFL.
            when (/^cf\d*\.i\./i) {$site = "I57"} # i.
            when (/^cf\d*\.k\./i) {$site = "KIS"} # k.
            when (/^cf\d*\.nc\./i) {$site = "NC"} # nc.
            when (/^cf\d*\.md\./i) {$site = "MD"} # md.
            when (/^cf\d*\.sc\./i) {$site = "SC"}
            default {die "unidentified site for dogwood fly $_\n"}
        }
        print OUT "${site}_dog\n";
    } elsif (/^MID/) {
        chomp;
        print OUT "$midhash{$_}\n";
    } else {
        chomp;
        die "$_ not identified\n";
    }
}
close(NAMES);
close(OUT);

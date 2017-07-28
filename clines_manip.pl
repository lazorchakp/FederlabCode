#!/usr/bin/perl

use strict;
use warnings;

print "proceed? ";
<STDIN> =~ /^y/i  or die "stopped\n";
print "\nrunning...\n\n";

foreach my $sig ('sigmap', 'sigall') {
    foreach my $site ( qw(
            Dow
            Fen
            Gra
            GU
            GU.app
            GU.haw
            Urb
        ) ) {

        my $set = "$sig.$site";
        print "$set\n";
        0 == system "clump_gen.pl -i entropy/$set -o clumpak/$set -n names/Clines2017.$site.txt -g Grant.app.haw.txt -k 1 -K 3 -c colors_file.txt -l labels_file.txt" or die "$!\n";
    }
}

=pod

3
Dow
Fen
Gra
GU
GU.app
GU.haw
Urb

8
all

=cut
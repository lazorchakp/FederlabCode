#!/usr/bin/perl

use strict;
use warnings;

print "proceed? ";
<STDIN> =~ /^y/i  or die "stopped\n";
print "\nrunning...\n\n";

foreach my $set ( qw(
        SB.KN.recode
        SFA.BB.recode.ld_filtered
        SFA.BB.recode
        SB.KN.recode.ld_filtered
        SB.KN.AHclineAll
        SFA.BB.AHclineAll
        SB.KN.AHclineMap
        SFA.BB.AHclineMap
    ) ) {

    print "$set\n";
    0 == system "clump_gen.pl -i entropy/$set -o clumpak/$set -n names/$set.txt -g Grant.app.haw.txt -k 1 -K 3 -c colors_file.txt -l labels_file.txt" or die "$!\n";
}

=pod
        6
        Pom.Dog.bi.clines.75.recode
        Pom.Dog.bi.clines.75.recode.ld_filtered
        Pom.Dog.bi.clines.75.AHclineAll
        Pom.Dog.bi.clines.75.AHclineMap

        3
        SB.KN.recode
        SFA.BB.recode.ld_filtered
        SFA.BB.recode
        SB.KN.recode.ld_filtered
        SB.KN.AHclineAll
        SFA.BB.AHclineAll
        SB.KN.AHclineMap
        SFA.BB.AHclineMap

=cut
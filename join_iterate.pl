#!/usr/bin/env perl

# Peter Lazorchak
# lazorchakp@gmail.com
# 6/16/17

use strict;
use warnings;
use Getopt::Long;

my ($datafile, $mapfile, $cp, $lodLim, $jump, $nsteps, $multioutput);
my $maleTheta = 0;
my $femaleTheta = 0.09;

GetOptions(
    'i=s'   =>  \$datafile,
    'm=s'   =>  \$mapfile,
    'c=s'   =>  \$cp,
    's=f'   =>  \$lodLim,
    'j=f'   =>  \$jump,
    'n=i'   =>  \$nsteps,
    'a=f'   =>  \$maleTheta,
    'f=f'   =>  \$femaleTheta,
    'o=i'   =>  \$multioutput
);

if (!defined $datafile || !defined $mapfile || !defined $lodLim || !defined $jump
        || !defined $nsteps) {
    die "Must specify:\n"
        ."\tdata.gz (-i), mapfile.txt (-m), starting LOD (-s),\n"
        ."\tjump size (-j), number of steps (-n)\n"
        ."Optional:\n"
        ."\tclasspath (-c), output intermediate maps every n steps (-o),\n"
        ."\tmale theta (-a) [0], female theta (-f) [0.09]\n";
}

if (defined $cp) {
    $cp = "-cp ".$cp;
} else {
    $cp = '';
}
my $mapbase = $mapfile;
$mapbase =~ s/\.txt//g;

for (my $step = 0; $step < $nsteps; ++$step) {
    my $nextmap = "$mapbase"."_join$lodLim"."by$jump".".txt";
    system(
        "zcat $datafile | "
        ."java $cp JoinSingles2All data=- map=$mapfile lodLimit=$lodLim "
        ."maleTheta=$maleTheta femaleTheta=$femaleTheta "
        ."> $nextmap"
    );
    # remove the temp map unless it's the first step or it's one of the maps
    # we are going to save
    if ($step != 0 && (!defined $multioutput || $step % $multioutput != 0)) {
        unlink glob $mapfile;
    }
    $lodLim -= $jump;
    $lodLim = substr $lodLim, 0, 5;
    $mapfile = $nextmap;
}

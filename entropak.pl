#!/usr/bin/perl


# Peter Lazorchak
# lazorchakp@gmail.com
# 6/13/2017


=pod

DESCRIPTION
    This script generates CLUMP files and runs them.


DEPENDENCIES
    Unless specified, clump_gen.pl will look for the following in $PATH.

Required programs:
    estpost_split
    CLUMPAK.pl
        note that this has many of its own dependencies - consult online manual
            one particularly frustrating one is the GD module. To install this,
            you first must install libgd (gdlib-config must be in your $PATH)
    BestKByEvanno.pl
Required if using input files option b (below):
    popfile_gen.pl

Required modules:
    hdf5
    gsl

CLUMPAK.pl and BestKByEvanno.pl rely on modules in the clumpak distribution.
Moving these scripts out of the clumpak directory is not recommended. Also,
all clumpak scripts must be run out of this directory, so this script cd's into
the specified clumpak directory immediately before executing.


FILES
Required input files:
    Entropy output - *_k2_1.hdf5, *_k2_2.hdf5, *_k3_1.hdf5, *_k3_2.hdf5 etc 
    ONE of the following (a OR b, recommended = a):
    a)  Populations file - single column, each line with the population
            corresponding to an individual. These must be in the same order
            as the original vcf file's individuals.
        {
            <pop_for_ind0>
            <pop_for_ind1>
            ...
        }
    b)  Individual names file
        {
            <ind0>
            <ind1>
            ...
        }
        AND
        Grant population map file - for flies from Grant
        {
            MID<ind0>   Grant<App|Haw>
            MID<ind1>   Grant<App|Haw>
            ...
        }
        note that option b relies on very specific regular expressions to
        determine populations based on individuals' names. Option a is thus
        recommended, but b is more convenient if you know popfile_gen.pl works.

Optional input files:
    Colors file - list of colors to be used by distruct (see distruct manual)

    drawparams file - visual specifications for distruct (see distruct manual)
        The default is contained within the distruct directory in clumpak.
        It may be simpler to modify this file than to specify a new one. Note
        when using this file that clumpak makes many modifications and ignores
        certain options (such as K). For INDIVWIDTH, use a negative value to
        allow clumpak to calculate page-fitting width.

    labels file - list of population order for the final distruct graphs.
        Excess populations in this file are permitted, and extras will be
        ignored. Thus, for a group of datasets (ex. Clines2017), only one
        labels file is needed. This file is NOT THE SAME as the labels_file
        from clumpak (or INFILE_LABEL_BELOW from distruct). It is used to
        generate this file.
    {
        <pop0>
        <pop1>
        ...
    }

=cut

use warnings;
use strict;
use Getopt::Long qw(:config no_ignore_case);
use File::Basename;
use File::Path qw(make_path);
use File::Spec;

&author;

my $indir;
my $outdir;
my $popfile;
my $namefile;
my $grantfile;

my $mink=2;
my $maxk=3;
my $burnin=0;

my $colorsfile;
my $drawparams;
my $labelsfile;

my $estpost = 'estpost_split';
my $clumpak = "$ENV{HOME}/bin/clumpak";
my $popfile_gen = 'popfile_gen.pl';

GetOptions(
    'i=s' => \$indir,
    'o=s' => \$outdir,
    'p=s' => \$popfile,
    'n=s' => \$namefile,
    'g=s' => \$grantfile,
    'k=i' => \$mink,
    'K=i' => \$maxk,
    'b=i' => \$burnin,
    'c=s' => \$colorsfile,
    'd=s' => \$drawparams,
    'l=s' => \$labelsfile,
    'e=s' => \$estpost,
    'u=s' => \$clumpak,
    'f=s' => \$popfile_gen,
    'h'   => \&usage,
);

&usage unless (defined $indir && defined $outdir);
&usage unless (defined $popfile || (defined $namefile && defined $grantfile));

my $nks = $maxk - $mink + 1;


# input dir
if (! -e $indir) {
    die ("\nCan't access input directory: $indir\n");
}
$indir = File::Spec->rel2abs($indir);

# output dir
if (! -e $outdir) {
    eval {make_path($outdir)}
        or die ("\nCan't create output directory: $outdir\n");
} else { # check for empty directory
    opendir DIR, $outdir or die("\nCan't open output directory: $outdir\n");
    my @contents = readdir DIR;
    die ("\nOutput directory must be empty\n") if (scalar @contents > 2); # . and ..
    closedir DIR;
}
$outdir = File::Spec->rel2abs($outdir);

# clumpak dir
die "Can't find clumpak directory: $clumpak\n" unless (-d $clumpak);
$clumpak = File::Spec->rel2abs($clumpak);

my $prerunDir = $outdir."/prerun";
my $clumpoutDir = $outdir."/clumpak_out";
my $bestKout = $outdir."/bestk_out";

-d $prerunDir or mkdir $prerunDir
    or die "Cannot make directory: $prerunDir\n";
-d $clumpoutDir or mkdir $clumpoutDir
    or die "Cannot make directory: $clumpoutDir\n";
-d $bestKout or mkdir $bestKout
    or die "Cannot make directory: $bestKout\n";

# check for popfile, or namefile and grantfile
if (defined $popfile) {
    if (! -e $popfile) {
        die "\nCan't find populations file: $popfile\n\n";
    }
} elsif (! -e $namefile) {
    die "\nCan't find names file: $namefile\n\n";
} elsif (! -e $grantfile) {
    die "\nCan't find Grant populations map file: $grantfile\n\n";
}

# set up the popfile if it doesn't exist, and get nind by counting its lines
if (defined $popfile) {
    $popfile = File::Spec->rel2abs($popfile);
} else {
    $namefile = File::Spec->rel2abs($namefile);
    $grantfile = File::Spec->rel2abs($grantfile);
    $popfile = "$outdir/populations_file.txt";
    0 == system "$popfile_gen -i $namefile -g $grantfile -o $popfile"
        or die "Failed to produce populations file\n";
}
my $nind = `wc -l < $popfile`;
chomp $nind;

# get the names of all input files
opendir DIR, $indir or die("\nCan't open input directory: $indir\n");
    my @hdf5files = grep(/\.hdf5$/, readdir DIR);
    foreach my $h (@hdf5files) {
        $h = "$indir/$h";
    }
closedir DIR;

decodeHDF5($mink, $maxk, $prerunDir, @hdf5files);

rearrange_output($mink, $maxk, $prerunDir, $nind);

0 == system "zip -r -j $outdir/prerun.zip $prerunDir/*_K*_*.dat" or die "$!\n";

# set up clumpak
my $clumpCommand = "./CLUMPAK.pl --id $$ --inputtype admixture "
    ."--dir $clumpoutDir --file $outdir/prerun.zip --indtopop $popfile ";
$clumpCommand .= "--colors ".File::Spec->rel2abs($colorsfile)." "
    if (defined $colorsfile);
$clumpCommand .= "--drawparams ".File::Spec->rel2abs($drawparams)." "
    if (defined $drawparams);

# convert the labels file to the proper format
if (defined $labelsfile) {
    my $labelsfile_conv = &convert_labelsfile($outdir, $popfile, $labelsfile);
    $clumpCommand .= "--labels ".File::Spec->rel2abs($labelsfile_conv)." ";
}

# cd into clumpak directory before running
# all filenames have absolute paths so this will not cause problems
print "Running CLUMPAK\n";
0 == system "cd $clumpak; $clumpCommand"
    or die "Failed to run CLUMPAK ($!)\n";

# execute BestKByEvanno
print "Running BestKByEvanno\n";
0 == system "cd $clumpak; ./BestKByEvanno.pl --id $$ --dir $bestKout "
    ."--file $outdir/log_prob_file.txt --inputtype lnprobbyk "
    or die "Failed to run BestKByEvanno ($!)\n";

# ============================== SUBROUTINES ==================================

# run estpost_split
sub decodeHDF5 {
    my $mink=shift;
    my $maxk=shift;
    my $prerunDir=shift;
    my @hdf5files=@_;

    # Write a single usable csv for each group of hdf5 files with the same k
    foreach my $k ($mink..$maxk) {
        my @hdf5k = grep(/\_k$k(\_[0-9]+)*\.hdf5$/, @hdf5files);
        my $input = join(' ', sort (@hdf5k));
        my $outfile = "$prerunDir/".basename($hdf5k[0]);
        $outfile =~ s/(\_[0-9]+)*\.hdf5$/\.q_postest_split.csv/g;
        $outfile =~ s/\.q_postest/\.b$burnin\.q_postest/g if ($burnin > 0);

        # generate intermediate files. Ideally estpost would output the correct
        # format (TODO future project), but for now we output them, then modify
        # Additionally, create the log probabilities file.
        system("$estpost -b $burnin -o $outfile $input "
            ."| awk '{print \"$k\t\"\$1}' >> $prerunDir/../log_prob_file.txt");
    }
}

=pod
Reshapes the data into multiple output files.
Originally, (after the header), the data is formatted as:


FILENAME: <info>_kP<moreinfo>_split.csv

q_ind_0_pop_0,mean1,mean2...meanM
...
q_ind_N-1_pop_0,mean1,mean2...meanM
q_ind_0_pop_1,mean1,mean2...meanM
...
q_ind_N-1_pop_P-1,mean1,mean2...meanM

for M reps, N individuals, and P populations (k = P)
The new format is (for example):


FILENAME: <info>_K3_2.dat

ind0_pop0_mean2 ind0_pop1_mean2 ind0_pop2_mean2
ind1_pop0_mean2 ind1_pop1_mean2 ind1_pop2_mean2
...
indN-1_pop0_mean2 indN-1_pop1_mean2 indN-1_pop2_mean2


Note that the output files are whitespace delimited
This is the format required by CLUMPAK
=cut

sub rearrange_output {
    my $mink = shift;
    my $maxk = shift;
    my $prerunDir = shift;
    my $nind = shift;

    opendir(DIR, $prerunDir) or die("Cannot open directory: $prerunDir\n");
    my @csvfiles=grep(/_k\d+.*_split.csv$/, readdir(DIR));
    foreach my $filename (@csvfiles) {
        $filename="$prerunDir/$filename";
    }
    closedir(DIR);

    foreach my $k ($mink..$maxk) {
        my ($currfile) = grep(/_k$k/, @csvfiles);
        die "$!\n(Did you forget to load the gsl or hdf5 modules?)\n"
            if (!defined $currfile);
        my $outfile ="$prerunDir/".basename($currfile);
        $outfile =~ s/_k$k.+_split.csv$/_K$k/g; #will eventually add _<rep>.dat

        # load all means into the memory
        open(FILE, "$currfile") or die ("\nCan't open file: $currfile\n");
        my $header=<FILE>;
        chomp($header);
        my ($nreps) = $header =~ /(\d+)$/;
        my @means;
        my $currind=0;
        while (<FILE>) {
            chomp;
            s/^q_ind_\d+_pop_\d+,//g; # remove individual id
            $means[$currind] .= "$_,";
            $currind = ++$currind % $nind;
        }
        close(FILE);
        unlink glob ("$currfile");

        foreach my $rep (1..$nreps) {
            open(OUT, ">$outfile"."_$rep.dat")
                or die("\nCan't open file: $outfile"."_$rep.dat\n");
            --$rep;
            foreach my $ind (0..($nind - 1)) {
                foreach my $pop (0..($k - 1)) {
                    my @indprobs = split(/,/,$means[$ind]);
                    print OUT "$indprobs[$rep + $nreps * $pop] ";
                }
                print OUT "\n";
            }
            close(OUT);
        }
    }
}

sub convert_labelsfile {
    my ($outdir, $popfile, $labelsfile) = @_;
    open(POPS, $popfile)
        or die "Failed to open pops file during labels file conversion\n";
    my %popids;
    my $currid = 1;
    while (<POPS>) {
        chomp;
        if (!exists $popids{$_}) {
            $popids{$_} = $currid;
            ++$currid;
        }
    }
    close POPS;
    open(LABELS, $labelsfile) or die
        "Failed to open original labels file during labels file conversion\n";
    my $labelsfile_conv = $outdir."/conv_".basename($labelsfile);
    open(CONV, ">$labelsfile_conv") or die
        "Failed to open new labels file during labels file conversion\n";
    my $linecount = 0;
    my $npops = keys %popids;
    while (<LABELS>) {
        chomp;
        if (exists $popids{$_}) {
            print CONV "$popids{$_}\t$_\n";
            delete $popids{$_};
            ++$linecount
        }
    }
    if ($linecount < $npops) {
        warn "Fatal: Too few valid populations in labels file\n"
            ."Expected: $npops\tFound $linecount\n";
        warn "The following individuals were not found in the labels file:\n";
        foreach (keys %popids) {
            warn "$_\n";
        }
        exit 1;
    }
    close LABELS;
    return $labelsfile_conv;
}

sub author {
    print "\n"
    ."#########################################\n"
    ."  ".basename($0)."\n"
    ."  Peter Lazorchak\n"
    ."  lazorchakp\@gmail.com\n"
    ."#########################################\n\n"
}

sub usage {
    print "  Usage:\n"
    ."    ".basename($0)."\n"
    ."      -i  <input directory (with *_k[0-9]+(\_[0-9]+){0,1}.hdf5 files)>\n"
    ."      -o  <output directory>\n"
    ."      -p  <populations file>\n"
    ."      -n  <names file> required if no populations file\n"
    ."      -g  <Grant population map file> required if no populations file\n"
    ."      -k  <min. K> [2]\n"
    ."      -K  <max. K> [3]\n"
    ."      -b  <extra burnin> [0]\n"
    ."      -c  <colors file> [none]\n"
    ."      -d  <drawparams file> [none]\n"
    ."      -l  <labels file> [none]\n"
    ."      -e  <path to estpost> [estpost_split]\n"
    ."      -u  <path to clumpak directory> [\$HOME/bin/clumpak]\n"
    ."      -f  <path to popfile_gen.pl> [popfile_gen.pl]\n"
    ."      -h  <show this help>\n"
    ."\n";
    exit;
}

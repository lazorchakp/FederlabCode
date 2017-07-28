#!/usr/bin/env perl

# Peter Lazorchak
# lazorchakp@jhu.edu
# 6/8/17

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case);
use File::Spec;
use File::Basename;

my $jobname;
my $infile;
my $hdf5file;
my $shell='bash';
my $email='plazorch@nd.edu';
my $ncores=1;
my $nsteps=10000;
my $burnin=1000;
my $thin=1;
my $mink=1;
my $maxk=6;
my $nreps=3;
my $submit;

GetOptions(
    'j=s'   =>  \$jobname, # base name of each .job, and .o file
    'i=s'   =>  \$infile, # in .gl format
    'o=s'   =>  \$hdf5file, # base name of hdf5 file to output
    's=s'   =>  \$shell, # the shell to use (zsh, csh, etc)
    'M=s'   =>  \$email, # email address where updates will be sent
    'c=i'   =>  \$ncores, # number of cores to use per job
    'l=i'   =>  \$nsteps, # number of mcmc steps
    'b=i'   =>  \$burnin, # number of mcmc steps to ignore initially
    't=i'   =>  \$thin, # only sample on every $thin step
    'k=i'   =>  \$mink, # minimum assumed k-value
    'K=i'   =>  \$maxk, # maximum assumed k-value
    'r=i'   =>  \$nreps, # number of runs for each k-value
    'x=s'   =>  \$submit, # actually submit jobs, and submit them from this dir
) or die "Error: unrecognized command line argument\n";

die "\nMandatory parameters:\n"
    ."\t-j  base name of job files: don't include .job extension\n"
    ."\t-i  input file name, in genotype likelihood format\n"
    ."Available options:\n"
    ."\t-o  base name of hdf5 files to output: don't include .hdf5 extension"
    ." [jobname]\n"
    ."\t-s  shell for job submission [bash]\n"
    ."\t-M  email address where job updates will be sent [plazorch\@nd.edu]\n"
    ."\t-c  number of cores to use per job [1]\n"
    ."\t-l  number of mcmc steps [10000]\n"
    ."\t-b  burnin: number of mcmc steps ignored [1000]\n"
    ."\t-t  thin: only save every nth sample [1]\n"
    ."\t-k  minimum k-value to run [1]\n"
    ."\t-K  maximum k-value to run [6]\n"
    ."\t-r  number of reps to perform at each k-value [3]\n"
    ."\t-x  submit jobs from this directory [do not submit any jobs]\n"
    ."\t    The .o log files will appear in the specified directory.\n"
    ."\t    entropy1.2 uses the process ID in the random seed, so simultaneous"
    ."\n\t    submissions will still be uniquely random.\n"
    ."Note:\n\tWhen specifying file base names, include the desired path.\n"
    ."\tDefault location is the current working directory.\n\n"
    if (!defined($jobname) || !defined($infile));

if (!defined($hdf5file)) {
    $hdf5file=$jobname;
}

$infile=File::Spec->rel2abs($infile);
$hdf5file=File::Spec->rel2abs($hdf5file);
$jobname=File::Spec->rel2abs($jobname);

my $origdir = `pwd`;
chomp $origdir;
if (defined $submit) {
    chdir $submit or die "Cannot cd into $submit\n";
}

for (my $k = $mink; $k <= $maxk; ++$k) {
    for (my $rep = 1; $rep <= $nreps; ++$rep) {
        my $currjob = $jobname.'_k'.$k."_$rep";
        my $basejob = basename($currjob);
        my $currhdf5 = $hdf5file.'_k'.$k."_$rep.hdf5";
        open(JOB, "> $currjob.job")
            or die "Error: could open $currjob.job\n";
        print JOB <<EOF;
#!/bin/$shell

#\$ -N $basejob
#\$ -pe smp $ncores
#\$ -M $email
#\$ -m ae

module load hdf5 gsl
entropy -i $infile -l $nsteps -b $burnin -t $thin -m 1 -k $k -o $currhdf5
EOF
        close(JOB);
        if (defined($submit)) {
            system("qsub $currjob.job");
        }
    }
}

chdir $origdir;

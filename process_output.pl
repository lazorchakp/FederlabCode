#!/usr/bin/env perl

# (c) Victor Soria-Carrasco
# victor.soria.carrasco@gmail.com
# Last modified: 02/12/2016 18:05:38

# Description:
# This script plots barplots for q
# estimated by entropy, as well as
# DIC plots.
# Entropy output files expected are:
# *_k2_1.hdf5, *_k2_2.hdf5, *_k3_1.hdf5, *_k3_2.hdf5 etc

# Samples ids file must have the following format:
# sp_loc_id

# WARNING: Beware this script may need large amounts of memory to run (R code)

use warnings;
use strict;
use Getopt::Long;
use File::Basename;
use File::Path qw(make_path);
use File::Spec;
use Statistics::R;
use Sort::Naturally;

my $estpost='estpost';

my $version='1.0-2016.12.02';
my $ncpu=1;

&author;
my $mink=2;
my $maxk=3;
my $burnin=0;
my ($indir, $outdir, $smpsfile);
GetOptions(
    'i|I=s'  => \$indir,
    'o|O=s'  => \$outdir,
    'e|E=s'  => \$estpost,
    's|S=s'  => \$smpsfile,
	'mink=i' => \$mink,
	'maxk=i' => \$maxk,
	'b|B=i'  => \$burnin,
    'h|help' => \&usage
);

&usage if (!defined($indir) || !defined($outdir));

my $nks=$maxk-$mink+1;


# Input directory
# -------------------------------------------------------------------------------------
if (! -e $indir){
	die ("\nCan't access input file: $indir\n\n");
}
$indir=File::Spec->rel2abs($indir);
# -------------------------------------------------------------------------------------

# Output directory
# -------------------------------------------------------------------------------------
if (! -e $outdir){
	eval {make_path($outdir)}
		or die ("\nCan't create output directory: $outdir\n\n");
}
$outdir=File::Spec->rel2abs($outdir);
# -------------------------------------------------------------------------------------

opendir (DIR, $indir)
	or die("\nCan't open input directory: $indir\n\n");
	my @hdf5files=grep(/\.hdf5$/, readdir(DIR));
	foreach my $h (@hdf5files) {
		$h="$indir/$h";
	}
closedir (DIR);

print "Processing... This could take a minute\n\n";

# Generate table with q posterior median estimates, and get number of loci
# and number of samples
my ($x, $y, $z)=make_q_table($mink, $maxk, $outdir, $smpsfile, @hdf5files);
my $nloci=$$x;
my $nsmps=$$y;
my @q_tables=@$z;

# Generate table with deviance, number of parameters, and
# DIC (Deviance Information Criterion)

my $bname=basename($hdf5files[0]);
$bname =~ s/\_k[0-9]+(\_[0-9]+){0,1}\.hdf5$//g;
$bname.=".b$burnin" if ($burnin > 0);
my $outfile="$outdir/$bname.deviance_DIC.csv";

make_DIC_table($mink, $maxk, $outfile, @hdf5files);

# Make plots using R
my $R = Statistics::R->new(shared=>1);

# plot Q plots
foreach my $k ($mink..$maxk){

	my @hdf5=grep(/.*\_k$k(\_[0-9]+){0,1}\.hdf5$/, sort @hdf5files);
	print "Processing files for K=$k (".scalar(@hdf5)." runs)...\n\n";

	my $output=basename($hdf5[0]);
	$output=~ s/(\_[0-9]+){0,1}\.hdf5/\.q_postest\.csv/g;
	$output=~ s/\.q_postest\.csv/\.b$burnin\.q_postest\.csv/g if ($burnin > 0);
	$output="$outdir/$output";

	my $ofile=basename($output);
	$ofile=~ s/\.csv//g;
	$ofile="$outdir/$ofile";

	plot_q($k,$output,$ofile,$nsmps);

	print "\n\tPlot saved as $ofile.pdf\n\n";
}

# plot DIC plot
my $ifiledic="$outdir/$bname.deviance_DIC.csv";
my $ofiledic="$outdir/$bname.deviance_DIC.pdf";
plot_DIC($ifiledic,$ofiledic);

print "\nDIC plot saved as $ofiledic\n\n";

$R->stop();


# ==============================================================================
# ==============================================================================
# ============================== SUBROUTINES ===================================
# ==============================================================================
# ==============================================================================

# Show copyright
# ==============================================================================
sub author{
    print "\n";
    print "#########################################\n";
    print "  ".basename($0)."\n";
	print "  version $version\n";
    print "  (c) Victor Soria-Carrasco             \n";
    print "  victor.soria.carrasco\@gmail.com      \n";
    print "#########################################\n";
	print "\n";
}
# ==============================================================================

# Show usage
# ==============================================================================
sub usage{
    print "\n";
	print "  Usage:\n";
    print "    ".basename($0)."\n";
	print "      -i <input directory (with *_k[0-9]+(\_[0-9]+){0,1}.hdf5 files)>\n";
	print "      -o <output directory>\n";
	print "      -e <path to estpost> (optional, default=estpost)\n";
	print "      -s <samples ids file (optional)>\n";
	print "      -mink <min. K (optional, default=2)>\n";
	print "      -maxk <max. K (optional, default=3)>\n";
	print "      -b <extra burnin (optional, default=0)>\n";
	print "      -h <show this help>\n";
    print "\n";
    exit;
}
# ==============================================================================

# Create table with deviance and DIC
# ==============================================================================
sub make_DIC_table{
	my $mink=shift;
	my $maxk=shift;
	my $outfile=shift;
	my @hdf5files=@_;

	system ("echo 'K,deviance,effective_no_parameters,DIC' > $outfile"); # header
	foreach my $k ($mink..$maxk){
		my @hdf5k=grep(/\_k$k(\_[0-9]+){0,1}\.hdf5$/, @hdf5files);
		my $input=join(' ', sort (@hdf5k));

		# Calculate DIC
		system("$estpost -b $burnin -p deviance -s 3 -o /dev/null $input | \\
			grep -P 'Model deviance|Effective number of parameters|Model DIC' | \\
			perl -pi -e 's/.*\: //g; s/\n/,/g;' | \\
			perl -pi -e 's/,\$/\n/g; s/^/$k,/g;' \\
			>> $outfile");
	}
}
# ==============================================================================

# Create table with posterior median estimates of q for each k
# ==============================================================================
sub make_q_table{
	my $mink=shift;
	my $maxk=shift;
	my $outdir=shift;
	my $smpsfile=shift;
	my @hdf5files=@_;

	# Samples file (if available)
	my @ids=();
	if (defined($smpsfile)){
		open (FILE, "$smpsfile")
			or die ("\nCan't open populations file: $smpsfile\n\n");
			while (<FILE>){
				s/\r$//g;# DOS/old Mac text files
				chomp;
				push(@ids,$_);
			}
		close (FILE);
	}


	# Get number of loci and number of samples
	my $tmp=`$estpost -b $burnin -p gprob $hdf5files[0] -o /dev/null | grep "parameter dimensions for gprob"`;
	$tmp=~ /loci \= ([0-9]+)\,/;
	my $nloci=$1;
	$tmp=~ /ind \= ([0-9]+)\,/;
	my $nsmps=$1;

	# Get admixture proportions
	my @q_tables=();
	foreach my $k ($mink..$maxk){
		my @hdf5k=grep(/\_k$k(\_[0-9]+)*\.hdf5$/, @hdf5files);
		my $input=join(' ', sort (@hdf5k));
		my $outfile="$outdir/".basename($hdf5k[0]);
		$outfile=~ s/(\_[0-9]+)*\.hdf5$/\.q_postest/g;
		$outfile=~ s/\.q_postest/\.b$burnin\.q_postest/g if ($burnin > 0);

		# Summarize q
		system("$estpost -b $burnin -p q -s 0 -o $outfile.tmp $input >& $outfile.log");

		open (FILE, "$outfile.tmp")
			or die ("\nCan't open file: $outfile.tmp\n\n");
			my $header=<FILE>; # Get rid of header
			my %output;
			while (<FILE>){
				my @aux=split(/[\,|\s]/,$_);
				$aux[0]=~ s/(\_[0-9]+)*$//g;
				$output{$aux[0]}.=$aux[1].','; # mean
				# $output{$aux[0]}.=$aux[2].','; # median
			}
		close (FILE);

		# Output file for this K
		# ------------------------------------------------------

		# Sort as in the original file
		my @smps=sort {
			my ($a1)=$a=~ /^q_ind_([0-9]+)_pop$/;
			my ($b1)=$b=~ /^q_ind_([0-9]+)_pop$/;
		$a1 <=> $b1;
		} keys %output;

		# Prepare output
		my %out=();
		foreach my $i (0..$#smps){
			my $out=$output{$smps[$i]};
			$out=~ s/\,$//g;
			my $id=$smps[$i];
			$id=~ s/\_pop//g;
			$id=$ids[$i] if (defined($smpsfile));
			$out{$id}="$id,$out\n";
		}

		open (FILEOUT, ">$outfile.csv")
			or die ("\nCan't write to file $outfile.csv\n\n");

			# New header
			$header="sample";
			foreach my $q (1..$k){ $header.=",k$q"; }
			print FILEOUT "$header\n";

			# Output median estimates
			# (re-sorting by sample id)
			foreach my $o (nsort keys %out){
				print FILEOUT $out{$o};
			}

		close (FILEOUT);

		unlink glob ("$outfile.tmp");
		push (@q_tables, "$outfile.csv");
	}

	return(\$nloci,\$nsmps,\@q_tables);
}
# ==============================================================================


# Plot q barplot
# ==============================================================================
sub plot_q{
	my $k=shift;
	my $infile=shift;
	my $output=shift;
	my $nsmps=shift;

	$R->run(qq`
		col.palette<-c(
			"#E31A1C", # red
			"green4",
			"dodgerblue2",
            "#6A3D9A", # purple
            "#FF7F00", # orange
            "black","gold1",
            "skyblue2","#FB9A99", # lt pink
            "palegreen2",
            "#CAB2D6", # lt purple
            "#FDBF6F", # lt orange
            "gray70", "khaki2",
            "maroon","orchid1","deeppink1","blue1","steelblue4",
            "darkturquoise","green1","yellow4","yellow3",
            "darkorange4","brown")

		q<-read.table("$infile", sep=",", header=T)

		sp<-as.factor(unlist(lapply (q[,1], function(x) unlist(strsplit(as.character(x),'_'))[1])))
		loc<-as.factor(unlist(lapply (q[,1], function(x) unlist(strsplit(as.character(x),'_'))[2])))
		host<-as.factor(unlist(lapply (q[,1], function(x) unlist(strsplit(as.character(x),'_'))[3])))
		id<-unlist(lapply (q[,1], function(x) unlist(strsplit(as.character(x),'_'))[4]))

		# png(filename="$output.png", width=(10*$nsmps), height=800, res=150)
		mult<-1
		if (length(q[,1]) < 100) mult<-2
		if (length(q[,1]) < 50) mult<-3
		if (length(q[,1]) < 10) mult<-5

		pdf(file="$output.pdf", width=(mult*(8*$nsmps)/72), height=(600/72), pagecentre=F)
		par(mar=c(10,3,3,2), oma=c(2,2,2,2))

		mp<-barplot(
			t(as.matrix(q[,-1])),
			main="K=$k",
			axes=F,
			xlim=c(0,length(q[,1])),
			ylim=c(-0.05,1.25),
			inside=F, beside=F,
			border=col.palette[1:$k],
			col=col.palette[1:$k],
			space=0, las=2, cex.main=1.5)
		mtext(text="(ids=$nsmps, loci=$nloci)",line=-0.5,cex=1)
		axis(side=2,at=c(0,0.2,0.4,0.6,0.8,1.0))

		# text(mp, par("usr")[3], labels=q[,1], font=4, cex=0.6, xpd=T, srt=45, adj=1)
		# text(mp, par("usr")[3], labels=q[,1], col=c(col.palette[as.numeric(loc)]), font=4, cex=0.6, xpd=T, srt=45, adj=1)
		text(mp, par("usr")[3], labels=id, font=4, cex=0.6, xpd=T, srt=45, adj=1)
		# legend(length(q[,1])+1, 1 ,legend=unique(loc), fill=c(unique(col.palette[as.numeric(loc)])), cex=1)


		# host
		a<-0; y.s<-1.01; y.t<-1.04
		for (i in 2:length(host)){
			if (host[i-1] != host[i] || loc[i-1] != loc[i] || sp[i-1] != sp[i]){
				segments(i-1,0,i-1,1, col="grey",lwd=1.5, lty=3)
				segments(a+0.5,y.s,i-1.5,y.s)
				text(a+((i-1)-a)/2, y.t, host[i-1], cex=0.8, adj=0.5)
				a<-i-1
				if (i==length(host)){
					segments(a+0.5,y.s,i-0.5,y.s)
					text(a+(i-a)/2, y.t, host[i], cex=0.8, adj=0.5)
				}
			}
			else if (i==length(host)){
				segments(a+0.5,y.s,i-0.5,y.s)
				text(a+(i-a)/2, y.t, host[i-1], cex=0.8, adj=0.5)
			}
		}

		# locality
		a<-0; y.s<-1.07; y.t<-1.10
		for (i in 2:length(loc)){
			if (loc[i-1] != loc[i] || sp[i-1] != sp[i]){
				segments(i-1,0,i-1,1, col="grey",lwd=1.5)
				segments(a+0.5,y.s,i-1.5,y.s)
				text(a+((i-1)-a)/2, y.t, loc[i-1], cex=0.8, adj=0.5)
				a<-i-1
				if (i==length(loc)){
					segments(a+0.5,y.s,i-0.5,y.s)
					text(a+(i-a)/2, y.t, loc[i], cex=0.8, adj=0.5)
				}
			}
			else if (i==length(loc)){
				segments(a+0.5,y.s,i-0.5,y.s)
				text(a+(i-a)/2, y.t, loc[i-1], cex=0.8, adj=0.5)
			}
		}

		# species
		a<-0; y.s<-1.13; y.t<-1.16
		for (i in 2:length(sp)){
			if (sp[i-1] != sp[i]){
				segments(i-1,0,i-1,1, col="black",lwd=1.5)
				segments(a+0.5,y.s,i-1.5,y.s)
				text(a+((i-1)-a)/2, y.t, sp[i-1], cex=0.8, adj=0.5)
				a<-i-1
				if (i==length(sp)){
					segments(a+0.5,y.s,i-0.5,y.s)
					text(a+(i-a)/2, y.t, sp[i], cex=0.8, adj=0.5)
				}
			}
			else if (i==length(sp)){
				segments(a+0.5,y.s,i-0.5,y.s)
				text(a+(i-a)/2, y.t, sp[i-1], cex=0.8, adj=0.5)
			}
		}

		dev.off()`);

}

# ==============================================================================
# Plot DIC plot
# ==============================================================================
sub plot_DIC{
	my $input=shift;
	my $output=shift;

	$R->run(qq`
		dic<-read.table("$input", header=T, sep=',')
		pdf(file="$output", width=1000/72, height=500/72)
		par (mfcol=c(2,1))
		plot(dic\$K, dic\$DIC, ylab="DIC", xlab="K", cex=1.25, type = "l")
		points(dic\$K, dic\$DIC, pch = 21, bg ="gray", col = "black", cex=1.25)
		plot(dic\$K, dic\$DIC, ylab="DIC", xlab="K", cex=1.25, type = "l", ylim=quantile(dic\$DIC,c(0,0.70)))
		points(dic\$K, dic\$DIC, pch = 21, bg ="gray", col = "black", cex=1.25)
		dev.off()`);
}

# ==============================================================================

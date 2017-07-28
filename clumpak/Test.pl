use strict;
use warnings;
use File::Slurp;
use File::Path qw(make_path remove_tree);
use Archive::Extract;
use Archive::Zip;
use File::Basename; 

my $cmd = "ssh lecs lss 2>&1";
my $cmd2 = "ssh bla ls 2>&1";

print "ssh to lecs\n";
my $output = `$cmd`;
my $exitVal = $? >> 8;

print "exit code: $exitVal\n";
print "output: $output\n";


print "ssh to bla\n";
$output = `$cmd2`;
$exitVal = $? >> 8;

print "exit code: $exitVal\n";
print "output: $output\n";

print "done\n";
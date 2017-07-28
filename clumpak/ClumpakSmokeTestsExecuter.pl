#!/usr/bin/perl -w

use strict;
use warnings;

use SmokeTestsMethods;

my $smokeTestsDir = "/bioseq/data/results/CLUMPAK/SmokeTests";


print "Creating .sh file\n";
my $qsub_script = "$smokeTestsDir/qsub.sh";
open (QSUB_SH,">$qsub_script");
  
print QSUB_SH "#!/bin/tcsh\n";
print QSUB_SH '#$ -N CLUMPAK_SMOKE_TESTS', "\n";
print QSUB_SH '#$ -S /bin/tcsh', "\n";
print QSUB_SH '#$ -cwd', "\n";
print QSUB_SH '#$ -l bioseq', "\n";
print QSUB_SH '#$ -e ', "$smokeTestsDir", '/$JOB_NAME.$JOB_ID.ER', "\n";
print QSUB_SH '#$ -o ', "$smokeTestsDir", '/$JOB_NAME.$JOB_ID.OU', "\n";
print QSUB_SH "cd /bioseq/CLUMPAK/;perl ClumpakSmokeTests.pl\n";

# this will send results ready email to user 
#my $cmdEmail = "cd /bioseq/$serverName/;perl sendLastEmail.pl --toEmail $email_to_address --id $jobId;";
#print QSUB_SH "$cmdEmail\n";

close (QSUB_SH);



my $cmd =  "ssh lecs2 qsub $qsub_script 2>&1";

print "qsub command: $cmd\n";

my $output = `$cmd`;
my $exitVal = $? >> 8;

if ($exitVal != 0) {
	print "error occurred. output:\n$output\n";
	
	my $testsFlagVal = 0;
	my $hasFlagChanged = &ChangeTestsSuccessFlag($testsFlagVal);
	
	my $mailMsgBody = "Error on clumpak smoke tests.\nFailed to submit tests job from ibis to lecs2 queue.\nFailed command: $cmd\nCommand ouptut: $output\n";
	
	if ($hasFlagChanged) {
		print "Smoked tests flag has changed. Disabled job submission";
		
		$mailMsgBody = $mailMsgBody."Changed Smoke Tests flag value. Disabled job submission.\n";	
	}
	
	&SendSmokeTestResultToEmail("CLUMPAK smoke test error", $mailMsgBody);
	
	
	# further erorr handling. change htmls, cgis, etc etc
}
else {
	print "Finished submitting smoke tests\n";
}








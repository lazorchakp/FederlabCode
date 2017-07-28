package SmokeTestsMethods;

use strict;
use warnings;
use CLUMPAK_CONSTS_and_Functions;

use vars qw(@ISA @EXPORT);

require Exporter;

@ISA = qw(Exporter);

@EXPORT = qw(SendSmokeTestResultToEmail ChangeTestsSuccessFlag ReadTestsSuccessFlag);

use constant SMOKE_TESTS_FLAG_FILENAME => "/bioseq/CLUMPAK/smokeTestsSuccessFlag.dat";



sub SendSmokeTestResultToEmail {
	my ($msgSubject, $msgBody) = @_;
	
	my $mailAddress = 'evolseq@post.tau.ac.il';

#	my $subject = 'CLUMPAK smoke test error';
#	my $message = "Error on clumpak smoke tests.\n$msgBody\n";
 
	open(MAIL, "|/usr/sbin/sendmail -t");
 
	# Email Header
	print MAIL "To: $mailAddress\n";
	print MAIL "From: $mailAddress\n";
	print MAIL "Subject: $msgSubject\n\n";
	
	# Email Body
	print MAIL $msgBody;

	close(MAIL);
}

sub ChangeTestsSuccessFlag {
	my ($flagValue) = @_;
	
	my $curFlagValue = &ReadFromFile(SMOKE_TESTS_FLAG_FILENAME, 0);
	
	my $changedFlagValue = 0;
	
	if ($curFlagValue != $flagValue) {
		my $shouldOverwrite = 1;
		&WriteToFile(SMOKE_TESTS_FLAG_FILENAME, $flagValue, $shouldOverwrite);
		$changedFlagValue = 1;
	}
	
	return $changedFlagValue;
}

sub ReadTestsSuccessFlag {
	my $curFlagValue = &ReadFromFile(SMOKE_TESTS_FLAG_FILENAME, 0);
	
	return $curFlagValue;
}

1;



use strict;
use warnings;
use File::Path qw(make_path remove_tree); 
use CLUMPAK_CONSTS_and_Functions;
use SmokeTestsMethods;

my $destDir = "/bioseq/data/results/CLUMPAK/SmokeTests/results";
my $logsDir = "/bioseq/data/results/CLUMPAK/SmokeTests/logs";


#my $destDir = "/groups/itay_mayrose/jonathanm/SmokeTests/results";
#my $logsDir = "/groups/itay_mayrose/jonathanm/SmokeTests/logs";
#my $filesDir = "/groups/itay_mayrose/jonathanm/structure_zip_files/test_files";

my $filesDir = "/bioseq/data/results/CLUMPAK/SmokeTests/files";
#my $programDir ="/bioseq/CLUMPAK/";

my $mainPipelineFile = "clumpak_manyK_structure_toy_data.zip";
my $distructFile = "distructForManyKs_structure_toy_data.zip";
my $labelFile = "distruct_labels_toy_data.txt";

print "CLUMPAK Smoke Tests\n";

make_path($destDir);
make_path($logsDir);
remove_tree( $destDir, {keep_root => 1} );

print "Calling Tests..\n";

my $testsFlagVal = &CallCLUMPAKStructure($^T, $mainPipelineFile, $labelFile);

if ($testsFlagVal) {
	my $hasFlagChanged = &ChangeTestsSuccessFlag($testsFlagVal);

	my $mailMsgBody = "CLUMPAK smoke tests finished successfully.\n";
		
	if ($hasFlagChanged) {
		print "Smoked tests flag has changed. Enabled job submission";
		
		$mailMsgBody = $mailMsgBody."Changed Smoke Tests flag value. Enabled job submission.\n";	
	}
	
	&SendSmokeTestResultToEmail("CLUMPAK smoke test success", $mailMsgBody);
}

print "Finished runnig tests\n";
exit(0);




sub CallExternalScript {
	my ($cmd, $jobDir, $jobId) = @_;
	
	print "Command:\n$cmd\nExecuting...\n";
	my $output = `$cmd`;
	my $exitVal = $? >> 8;
	
	if ($exitVal) {
		my $errLog = "$jobDir/" . CLUMPAK_CONSTS_and_Functions::ERROR_STATUS_LOG_FILE;

		my $error = &ReadFromFile($errLog, "");
		print "error occurred:\n$error\n";
		
		if (index($error, "Error occurred running") != -1) {
			# error was in external program. running cmd again
			print "Error on external program. Running command again.\nCommand:\n$cmd\nExecuting...\n";
			remove_tree( $jobDir, {keep_root => 1} );
			
			my $output = `$cmd`;
			my $exitVal = $? >> 8;
			
			if ($exitVal) {
				my $error = &ReadFromFile($errLog, "");
				print "error occurred:\n$error\n";
				return &HandleError($jobDir, $jobId, $error);
				
			}
			else {
				print "command finished successfully.\n";
		
				my $testsFlagVal = 1;
				return $testsFlagVal;
			}
			
		}
		else {
			return &HandleError($jobDir, $jobId, $error);
		}
	}
	else {
		print "command finished successfully.\n";
		
		my $testsFlagVal = 1;
		return $testsFlagVal;
	}
}

sub HandleError {
	my ($jobDir, $jobId, $error) = @_;
	my $testsFlagVal = 0;
	my $hasFlagChanged = &ChangeTestsSuccessFlag($testsFlagVal);
	
	my $mailMsgBody = "Error on clumpak smoke tests.\nerror occurred on job $jobId.\nJob dir: $jobDir\nError: $error\n";
	
	if ($hasFlagChanged) {
		print "Smoked tests flag has changed. Disabled job submission";
		$mailMsgBody = $mailMsgBody."Changed Smoke Tests flag value. Disabled job submission.\n";	
	}
	
	print "sending error to mail..\n";
	&SendSmokeTestResultToEmail("CLUMPAK smoke test error", $mailMsgBody);
	
	return $testsFlagVal;
}


sub CallCLUMPAKStructure{
	my ($jobId, $zipFile, $labelFile) = @_;
	
	$jobId = "clumpak_stru_$jobId";
	my $jobDir = "$destDir/$jobId";
	
	my $cmd = "perl CLUMPAK.pl --id $jobId --dir $destDir/$jobId --file $filesDir/$zipFile";
	
	if (defined $labelFile) {
		$cmd = $cmd." --label $filesDir/$labelFile";
	}
	 
	$cmd = $cmd." > $logsDir/$jobId.log";
	
	&CallExternalScript($cmd, $jobDir, $jobId);
}

sub CallCLUMPAKAdmixture{
	my ($jobId, $zipFile, $indToPopFile) = @_;

	$jobId = "clumpak_admix_$jobId";
	my $jobDir = "$destDir/$jobId";
	
	
	my $cmd = "perl CLUMPAK.pl --id $jobId --dir $destDir/$jobId --file $filesDir/$zipFile --indtopop $filesDir/$indToPopFile --inputtype admixture";
	 
	$cmd = $cmd." > $logsDir/$jobId.log";
	
	&CallExternalScript($cmd, $jobDir);
	
}

sub CallBestKByEvanno{
	my ($jobId, $zipFile) = @_;

	$jobId = "bestK_$jobId";
	my $jobDir = "$destDir/$jobId";
	
	my $cmd = "perl BestKByEvanno.pl -i $jobId -d $destDir/$jobId -f $filesDir/$zipFile";

	 
	$cmd = $cmd." > $destDir/$jobId.log";
		
	&CallExternalScript($cmd, $jobDir);

}

sub CallCompareDifferentPrograms{
	my ($jobId, $firstZipType, $firstZipFile, $secondZipType, $secondZipFile, $labelFile, $indToPop) = @_;

	$jobId = "compare_$jobId";
	my $jobDir = "$destDir/$jobId";
	
	my $cmd = "perl CompareDifferentPrograms.pl --id $jobId --dir $destDir/$jobId --firstinputtype $firstZipType --firstfile $filesDir/$firstZipFile
					 --secondinputtype $secondZipType --secondfile $filesDir/$secondZipFile";
	
	if (defined $labelFile) {
		$cmd = $cmd." --label $filesDir/$labelFile";
	}
	
	if (defined $indToPop) {
		$cmd = $cmd." --indtopop $filesDir/$indToPop";
	}
	 
	$cmd = $cmd." > $destDir/$jobId.log";
		
	&CallExternalScript($cmd, $jobDir);

}

sub CalldistructForManyKsStructure{
	my ($jobId, $zipFile, $labelFile) = @_;
	
	$jobId = "distruct_stru_$jobId";
	my $jobDir = "$destDir/$jobId";

	my $cmd = "perl distructForManyKs.pl --id $jobId --dir $destDir/$jobId --file $filesDir/$zipFile";
	
	if (defined $labelFile) {
		$cmd = $cmd." --label $filesDir/$labelFile";
	}
	 
	$cmd = $cmd." > $destDir/$jobId.log";
		
	&CallExternalScript($cmd, $jobDir);

}

sub CalldistructForManyKsAdmixture{
	my ($jobId, $zipFile, $indToPop) = @_;

	$jobId = "distruct_admix_$jobId";
	my $jobDir = "$destDir/$jobId";
	
	my $cmd = "perl distructForManyKs.pl --id $jobId --dir $destDir/$jobId --file $filesDir/$zipFile --inputtype admixture --indtopop $filesDir/$indToPop";

	$cmd = $cmd." > $destDir/$jobId.log";
		
	&CallExternalScript($cmd, $jobDir);
}



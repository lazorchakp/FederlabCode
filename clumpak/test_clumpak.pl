use strict;
use warnings;
use File::Path qw(make_path remove_tree); 

my $destDir = "/groups/itay_mayrose/jonathanm/test_results";
my $filesDir = "/groups/itay_mayrose/jonathanm/structure_zip_files/test_files";
my $programDir ="/groups/itay_mayrose/jonathanm/StructureServer";

make_path($destDir);
remove_tree( $destDir, {keep_root => 1} );

######## CLUMPAK structure ##########
&CallCLUMPAKStructure("clumpak_singleK_structure_short_stru_Fig5_res_5", "clumpak_singleK_structure_short_stru_Fig5_res_5.zip", "distruct_labels_short_stru.txt");
&CallCLUMPAKStructure("clumpak_singleK_structure_short_stru_Fig5_res_5_no_label", "clumpak_singleK_structure_short_stru_Fig5_res_5.zip");
&CallCLUMPAKStructure("clumpak_manyK_structure_short_stru_Fig5", "clumpak_manyK_structure_short_stru_Fig5.zip", "distruct_labels_short_stru.txt");
&CallCLUMPAKStructure("clumpak_manyK_structure_short_stru_Fig5_no_label", "clumpak_manyK_structure_short_stru_Fig5.zip");
&CallCLUMPAKStructure("clumpak_manyK_structure_toy_data", "clumpak_manyK_structure_toy_data.zip", "distruct_labels_toy_data.txt");
&CallCLUMPAKStructure("clumpak_manyK_structure_toy_data_no_label", "clumpak_manyK_structure_toy_data.zip");

####### CLUMPAK admixture ##########
&CallCLUMPAKAdmixture("clumpak_manyK_admixture", "clumpak_manyK_admixture.zip", "admixture_lables.txt");
&CallCLUMPAKAdmixture("clumak_singleK_admixture_K5", "clumak_singleK_admixture_K5.zip", "admixture_lables.txt");

####### bestKByEvanno ##########
&CallBestKByEvanno("best_k_by_evanno_structure_toy_data", "best_k_by_evanno_structure_toy_data.zip");
&CallBestKByEvanno("best_k_by_evanno_structure_short_stru_Fig5", "best_k_by_evanno_structure_short_stru_Fig5.zip");

####### CompareDifferentPrograms structure ##########
&CallCompareDifferentProgramsStructure("compare_structure_short_stru", "compare_diff_programs_structure_first_short_stru_Fig5_res_5.zip", "compare_diff_programs_structure_second_short_stru_Fig5_res_5.zip", "distruct_labels_short_stru.txt");
&CallCompareDifferentProgramsStructure("compare_structure_short_stru_no_label", "compare_diff_programs_structure_first_short_stru_Fig5_res_5.zip", "compare_diff_programs_structure_second_short_stru_Fig5_res_5.zip");

####### CompareDifferentPrograms admixture ##########
&CallCompareDifferentProgramsAdmixture("compare_diff_programs_admixture", "compare_diff_programs_admixture_first_K5.zip", "compare_diff_programs_admixture_second_K5.zip", "admixture_lables");

####### CompareDifferentPrograms structure to admixture ##########
&CallCompareDifferentProgramsStructure("compare_structure_to_admixture", "compare_diff_programs_structure_runs_K5.zip", "compare_diff_programs_admixture_runs_K5.zip");

####### CalldistructForManyKs structure ##########
&CalldistructForManyKsStructure("distructForManyKs_structure_Fig5", "distructForManyKs_structure_Fig5.zip", "distruct_labels_short_stru.txt");
&CalldistructForManyKsStructure("distructForManyKs_structure_Fig5_no_label", "distructForManyKs_structure_Fig5.zip");
&CalldistructForManyKsStructure("distructForManyKs_structure_toy_data", "distructForManyKs_structure_toy_data.zip", "distruct_labels_toy_data.txt");
&CalldistructForManyKsStructure("distructForManyKs_structure_toy_data_no_label", "distructForManyKs_structure_toy_data.zip");


####### CalldistructForManyKs admixture ##########
&CalldistructForManyKsAdmixture("distructForManyKs_admixture", "distructForManyKs_admixture.zip", "admixture_lables.txt");


#########################################


sub CallCLUMPAKStructure{
	my ($jobId, $zipFile, $labelFile) = @_;

	my $cmd = "perl CLUMPAK.pl --id $jobId --dir $destDir/$jobId --file $filesDir/$zipFile";
	
	if (defined $labelFile) {
		$cmd = $cmd." --label $filesDir/$labelFile";
	}
	 
	$cmd = $cmd." > $destDir/$jobId.log";
	
	print "\n\nCLUMPAK command:\n";
	print $cmd, "\n"; 
	print `$cmd`;
}

sub CallCLUMPAKAdmixture{
	my ($jobId, $zipFile, $indToPopFile) = @_;

	my $cmd = "perl CLUMPAK.pl --id $jobId --dir $destDir/$jobId --file $filesDir/$zipFile --indtopop $filesDir/$indToPopFile --inputtype admixture";
	 
	$cmd = $cmd." > $destDir/$jobId.log";
	
	print "\n\nCLUMPAK command:\n";
	print $cmd, "\n"; 
	print `$cmd`;
}

sub CallBestKByEvanno{
	my ($jobId, $zipFile) = @_;

	my $cmd = "perl BestKByEvanno.pl -i $jobId -d $destDir/$jobId -f $filesDir/$zipFile";

	 
	$cmd = $cmd." > $destDir/$jobId.log";
		
	print "\n\nBestKByEvanno command:\n";
	print $cmd, "\n"; 
	print `$cmd`;
}

sub CallCompareDifferentProgramsStructure{
	my ($jobId, $firstZipFile, $secondZipFile, $labelFile) = @_;

	my $cmd = "perl CompareDifferentPrograms.pl --id $jobId --dir $destDir/$jobId --firstfile $filesDir/$firstZipFile --secondfile $filesDir/$secondZipFile";
	
	if (defined $labelFile) {
		$cmd = $cmd." --label $filesDir/$labelFile";
	}
	 
	$cmd = $cmd." > $destDir/$jobId.log";
		
	print "\n\nCompareDifferentPrograms command:\n";
	print $cmd, "\n"; 
	print `$cmd`;
}

sub CallCompareDifferentProgramsAdmixture{
	my ($jobId, $firstZipFile, $secondZipFile, $indToPop) = @_;

	my $cmd = "perl CompareDifferentPrograms.pl --id $jobId --dir $destDir/$jobId --firstfile $filesDir/$firstZipFile --secondfile $filesDir/$secondZipFile";
	$cmd = $cmd." --firstinputtype admixture --secondinputtype admixture --indtopop $filesDir/$indToPop";
	
	 
	$cmd = $cmd." > $destDir/$jobId.log";
		
	print "\n\nCompareDifferentPrograms command:\n";
	print $cmd, "\n"; 
	print `$cmd`;
}

sub CallCompareDifferentProgramsStructureToAdmixture{
	my ($jobId, $structureZipFile, $admixtureZipFile, $labelFile) = @_;

	my $cmd = "perl CompareDifferentPrograms.pl --id $jobId --dir $destDir/$jobId --firstfile $filesDir/$structureZipFile --secondfile $filesDir/$admixtureZipFile";
	$cmd = $cmd." --secondinputtype admixture";
	
	if (defined $labelFile) {
		$cmd = $cmd." --label $filesDir/$labelFile";
	}
	 
	$cmd = $cmd." > $destDir/$jobId.log";
		
	print "\n\nCompareDifferentPrograms command:\n";
	print $cmd, "\n"; 
	print `$cmd`;
}

sub CalldistructForManyKsStructure{
	my ($jobId, $zipFile, $labelFile) = @_;

	my $cmd = "perl distructForManyKs.pl --id $jobId --dir $destDir/$jobId --file $filesDir/$zipFile";
	
	if (defined $labelFile) {
		$cmd = $cmd." --label $filesDir/$labelFile";
	}
	 
	$cmd = $cmd." > $destDir/$jobId.log";
		
	print "\n\ndistructForManyKs command:\n";
	print $cmd, "\n"; 
	print `$cmd`;
}

sub CalldistructForManyKsAdmixture{
	my ($jobId, $zipFile, $indToPop) = @_;

	my $cmd = "perl distructForManyKs.pl --id $jobId --dir $destDir/$jobId --file $filesDir/$zipFile --inputtype admixture --indtopop $filesDir/$indToPop";

	$cmd = $cmd." > $destDir/$jobId.log";
		
	print "\n\ndistructForManyKs command:\n";
	print $cmd, "\n"; 
	print `$cmd`;
}

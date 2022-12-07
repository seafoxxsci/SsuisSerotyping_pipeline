#!/usr/bin/perl
use strict;
use warnings;
use File::Glob;
use POSIX;
use Getopt::Long;
use Cwd 'abs_path';

#: Make sure child dies if keyboard interrupts or kills the signal
$SIG{CHLD} = 'IGNORE';
$SIG{INT} = \&signal_handler;
$SIG{TERM} = \&signal_handler;

#First argument is the Fasta database file, second argument is the results from SRST2, should start in the same directory as the scores file.

#: Get path to script and user inputs
my $script_path = abs_path($0);
my $script_dir = dirname($script_path);
my $working_dir = getcwd();
my $fastq_directory;
my $scoresName;
my $fasta_input;
my $definitions_file;
my $gene_fasta;
my $MLST_input;
my $MLST_definitions_file;
my $recN_input;
my $virulence_input;
my $forward;
my $reverse;
my $SingleOrPaired;

my %opt=qw();
GetOptions(\%opt, 
	"help|h!", "fastq_directory:s", "scoreName:s", "serotype_db:s", "serotype_definitions:s", "cps2K:s", 
	"MLST_db:s", "MLST_definitions:s", "recN_db:s", "Virulence_db:s", "forward:s", "reverse:s", "ends:s");


if (exists($opt{"help"})) {
    print <<EOF;
**********************************************
Implementation of the Streptococcus suis serotyping pipeline:
**********************************************

perl Ssuis_serotypingPipeline.pl --fastq_directory /path/to/fastq/directory --scoreName Scores_output_name --serotype_db serotype.fasta -- serotype_definitions serotype_definitions.txt --cps2K cps2K.fasta


--fastq_directory   Path to directory containing paired-end or single-end fastq files. Must be full path to 
			directory. Please do not use '.' or '..' to declare path. 
			[default: your working directory]
					
--scoreName         Name of SRST2 results file. 
			[default: 'Results']

--serotype_db       Multifasta file containing the serotype database.
                    	[default: 'Ssuis_Serotyping.fasta' in the directory containing the script]
					
--serotype_definitions   Text file containing the definitions for the serotype database file.
                    	[default: 'Ssuis_Serotyping_Definitions.txt' in the directory containing the script]
					
--cps2K             Multifasta file containing the cpsH confirmation database.
                    	[default: 'Ssuis_cps2K.fasta' in the directory containing the script]
					
--MLST_db           Multifasta file containing the MLST database.
                    	[default: 'Streptococcus_suis.fasta' in the directory containing the script]
					
--MLST_definitions  Text file containing the definitions for the MLST database file.
                    	[default: 'ssuis.txt' in the directory containing the script]

--recN_db           Fasta file containing the recN species specfic gene.
                    	[default: 'recN_Full.fasta' in the directory containing the script]
					
--Virulence_db      Multifasta file containing the Virulence genes.
                  	[default: 'Virulence.fasta' in the directory containing the script]
					
--forward           Indicator delimiting the forward reads file for paired-end read fastq files. This 
			option is ignored if 'se' is selected.
			[default: '_R1']

--reverse	    Indicator delimiting the reverse reads file for paired-end read fastq files. This 
			option is ignored if 'se' is selected.
			[default: '_R2']

--ends		    Indicates whether the reads are paired-end 'pe' or single-end 'se' fastq files. Note: 
			We recommend using paired-end reads of at least 100nt in length and 30X coverage. We 
			have not tested the efficiency of this pipeline with reads shorter than 80nt.
			[default: 'pe']
EOF
}


#: Set values to defaults in an array
my %inputs = (
  "fastq_directory" => $working_dir,
  "serotype_db" => $script_dir."/Ssuis_Serotyping.fasta",
  "serotype_definitions" => $script_dir."/Ssuis_Serotyping_Definitions.txt",
  "cps2K" => $script_dir."/Ssuis_cps2K.fasta",
  "MLST_db" => $script_dir."/Streptococcus_suis.fasta",
  "MLST_definitions" => $script_dir."/ssuis.txt",
  "recN_db" => $script_dir."/recN_full.fasta",
  "Virulence_db" => $script_dir."/Virulence.fasta",
  "forward" => "_R1",
  "reverse" => "_R2",
  "ends" => "pe"
);
# Double-check that defaults or user-inputs are present for each flag; print an error message and exit if values are not provided
foreach my $key (keys %inputs) {
	if(exists($opt{$key})) {
    	$inputs{$key} = $opt{$key};
	} else {
	if ($key eq "fastq_directory") {
		$inputs{$key} = $workingdir;
		my $any_fastqs = glob("*.fastq");
		if (!defined($any_fastqs)) {
			print "Please provide the full path to the directory containing fastq files to use. See help 
			file for additional details [--help]."
			exit;
			}
		} else {
  			my $file = $inputs{$key};
			unless(-e $file) {
    			print "Please provide the $key file(s). See help file for additional details [--help].";
      			exit;
			}
		}
	}
	if ($key eq "ends") {
	my $SingleOrPaired = $inputs{$key};
	if(($SingleOrPaired ne "pe") and ($SingleOrPaired ne "se")){
		die "Ends must be either pe or se";
}


mkdir "recN" or die "Failed to create recN directory: $!";
chdir "recN" or die "Failed to change recN directory: $!";

#: Set up arguments for SRST2 based on whether input is paired-end reads or single-end reads
my $srst2_args = "--log --gene_db $recN_input --save_scores";
if ($SingleOrPaired eq "pe") {
  $srst2_args .= " --input_pe $fastq_directory/*.fastq --forward $forward --reverse $reverse";
} elsif ($SingleOrPaired eq "se") {
  $srst2_args .= " --input_se $fastq_directory/*.fastq";
} else { # Throw an error if input is neither paired-end 'pe' or single-end 'se'
  die "Invalid value for SingleOrPaired: input must be 'se' or 'pe'. See help file for additional details [--help].";
}
# Run SRST2 to check if the sequence file actually belongs to Streptococcus suis
system("srst2.py $srst2_args --output $scoresName\_recN");


#: Open the recN results file to verify that the SRST2 pipeline ran correctly
my $recNResults = glob("$scoresName\_recN__genes*.txt");
if (!defined($recNResults)) {
  die "Failed to find recN results file";
}
open my $recNs, "<$recNResults" or die "Can't open recN results file: $!";

#: Initialize arrays to store filenames of Streptococcus suis and non-Streptococcus suis samples
my @ssuis;
my @nonssuis;

#: Read the recN results file line-by-line:
my $recCount = 0;
while (my $recN = <$recNs>) {
	$recN =~ s/\r?\n//; # Remove any line breaks from the end of each line

	if ($recCount > 0) { 
		# Split lines on whitespace characters to get sample and gene names; skip the header
		my @info = split(/\s+/, $recN);
		next if $recN == 0;
		
		#: For all other lines, which contain sample and gene names (more than one field)
		if (scalar(@info) > 1) {
			# Check if the gene name starts with recN
			if(substr($info[1], 0, 4) eq "recN"){
				# Add both fastq files to the list of S. suis samples if paired-end
				if($SingleOrPaired eq "pe"){
					push(@ssuis, "$fastq_directory/$info[0]$forward.fastq");
					push(@ssuis, "$fastq_directory/$info[0]$reverse.fastq");
				} elsif($SingleOrPaired eq "se"){ # Add the fastq file if single-end
					push(@ssuis, "$fastq_directory/$info[0]*.fastq");
				}
			}
			# If the gene name doesn't start with recN, add sample to the non-S. suis list
			else{push(@nonssuis, $info[0])}; 
		}
	}
	$recCount++;
}

my $ss = join(' ', @ssuis);

#: Organize species outputs from the SRST2 analysis of recN gene hits in each sample
system("mkdir pileups sorted_bam scores")
system("mv *.pileup ./pileups && mv *.sorted.bam ./sorted_bam && mv *.scores ./scores");

close($recNs);
system("mv $recNResults $scoresName\_speciesConfirmation.txt");

# Check if there are any Streptococcus suis samples in the dataset
if(scalar(@ssuis) < 1){
	print "No Streptococcus suis samples are in the dataset";
	exit;
}

#: Create MLST, Virulence, and Serotype directories; organize system commands now that you have SRST2 ran
my @commands;
if($SingleOrPaired eq "pe"){
	@commands = ("srst2.py --input_pe $ss --forward $forward --reverse $reverse --output $scoresName\_MLST --log --mlst_db $MLST_input --mlst_definitions $MLST_definitions_file --save_scores", 
	"srst2.py --input_pe $ss --forward $forward --reverse $reverse --output $scoresName\_VirulenceFactors --log --gene_db $virulence_input --save_scores", 
	"srst2.py --input_pe $ss --forward $forward --reverse $reverse --output $scoresName --log --mlst_db $fasta_input --mlst_definitions $definitions_file --save_scores");
} elsif ($SingleOrPaired eq "se") {
	@commands = ("srst2.py --input_se $ss --output $scoresName\_MLST --log --mlst_db $MLST_input --mlst_definitions $MLST_definitions_file --save_scores", 
	"srst2.py --input_se $ss --output $scoresName\_VirulenceFactors --log --gene_db $virulence_input --save_scores", 
	"srst2.py --input_se $ss --output $scoresName --log --mlst_db $fasta_input --mlst_definitions $definitions_file --save_scores");
}

#: For each command in the @commands array, change into that directory and execute the respective command
my @dirs = ("MLST", "Virulence", "Serotype");
my $prev_pid;
chdir "..";
foreach my $i (0 .. $#commands) {
	mkdir $dirs[$i];
    	chdir $dirs[$i];

   	my $pid = fork();
    	if ($pid == 0) {
        	setpgrp;
		system($commmands[$i]);
		
		#: Files in the Serotype directory are handled slightly differently
		if ($i == 2) {
			my $ResultsName = glob("$scoresName\__mlst__*.txt");
			# Checks to see if the output file ResultsName has any results
			if ((-s $ResultsName) < 50){
				print "Can't run Serotyping!";
				&signal_handler;
			}
		} 
		exit 0;
    	}
	# If the previous command is still running, kill it; else run the next command
	if ($prev_pid) {
		sleep(1);
    		my $status = waitpid($prev_pid, WNOHANG);
		if ($status == 0) { 
			print "Unable to run $dirs[$i] analysis, killing this command..."
			kill "TERM", $prev_pid;
			waitpid($prev_pid, 0);
			exit;
		} else {
			print "Done with $dirs[$i] analysis, running the next command..."
		}
	}
	$prev_pid = $pid
}


#Read in results of SRST2
open my $srst2_results, "<$ResultsName" or die "Can't open Serotyping results file!";

mkdir "Pipeline";
open my $EndResults, ">", "$scoresName\_FinalSerotypingResults.txt" or die $!;

print $EndResults "Strain\tSerotype\n";

#READING IN SEROTYPE CALLING RESULTS#
my $resultCount = 0;
my @furtherAnalysis_1;
my @furtherAnalysis_2;
my @notTesting;
foreach my $result (<$srst2_results>){
	$result =~ s/\r?\n//;

	if($resultCount > 0){
		my @reading = split(/\t/, $result);
		#IF SEROTYPE IS CALLED AS A 1, PUT IN ARRAY TO CHECK IF IT IS A 1 OR A 14#
		if(($reading[1] eq "1") || ($reading[1] eq "1*") || ($reading[1] eq "1?") || ($reading[1] eq "1*?")){
			push(@furtherAnalysis_1, $reading[0]);
		}
		#IF SEROTYPE IS CALLED AS A 2, PUT IN ARRAY TP CHECK IF IT IS A 2 OR A 1/2#
		elsif(($reading[1] eq "2") || ($reading[1] eq "2*") || ($reading[1] eq "2?") || ($reading[1] eq "2*?")){
			push(@furtherAnalysis_2, $reading[0]);
		}
		#IF SEROTYPE IS NEITHER A 1 OR A 2, PRINT RESULTS TO FILE#
		else{
			push(@notTesting, ($reading[0] . "\t" . $reading[1]));
			print $EndResults "$reading[0]\t$reading[1]\n";
		}
	}
	$resultCount++;
}

close($srst2_results);
system("mv $ResultsName $scoresName\_InitialCapsuleResults.txt");

chdir "Pipeline";

#CREATE BOWTIE AND SAMTOOLS INDEX FILES FOR CPSK FASTA FILE#
system("bowtie2-build $gene_fasta $gene_fasta");
system("samtools faidx $gene_fasta");

#FOR ALL OF THE 1s ALIGN TO CPSK FASTA TO SEE IF IT IS A 1 OR A 14#
foreach my $one (@furtherAnalysis_1){
	#Run Bowtie2 against cpsK
	if($SingleOrPaired eq "pe"){
		system("bowtie2 -1 $fastq_directory/$one$forward.fastq $fastq_directory/$one$reverse.fastq -S $one\_vs_cpsK.sam --very-sensitive-local --no-unal -a -x $gene_fasta");
	}
	elsif($SingleOrPaired eq "se"){
		system("bowtie2 -U $fastq_directory/$one*.fastq -S $one\_vs_cpsK.sam --very-sensitive-local --no-unal -a -x $gene_fasta");
	}
	#Run Samtools to call SNPs
	system("samtools view -b -o $one\_vs_cpsK.bam -q 1 -S $one\_vs_cpsK.sam");
	system("samtools sort $one\_vs_cpsK.bam $one\_vs_cpsK.sorted");
	system("samtools mpileup -u -L 1000 -f $gene_fasta -Q 20 -q 1 $one\_vs_cpsK.sorted.bam > $one\_vs_cpsK.pileup");
	system("bcftools view -vcg $one\_vs_cpsK.pileup > $one\_vs_cpsK.raw.vcf");
	system("vcfutils.pl varFilter -Q 20 $one\_vs_cpsK.raw.vcf > $one\_vs_cpsK.vcf");

	open my $vcf, "<$one\_vs_cpsK.vcf" or die "Can't find input file!";
	open my $vcf_out, ">", "$one\_vs_cpsK.snps" or die $!;

	#Re-write SNP file to a more readable format
	foreach my $vcf_line (<$vcf>){
		if(substr($vcf_line, 0, 2) ne "##"){
			print $vcf_out $vcf_line;
		}
	}

	open my $SNP_out, ">", "$one\_SNPeffect.txt" or die $!;

	#Check AminoAcid changes caused by SNPs
	system("perl $script_dir/SNP_AminoAcidChange.pl $one\_vs_cpsK.snps $gene_fasta > $one\_SNPeffect.txt");

	#READ IN SNP EFFFECT FILE TO CHECK IF THERE IS A TRP OR CYS AT AMINO ACID POSITION 161#
	open my $SNP_read, "<$one\_SNPeffect.txt" or die "Can't open input file!";

	my @SNP_array;
	my $SNP_line = 0;
	my $Target = 0;
	my $AminoAcid;
	foreach my $snp(<$SNP_read>){
		if($SNP_line > 0){
		        @SNP_array = split(/\t/, $snp);
			
			if($SNP_array[5] == 161){
				$Target = 1;
				$AminoAcid = $SNP_array[7];
			}
		}
		$SNP_line++;
	}

	#IF THERE WAS A SNP AT POSITION 161, CHECK IF IT IS A C, PRINT TO FILE#	
	if($Target == 1){
		if($AminoAcid eq "C"){
			print $EndResults "$one\t1\n";
			push(@notTesting, ($one . "\t1"));
		}
		else{
			print $EndResults "$one\t14\n";
			push(@notTesting, ($one . "\t14"));
		}
	}
	else{
		print $EndResults "$one\t14\n";
		push(@notTesting, ($one . "\t14"));
	}
	
}

#FOR ALL OF THE 2s ALIGN TO CPSK FASTA TO SEE IF IT IS A 2 OR A 1/2#
foreach my $two (@furtherAnalysis_2){
	#Run Bowtie against cpsK
	if($SingleOrPaired eq "pe"){
		system("bowtie2 -1 $fastq_directory/$two$forward.fastq $fastq_directory/$two$reverse.fastq -S $two\_vs_cpsK.sam --very-sensitive-local --no-unal -a -x $gene_fasta");
	}
	elsif($SingleOrPaired eq "se"){
		system("bowtie2 -U $fastq_directory/$two*.fastq -S $two\_vs_cpsK.sam --very-sensitive-local --no-unal -a -x $gene_fasta");		
	}

	#Run Samtools to call SNPs
	system("samtools view -b -o $two\_vs_cpsK.bam -q 1 -S $two\_vs_cpsK.sam");
	system("samtools sort $two\_vs_cpsK.bam $two\_vs_cpsK.sorted");
	system("samtools mpileup -u -L 1000 -f $gene_fasta -Q 20 -q 1 $two\_vs_cpsK.sorted.bam > $two\_vs_cpsK.pileup");
	system("bcftools view -vcg $two\_vs_cpsK.pileup > $two\_vs_cpsK.raw.vcf");
	system("vcfutils.pl varFilter -Q 20 $two\_vs_cpsK.raw.vcf > $two\_vs_cpsK.vcf");

	open my $vcf, "<$two\_vs_cpsK.vcf" or die "Can't find input file!";
	open my $vcf_out, ">", "$two\_vs_cpsK.snps" or die $!;

	#Re-write SNP file to a more readable format
	foreach my $vcf_line (<$vcf>){
		if(substr($vcf_line, 0, 2) ne "##"){
			print $vcf_out $vcf_line;
		}
	}

	open my $SNP_out, ">", "$two\_SNPeffect.txt" or die $!;

	#Check AminoAcid changes caused by SNPs
	system("perl $script_dir/SNP_AminoAcidChange.pl $two\_vs_cpsK.snps $gene_fasta > $two\_SNPeffect.txt");

	#READ IN SNP EFFFECT FILE TO CHECK IF THERE IS A TRP OR CYS AT AMINO ACID POSITION 161#
	open my $SNP_read, "<$two\_SNPeffect.txt" or die "Can't open input file!";

	my @SNP_array;
	my $SNP_line = 0;
	my $Target = 0;
	my $AminoAcid;
	foreach my $snp(<$SNP_read>){
		if($SNP_line > 0){
		        @SNP_array = split(/\t/, $snp);
			if($SNP_array[5] == 161){
				$Target = 1;
				$AminoAcid = $SNP_array[7];
			}
		}
		$SNP_line++;
	}

	#IF THERE WAS A SNP AT POSITION 161, CHECK IF IT IS A C, PRINT TO FILE#		
	if($Target == 1){
		if($AminoAcid eq "C"){
			print $EndResults "$two\t1/2\n";
			push(@notTesting, ($two . "\t1/2"));
		}
		else{
			print $EndResults "$two\t2\n";
			push(@notTesting, ($two . "\t2"));
		}
	}
	else{
		print $EndResults "$two\t2\n";
		push(@notTesting, ($two . "\t2"));
	}
	
}

#SORT OUTPUT FILES#
system("mkdir sam");
system("mkdir unsorted_bam");
system("mkdir sorted_bam");
system("mkdir pileup");
system("mkdir raw_vcf");
system("mkdir filtered_vcf");
system("mkdir snps");
system("mkdir snp_effect");
system("mv *.sam ./sam");
system("mv *.sorted.bam ./sorted_bam");
system("mv *.bam ./unsorted_bam");
system("mv *.pileup ./pileup");
system("mv *.raw.vcf ./raw_vcf");
system("mv *.vcf ./filtered_vcf");
system("mv *.snps ./snps");
system("mv *SNPeffect.txt ./snp_effect");

#WAIT FOR MLST TO FINISH BEFORE MOVING ON TO ORGANIZATION STEP
wait();

chdir "../../MLST";

#Organize SRST2 output
system("mkdir pileups");
system("mkdir sorted_bam");
system("mkdir scores");
system("mv *.pileup ./pileups");
system("mv *.sorted.bam ./sorted_bam");
system("mv *.scores ./scores");

my $mlst_Name = glob("$scoresName\_MLST__mlst__*.txt");

#Read in results of SRST2
open my $mlst_results, "<$mlst_Name" or die "Can't open MLST results file!";

#READING IN SEROTYPE CALLING RESULTS#
my $mlstresultCount = 0;
my @mlst;
foreach my $stresult (<$mlst_results>){
	$stresult =~ s/\r?\n//;

	if($mlstresultCount > 0){
		my @reading = split(/\t/, $stresult);
		push(@mlst, \@reading);
	}
	$mlstresultCount++;
}

#Change name of MLST output
close $mlst_results;
system("mv $mlst_Name $scoresName\_MLSTResults.txt");

chdir "../Virulence";

#Organize Virulence output
system("mkdir pileups");
system("mkdir sorted_bam");
system("mkdir scores");
system("mv *.pileup ./pileups");
system("mv *.sorted.bam ./sorted_bam");
system("mv *.scores ./scores");

my $vir_name = glob("$scoresName\_VirulenceFactors__genes*.txt");

#Read in Virulence factor gene results
open my $vir_results, "<$vir_name" or die "Can't open Virulence results files!";

#READING IN VIRULENCE RESULTS
my $virulenceResultsCount = 0;
my @virulence;
my @virHeader;
my $virNames;
my @virRead;
my @virStrain;
my @virWrite;
foreach my $virResult (<$vir_results>){
	$virResult =~ s/\r?\n//;

	if($virulenceResultsCount == 0){
		@virHeader = split(/\t/, $virResult);
		shift(@virHeader);
		$virNames = join("\t", @virHeader);
	}
	else{
		@virRead = split(/\t/, $virResult);
		push(@virStrain, $virRead[0]);
		shift(@virRead);
		push(@virWrite, join("\t", @virRead));
	}

	$virulenceResultsCount++;
}

my $factorNum = scalar(split(/\t/, $virNames));

close($vir_results);
system("mv $vir_name $scoresName\_VirulenceFactorResults.txt");

chdir "..";

#OPEN OUTPUT FILE
open my $combined_results, ">", "$scoresName\_FinalResults.txt" or die $!;
print $combined_results "Strain\tSerotype\tST\taroA\tcpn60\tdpr\tgki\tmutS\trecA\tthrA\trecN\t$virNames\n";

#NEW READ IN STRAINS
my $strainNum = 0;
foreach my $readings (@virStrain) {
	print $combined_results $readings;

	#FIND SEROTYPE RESULT MATCHING VIRULENCE STRAIN NAME
	my $seroMatch = 0;
	my $SeroCount = 0;
	my @sero_split;

	while(($seroMatch == 0) && ($SeroCount < scalar(@notTesting))){
		@sero_split = split(/\t/, $notTesting[$SeroCount]);
		if($sero_split[0] eq $readings){
			$seroMatch = 1;
			print $combined_results "\t$sero_split[1]";
		}
		else{$SeroCount++};
	}
	if($seroMatch == 0){print $combined_results "\t-"};

	#FIND MLST RESULT MATCHING VIRULENCE STRAIN NAME
	my $mlst_match = 0;
	my $mlst_count = 0;
	my $mlst_scalar = scalar(@mlst);

	while(($mlst_match == 0) && ($mlst_count < $mlst_scalar)){
		if($mlst[$mlst_count][0] eq $readings){
			$mlst_match = 1;
			print $combined_results "\t$mlst[$mlst_count][1]\t$mlst[$mlst_count][2]\t$mlst[$mlst_count][3]\t$mlst[$mlst_count][4]\t$mlst[$mlst_count][5]\t$mlst[$mlst_count][6]\t$mlst[$mlst_count][7]\t$mlst[$mlst_count][8]";
		}
		else{$mlst_count++};
	}
	if($mlst_match == 0){print $combined_results "\tNF\tNF\tNF\tNF\tNF\tNF\tNF\tNF"};

	#PRINT recN+ and Virulence results
	if(substr($virWrite[$strainNum], 0, 1) eq "?"){
		my $vir_rep = "\t-" x $factorNum;
		print $combined_results "\t+$vir_rep\n";
	}
	else{
		print $combined_results "\t+\t$virWrite[$strainNum]\n";
	}

	$strainNum++;
}

if(scalar(@nonssuis) > 0){
	my $vir_rep = "\tN/A" x $factorNum;
	foreach my $non (@nonssuis){
		print $combined_results "$non\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\t-$vir_rep\n";
	}
}

#MAKES SURE CHILD DIES IF PARENT IS INTERUPTED
sub signal_handler {
	my $exists_kill = kill 0, $pid1;
	if ($exists_kill == 1){
		kill -9, $pid1;
	}
	die;
};

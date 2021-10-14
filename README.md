# Running Ssuis_serotypingPipeline.pl

Dependencies for this pipeline are the same as those required for [SRST2](https://github.com/katholt/srst2). 
As well, the provided SNP_AminoAcidChange.pl scripts must be located in the same directory as Ssuis_serotypingPipeline.pl.

## Installation
This installation guide assumes that you have python and conda already installed on your system.

1. Create a new environment and installation tools
```
conda create -n ssuis_sero
conda activate ssuis_sero
conda install git pip
```

2. Install dependencies for SRST2
```
conda install -c conda-forge scipy
conda install -c bioconda bowtie2=2.2.4
conda install -c bioconda samtools=0.1.18-11
git clone https://github.com/katholt/srst2
pip install srst2/
```

4. Test the SRST2 installation
```
srst2 --version
getmlst.py -h
slurm_srst2.py -h
```

5. Install ssuis_sero directory from Github:
```
git clone https://github.com/seafoxxsci/SsuisSerotyping_pipeline
```

## Main executable and arguments for ssuis_sero
The execution line for the pipeline is as follows:
```
perl Ssuis_serotypingPipeline.pl --fastq_directory /path/to/fastq/directory --scoreName Scores_output_name
```

where Ssuis_serotypingPipeline.pl is the provided pipeline script, /path/to/fastq/directory is the path to the directory of the fastq files to be analyzed (Note: must be a full directory path, all fastq files in the folder will be analyzed). The pipeline runs with either single or paired-end reads, but looks for paired-end reads by default. There are optional commands --forward and --reverse used to indicate the names of the forward and reverse pairs, otherwise the program will assume pairs are named as _R1 and _R2 by default.

All parameters are optional.  If the fastq directory is not provided then the current working directory will be used.  If the current working directory does not contain fastq files, then the user will be prompted to provide an appropriate working directory.  If the database and fasta files are not provided, then they will be looked for in the directory containing the Ssuis_serotypingPipeline.pl script.  If the directory containing the script does not contain the database and fasta files, then the user will be prompted to provide these files.  To access the help file, the user can enter Ssuis_serotypingPipeline.pl --help.

Optional parameteres:
--fastq_directory	Path to directory containing paired-end fastq files. Must be full path to directory, please do not use '.' or '..' to declare path. If no path is given, the current working directory is used.

| Argument | Explanation |
|:-----|:------:|
|`--scoreName`            |Name of SRST2 results file [optional: default name 'Results']|
|`--serotype_db`          |Multifasta file containing the serotype database (Ssuis_Serotyping.fasta) [If none is provided, Ssuis_Serotyping.fasta is looked for in the directory containing the script].|
|`--serotype_definitions` |Text file containing the definitions for the serotype database file (Ssuis_Serotyping_Definitions.txt) [If none is provided, Ssuis_Serotyping_Definitions.txt is looked for in the directory containing the script].|
|`--cps2K`                |Multifasta file containing the cpsH confirmation database (Ssuis_cps2K.fasta) [If none is provided, Ssuis_cps2K.fasta is looked for in the directory containing the script].|
|`--MLST_db`              |Multifasta file containing the MLST database (Streptococcus_suis.fasta) [If none is provided, Streptococcus_suis.fasta is looked for in the directory containing the script].|
|`--MLST_definitions`     |Text file containing the definitions for the MLST database file (ssuis.txt) [If none is provided, ssuis.txt is looked for in the directory containing the script].|
|`--recN_db`              |Fasta file containing the recN species specfic gene (recN_Full.fasta) [If none is provided, recN_full.fasta is looked for in the directory containing the script].|
|`Virulence_db`		        |Multifasta file containing the Virulence genes (Virulence.fasta) [If none is provided, Virulence.fasta is looked for in the directory containing the script].|
|`--forward`              |Indicator delimiting the forward reads file for paired read fastq files. This option is ignored if single-end reads is selected [optional: default '_R1'].|
|`--reverse`              |Indicator delimiting the reverse reads file for paired read fastq files. This option is ignored if single-end reads is selected [optional: default '_R2'].|
|`--ends`			            |Indicates whether the reads are paired-end (pe) or single-end (se) [optional: default 'pe'].|

The pipeline checks for the presence of the gene recN first, and all strains containing recN move on to the next steps of the SRST2 pipeline: Serotyping, MLST, and virulence genes. Once SRST2 has run, the output is moved into separate directories.  

For all strains with the putative serotype 1 or 2, reads are re-mapped to the gene cps2K.  Based on SNPs in the cps2K, putative serotype 1s are either called as 1 or 14 and putative serotype 2s are either called as 2 or 1/2. 

Output from this pipeline analysis is contained in four directories:
- “recN”
- "Serotyping"
- "MLST"
- "Virulence." 

These directories contain the samtools pileups, sorted bam files, and score files, as well as the original results of these analyses.  The Serotyping directory also contains the directory pipeline, containing the files for serotype 1, 14, 2, and 1/2 confirmation.  Final results are located in the 'Results'_FinalResults.txt file.

Please note: We recommend using paired end reads of at least 100nt in length and at least 30X coverage. We have not tested the efficiency of the pipeline with reads shorter than 80nt.
The MLST database provided was downloaded from PubMLST (https://pubmlst.org/bigsdb?db=pubmlst_ssuis_seqdef) on September 29, 2021.  We recommend updating this database from time to time to guarantee accurate MLST assignment.

Provided files: Ssuis_serotypingPipeline.pl, recN_full.fasta, Ssuis_Serotyping.fasta, Ssuis_Serotyping_Definitions.txt, Ssuis_cps2K.fasta, Streptococcus suis.fasta, ssuis.txt, Virulence.fasta, SNP_AminoAcidChange.pl

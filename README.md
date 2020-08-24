# RNA-seq TCR/BCR Extraction Pipeline - Olopade Lab

## Introduction:

This pipeline is intended to automate the T-cell and B-cell CDR3 sequence extraction from paired-end bulk RNA sequencing. analysis around the breakpoints of structural variants (SVs). The program first extracts SV breakpoints from the input bedpe files: for each SV breakpoint, the program first checks whether it falls within the mappable regions of the chromosome, and then generates random positions within that same chromosome's mappable regions. The pipeline then extracts the 100 bp length sequence from either side of the breakpoint from the reference genome for each real and random SV breakpoint. The program then runs both the AME and FIMO algorithm to determine whether a specific sequence motif is enriched in the real SVs breakpoint sites versus the randomly generated ones.

This pipeline is intended to automate the use of TRUST for neoantigen analysis. The pipeline may also call upon the programs STAR (if alignment or genome indices generation is required) and Samtools (if an index file for genome FASTA file or the alignment BAM file has not already been generated).  
  
The Python script should be run as follows:  
**`$ python TRUST_script.py <options>`**  

## Requirements:

Python 3: https://www.python.org/downloads/release/python-373/  
TRUST: https://bitbucket.org/liulab/trust/src/master/  
STAR: https://github.com/alexdobin/STAR  
Samtools: http://www.htslib.org/  

## Options:  
**`-h` or `--help` :** displays the list of options and exits the program  
**`-d` or `--g_dir <path/to/directory>` :** specify the hg19 genome directory path (REQUIRED)  
**`-f` or `--g_fa <path/to/file>` :** specify the hg19 genome FASTA file path for STAR genomeGenerate (this will automatically run STAR genomeGenerate, which also requires -g/--g_gtf) *[default: the program will skip the STAR genomeGenerate step]*  
**`-g` or `--g_gtf <path/to/file>` :** specify the hg19 genome GTF file for STAR genomeGenerate (this will automatically run STAR genomeGenerate, which also requires -f/--g_fa) *[default: the program will skip the STAR genomeGenerate step]*  
**`-r` or `--r1r2_dir <path/to/directory>` :** specify the directory that contains the R1 and R2 FASTQ files for the sample (this will automatically run STAR alignment for the FASTQ files) *[default: the program will skip the STAR alignment step]*  
**`-a` or `--align_dir <path/to/directory>` :** specify the directory of the alignment BAM file (REQUIRED; if -r/--r1r2_dir is indicated, this directory will be used for the output of the STAR alignment)  
**`-o` or `--out <path/to/directory>` :** specify the final output directory for TRUST (REQUIRED)  
**`-p` or `--param <list parameters>` :** specify additional parameters for TRUST (the only allowed options are "B" and "BL" and "none" separated by commas, e.g. "none,BL") *[default: the program will run with no additional paramaters, e.g. "none"]*  
**`-s` or `--log_suffix <suffix>` :** specify the suffix to be added to the log file name, e.g. "TRUST_pipeline_log&lt;suffix>.txt" *[default: no suffix is added]*  
**`-l` or `--loop`:** runs the program in loop mode: the program will run for all subdirectories within the specified FASTQ and alignment BAM directories *[default: the program will run for only the specified directory]*  

## Input:

**_The paths for all files and directories in the command line should always be absolute paths._**  
The genome directory should contain the genome FASTA file as well as the GTF file of the hg19 genome. If the user is not running STAR genomeGenerate, the program will assume those index files are already generated and located in the genome directory. If there is no index file named &lt;genome FASTA file name>.fai, the program will automatically generate one.  
The directory containing the R1 and R2 FASTQ files should only contain the files specific to that sample. It should be named exactly like the sample name. It is assumed those files follow the typical naming conventions for FASTQ (https://support.illumina.com/content/dam/illumina-support/help/BaseSpaceHelp_v2/Content/Vault/Informatics/Sequencing_Analysis/BS/swSEQ_mBS_FASTQFiles.htm).  
The alignment directory should also be named exactly like the sample name. The alignment BAM file should be located within that directory, with the name &lt;sample name>Aligned.sortedByCoord.out.bam (this will be automatically generated if the user choses to run the STAR alignment). If there is no index file named &lt;alignment BAM file>.bai, the program will automatically generate one. If the alignment directory already exists, and the user is running the STAR alignment, the alignment directory and all of its content will be deleted, and a new one will be created.  
If the program is run in loop mode, all the alignment directories and the FASTQ directories for a specific sample should be located in a directory containing all of the samples of that file type, e. g. "fastq_directory/&lt;FASTQ directory for specific sample>/" and "alignment_directory/&lt;alignment directory for specific sample>/". The input file in the command line should be that directory containing all the subdirectories (continuing with the previous example, this would be: "fastq_directory/" and "alignment_directory/").  
The output directory should be a directory which will contain all the parameter-specific output directories that the program will generate (for all the specific samples, if the program is running in loop mode).  
The program will also generate a log file in the output directory named "TRUST_pipeline_log&lt;suffix>.txt" (or append to the file if it already exists). This log file will keep track of which samples have been aligned, which index files have been generated and for which samples TRUST has been run.

## Ouput:

If the user is running STAR genomeGenerate, all the index files will be generated within the genome input directory indicated in the command line.  
If the user is running the STAR alignment without the loop mode, the alignment files will be stored in the alignment directory specified by the user. If the user is running the alignment in loop mode, then the program will create a directory for each sample (i.e. each sample directory within the input FASTQ directory) within the indicated alignment directory, and output the result for a specific sample within that subdirectory. The subdirectory name will be the same as the FASTQ subdirectory name, which should be the sample name.  
For each parameter indicated in the command line (i.e. "none", "B" and/or "BL"), the program will output a directory with the corresponding suffix, i.e. "&lt;sample name>", "&lt;sample name>_B" and "&lt;sample name>_B_L". If the program is running in loop mode, each sample processed will have its own set of subdirectories within the main output directory.


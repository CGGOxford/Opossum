# Opossum README file


## 1. Requirements


Opossum requires Python 2.7 or greater. It has not been tested with the Python 3.X series.

It also requires the following Python packages: pysam v0.10.0, itertools, argparse, os and sys. 



## 2. Running Opossum


Opossum is a tool to pre-process RNA-seq reads prior to variant calling. Currently, variant callers do not generally behave well with reads encompassing intronic regions. Opossum splits the reads to get rid of the intronic parts. At the same time, it performs several quality control measures such as discarding secondary alignments, poorly mapped reads and read pairs where the two reads are pointing outwards or in the same direction, or that have been mapped in different chromosomes. Opossum also discards duplicate reads based on the start and end positions of the read pair. Furthermore, it merges overlapping paired-end reads.

Opossum has been designed to work specifically with Platypus, but it can be used equally well with other variant callers such as GATK.


Opossum can be run as follows:

    python Opossum.py --BamFile=input.bam --OutFile=output.bam

Opossum has been tested to work with BAM files that have been aligned with TopHat and STAR. Whether or not the aligner uses soft clips should be specified with the parameter *SoftClipsExist* (set *False* for TopHat, *True* for STAR).

The BAM file given as input to Opossum should be sorted according to read position. The reads should include the MD tag in the optional field. This can be done during alignment with STAR using `--outSAMattributes NH HI AS nM MD`. 


If a BAM file without MD tags already exists, they can be added as follows:

```bash
samtools calmd -b input.bam ref.fasta > opossum_input.bam
```

where input.bam is the file missing MD tags, ref.fasta is the reference fasta file, and opossum_input.bam is a BAM file Opossum can use.

Opossum does not currently handle reads with hard clips and these will be discarded.

Opossum expects that base qualities have been expressed in the standard Illumina encoding starting at 33. If another encoding has been used (such as the encoding used in earlier Illumina platforms, which starts at 64), the base qualities should be first converted to the standard encoding with e.g. GATK.

Opossum sorts the output file and creates a corresponding index file .bai.


You can see a list of all the possible input options by running the following command:

    python Opossum.py --help

You should adjust the input parameters according to which mapper has been used for aligning the reads in the input BAM file. The appropriate values for *MinFlankStart* and *MinFlankEnd* can be determined with the help of **compute_mismatch_rates.py** and **plot_mismatch_rates.py**, which are provided alongside the Opossum code.



## 3. Main command-line arguments to Opossum


*--BamFile* : (Required) BAM file name, or entire path if the file is not located in the same directory as where the code is run.

*--MinFlankStart* : Ignore base-changes closer than MinFlankStart bases to the start of the read. The corresponding base qualities will be set to 0. Default = 0

*--MinFlankEnd* : Ignore base-changes closer than MinFlankEnd bases to the end of the read. The corresponding base qualities will be set to 0. Default = 0

*--SoftClipsExist* : If set to False (default), no soft clips are expected in the reads, and MinFlank and MinFlankStart are applied to all reads. If set to True, then MinFlank and MinFlankStart are only applied to reads with soft clips. The value typically depends on the aligner used, e.g. for TopHat, it should be False, and for Star, it should be True.

*--MapCutoff* : Minimum mapping quality of a read pair. If either read in the pair has a lower mapping quality, both will be discarded. Default = 40

*--ProperlyPaired* : If set to True (default), only properly paired reads are considered.

*--KeepMismatches* : If set to False (default), overlapping paired-end reads having at least one base mismatch within the overlap region are discarded.

*--OutFile* : (Required) Output BAM file name, or entire path if the file will not be created in the same directory as where the code is run.



## 4. Suggested parameters when calling variants with Platypus


If, after running Opossum, variants are called using Platypus, we suggest keeping the Platypus settings at a minimum. This is because the reads have already been formatted in a way to take into account the various options included in Platypus.

    python Platypus.py --callVariants --bamFiles input.bam --refFile ref.fasta --filterDuplicates 0 --minMapQual 0 --minFlank 0 --maxReadLength 500 --minGoodQualBases 10 --minBaseQual 20 -o variants.vcf

It should be noted that *minGoodQualBases* also defines the minimum acceptable read length to be considered for variant calling. Using a value below the suggested 10 here may cause segmentation fault.

By default, Platypus flags variants that do not fulfill all of its filtering criteria (please see Platypus publication). These criteria have been designed to make the most out of DNA data. The same criteria can well be used with RNA-seq data if the user wants to maximize precision at the cost of sensitivity. However, if the user seeks a greater balance between precision and sensitivity, it would be advisable to include variants flagged as 'badReads', 'SC', and 'Q20' among the final variants.


## 5. Release History


### 0.2

Released on February 23, 2017. Updated dependency to Pysam v0.10.0. Opossum now supports unpaired data and filters out unmapped reads. Bug fixes.

### 0.1


Released in January 2016. First stable release of Opossum.



### Contact

Laura Oikkonen: firstname.surname (AT) well.ox.ac.uk


### License

Opossum is available under the GPL v3 license.


Please cite: Oikkonen L and Lise S. Making the most of RNA-seq: Pre-processing sequencing data with Opossum for reliable SNP variant detection. Wellcome Open Res 2017, 2:6 (doi: 10.12688/wellcomeopenres.10501.1) 

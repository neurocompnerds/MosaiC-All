# MosaiC-All

## Comprehensive analysis of somatic and parental gonosomal mosaicism using trio and singleton design.

### Step 1: Pre-requisites/ Softwares Installation

The following are softwares and resources that can be dowloaded/prepared based specified instructions.

#### 1.1 Downloadable sources<br>
MosaicHunter: https://github.com/zzhang526/MosaicHunter<br>
MosaicForecast: https://github.com/parklab/MosaicForecast<br>
GATK: https://github.com/broadinstitute/gatk/releases

#### 1.2 General Resources

|  Resources                    |     Example/Sources/Notes          | 
|-------------------------------|------------------------------------|  
|  Reference genome             |     e.g hs37d5.fa                  |
|  Common variant databases     |     dbSNP (e.g b37_dbsnp_138.b37.vcf)<br> gnomAD (e.g somatic-b37_af-only-gnomad.raw.sites.vcf)      |
|  Repeats                      |     e.g all_repeats.b37.bed;<br> can be found in MosaicHunter installation (i.e $MHDIR/MosaicHunter-master/resources) |
|  Exome errors databases       |     e.g WES_Agilent_71M.error_prone.b37.bed;<br> can be found in MosaicHunter installation (i.e $MHDIR/MosaicHunter-master/resources)                |
|  Panel of Normal (PON)        |     Should be prepared based on samples that are not part of the analysis.<br>As a suggestion for large cohort analysis, samples can be divided into two cohorts to create two Panel of Normals (PON_A and PON_B).<br>This PON can be prepared based on GATK option-CreateSomaticPanelOfNormals (i.e https://gatk.broadinstitute.org/hc/en-us/articles/4405451431963-CreateSomaticPanelOfNormals-BETA)          |

### Step 2: Config-files
There are two config-files, in which directories of softwares/resources (as prepared in Step 1) should be specified.

#### 2.1 Mosaic-All.config
#### 2.2 BWA-GATKHC.hs37d5_phoenix.cfg

### Step 3: Variant calling using four tools

#### 3.1 Summary: 
Variants are either detected in singletons and/or trios using three mosaic variant callers (MosaicHunter, MosaicForecast, Mutect2) and GATK-HC. 
A wrapper-script (MasterScript_MosaiC-All.sh) was developed for this purpose, and can be executable in University of Adelaide HPC environment using the following command.

#### 3.2 Command:

> MosaiC-All/MasterScript_MosaiC-All.sh -s $SampleID.list -o $Outputs -c MosaiC-All/Mosaic-All.config

#### 3.3 Requirements:-

1. $SampleID.list: A tab-separated-file as following format based on the Bam.files of each sample (e.g 001P.realigned.bam)

|  Directory of Bam files  | ProbandID | Gender   | MotherID | FatherID | 
|--------------------------|-----------|----------|----------|----------|
|   ./path                 |   001P    |   F      |  001M    |   001F   |

2. $Outputs: An Output directory to store all final outputs
   
3. $DIR/MosaiC-All/MosaiC-All.config: A config file prepared as mentioned in Step 1 and 2.
   
### Step 4: Analysis for somatic mosaicism (M3 pipeline)

#### 4.1 Processing all outputs 
Aims:
- to filter MFcalls manually and
- Merge all variants based on each tool.

> sbatch $SCRIPTDIR/CombineCalls.sh -s $sampleID -d $Outputs

Requirements:-
- sampleID (i.e 001P)
- Outputs (Output directory as specified in Step 3)

#### 4.2 Find overlaps using R script (M3pipeline.R)
Aims:
- To flag variants that were found in same sample, using 1-3 mosaic variant calling tools
- Followed by filtering out variants that were found only by one tool

### Step 5: Analysis for parental gonosomal mosaicism (pGoM)

> sbatch $SCRIPTDIR/pGoM.sh -v $VCFDIR -s $OUTDIR/FamilyID.txt -o $OUTDIR
- VCFDIR					(directory for family VCFs)
- FAMILYID   			(familyID.txt:list of familyIDs in a txt file; $FAMILYID*.vcf)
- OUTDIR					(output directory) 


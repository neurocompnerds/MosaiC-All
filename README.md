# MosaiC-All

## Comprehensive analysis of somatic and parental gonosomal mosaicism using trio or singleton design.

### Step 1: Pre-requisites/ Softwares Installation

The following packages and resources are required to run the pipelines. Packages should be installed according to developer instructions.

#### 1.1 Downloadable sources<br>
MosaicHunter: https://github.com/zzhang526/MosaicHunter<br>
MosaicForecast: https://github.com/parklab/MosaicForecast<br>
GATK: https://github.com/broadinstitute/gatk/releases

#### 1.2 General Resources

|  Resources                    |     Example/Sources/Notes          | 
|-------------------------------|------------------------------------|  
|  Reference genome             |     e.g hs37d5.fa                  |
|  Variant databases            |     dbSNP (e.g b37_dbsnp_138.b37.vcf)<br> gnomAD (e.g somatic-b37_af-only-gnomad.raw.sites.vcf)      |
|  Repeats                      |     e.g all_repeats.b37.bed;<br> can be found in MosaicHunter installation (i.e $MHDIR/MosaicHunter-master/resources) |
|  Exome errors databases       |     e.g WES_Agilent_71M.error_prone.b37.bed;<br> can be found in MosaicHunter installation (i.e $MHDIR/MosaicHunter-master/resources)                |
|  Panel of Normals (PON)        |     Should be prepared based on samples that are not part of the analysis.<br>As a suggestion for large cohort analysis, samples can be divided into two batches to create two PONs (PON_A and PON_B).<br> PONs are prepared based on GATK option-CreateSomaticPanelOfNormals (i.e https://gatk.broadinstitute.org/hc/en-us/articles/4405451431963-CreateSomaticPanelOfNormals-BETA)          |

### Step 2: Config-file
The config-file (Mosaic-All.config) is used to specify locations of required software and resources (as prepared in Step 1). 

MosaiC-ALL/config/Mosaic-All.config

### Step 3: Variant calling using four tools (for M3 and pGoM pipelines)

#### 3.1 Summary: 
Variant detection from either singleton or trio WES data is performed using three mosaic variant callers (MosaicHunter, MosaicForecast, Mutect2). Germline variants are called using GATK4 best practices workflow which should be run separately (https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels). 

The wrapper-script (MasterScript_MosaiC-All.sh) will run all tools for mosaic variant calling using the following command:

#### 3.2 Command:

> MosaiC-All/MasterScript_MosaiC-All.sh -s $SampleID.list -o $Outputs -c MosaiC-All/config/Mosaic-All.config

This script was designed for slurm workload manager in an HPC environment, but can run locally with some adjustments.

#### 3.3 Requirements:-

1. $SampleID.list: A tab-separated-file as following format based on the bam.files of each sample (e.g 001P.realigned.bam)

|  Directory of bam files  | ProbandID | Gender   | MotherID | FatherID | 
|--------------------------|-----------|----------|----------|----------|
|   ./path                 |   001P    |   F      |  001M    |   001F   |

2. $Outputs: An Output directory to store all final outputs
   
3. MosaiC-ALL/config/Mosaic-All.config: A config file prepared as described in Step 1 and 2.
   
### Step 4: Analysis for somatic mosaicism (M3 pipeline)

#### 4.1 Merging mosaic variant calls 
Aims:
- to filter MFcalls and
- Merge callsets from all tools.

> sbatch MosaiC-ALL/postprocessing/M3_CombineCalls.sh -s $sampleID -d $Outputs

Requirements:-
- sampleID (i.e 001P)
- Outputs (Output directory as specified in Step 3)

#### 4.2 Example R script for finding overlaps (M3pipeline.R)
Aims:
- To count how many tools called each variant
- Followed by filtering out variants that were called by only by one tool

The script MosaiC-ALL/postprocessing/M3pipeline.R is an example script for performing these filtering steps in R.


### Step 5: Analysis for parental gonosomal mosaicism (pGoM)

Aims: 
- To filter parental mosaic variant calls based on transmission to children

> sbatch /MosaiC-ALL/postprocessing/pGoM.sh -v $VCFDIR -s $OUTDIR/FamilyID.txt -o $OUTDIR
- VCFDIR					(directory for family VCFs)
- FAMILYID   			(familyID.txt:list of familyIDs in a txt file; $FAMILYID*.vcf)
- OUTDIR					(output directory) 


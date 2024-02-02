# MosaiC-All

## Comprehensive analysis of somatic and parental gonosomal mosaicism using trio or singleton design.

### Step 1: Pre-requisites/ Software Installation

The following packages and resources are required to run the pipelines. Packages should be installed according to developer instructions.

#### 1.1 Downloadable sources<br>
MosaicHunter: https://github.com/zzhang526/MosaicHunter<br>
MosaicForecast: https://github.com/parklab/MosaicForecast<br>
GATK: https://github.com/broadinstitute/gatk/releases

#### 1.2 General Resources

|  Resources                    |     Example/Sources/Notes          | 
|-------------------------------|------------------------------------|  
|  Reference genome             |     e.g hs37d5.fa                  |
|  Variant databases            |     dbSNP (e.g b37_dbsnp_138.b37.vcf)<br> gnomAD (e.g somatic-b37_af-only-gnomad.raw.sites.vcf). Obtain these from the [GATK resource bundle](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle)     |
|  Repeats                      |     e.g all_repeats.b37.bed;<br> can be found in the [MosaicHunter repository](https://github.com/zzhang526/MosaicHunter/tree/master/resources) |
|  Exome errors databases       |     e.g WES_Agilent_71M.error_prone.b37.bed;<br> can be found in the [MosaicHunter repository](https://github.com/zzhang526/MosaicHunter/tree/master/resources)                |
|  Panel of Normals (PON)        |     Should be prepared based on samples that are not part of the analysis.<br>As a suggestion for large cohort analysis, samples can be divided into two batches to create two PONs (PON_A and PON_B).<br> PONs are prepared based on GATK option-CreateSomaticPanelOfNormals (i.e https://gatk.broadinstitute.org/hc/en-us/articles/4405451431963-CreateSomaticPanelOfNormals-BETA)          |

### Step 2: Config-file
The config-file (Mosaic-All.config) is used to specify locations of required software and resources (as prepared in Step 1). 

MosaiC-ALL/config/Mosaic-All.config

### Step 3: Variant calling using three tools (for M3 and pGoM pipelines)

#### 3.1 Summary: 
Variant detection from either singleton or trio WES data is performed using three mosaic variant callers (MosaicHunter, MosaicForecast, Mutect2). Germline variants are called using GATK4 best practices workflow which should be run separately (https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels). 

The wrapper-script (MasterScript_MosaiC-All.sh) will run all tools for mosaic variant calling using the following command:

#### 3.2 Command:

`MosaiC-All/MasterScript_MosaiC-All.sh -s SampleID.list -o /path/to/output/directory -c MosaiC-All/config/Mosaic-All.config`

This script was designed for slurm workload manager in an HPC environment, but it can be adapted to run locally with some adjustments.

#### 3.3 Requirements:-

1. SampleID.list: A tab-separated-file as following format based on the bam.files of each sample (e.g 001P.realigned.bam)

|  Directory of bam files  | ProbandID | Gender   | MotherID | FatherID | 
|--------------------------|-----------|----------|----------|----------|
|   ./path                 |   001P    |   F      |  001M    |   001F   |
|   ./path                 |   003P    |   M      |  003M    |   003F   |

2. /path/to/output/directory: An output directory to store all final outputs
   
3. MosaiC-ALL/config/Mosaic-All.config: A config file prepared as described in Step 1 and 2.
   
### Step 4: Analysis for somatic mosaicism (M3 pipeline)

#### 4.1 Merging mosaic variant calls 
Aims:
- to filter MFcalls and
- Merge callsets from all tools.

Command

`sbatch MosaiC-ALL/postprocessing/M3_CombineCalls.sh -s sampleID -o /path/to/output/directory`

Requirements:

1. sampleID (i.e 001P)
2. /path/to/output/directory (Output directory as specified in Step 3)

#### 4.2 Example R script for finding overlaps (M3pipeline.R)
Aims:
- To count how many tools called each variant
- Followed by filtering out variants that were called by only by one tool

The script MosaiC-ALL/postprocessing/M3pipeline.R is an example script for performing these filtering steps in R.


### Step 5: Analysis for parental gonosomal mosaicism (pGoM)

Aims: 
- To filter parental mosaic variant calls based on transmission to children

5.1 Prefilter
- Identify inherited variants that expected to be mosaic based on AAF and GT, using GATKHC outputs

Command:

`sbatch /MosaiC-ALL/postprocessing/pGoM.sh -v /path/to/directory_of_vcf -s FamilyID.txt -o /path/to/output/directory`

Requirements:

1. input directory (where can we find the family.vcf).

2. sampleID list (one header row and then tab-delimited columns \$BAMdir,\$ProbandID,\$Gender,\$Mother,\$Father).
   
|  Directory of Bam files  | ProbandID | Gender   | MotherID | FatherID | FamilyVCF | 
|--------------------------|-----------|----------|----------|----------|-----------|
|   ./path                 |   001P    |   F      |  001M    |   001F   | Trio001.vcf |
|   ./path                 |   004P    |   F      |  004M    |   004F   | 004.family.vcf |

3. /path/to/output/directory	(A location for the output files).

5.2 Postfilter
- Mosaic variants are identified among prefiltered pGoM variants using one or more mosaic variant calling tools
- Example script: MosaiC-ALL/postprocessing/pGoMpipeline.R

Requirements:
1. Three Output files from MosaiC-ALL/postprocessing/M3_CombineCalls.sh
2. pGoM.sh output file
3. Amend the working directory and output_file prefix in the R.script

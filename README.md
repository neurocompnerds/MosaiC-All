# MosaiC-All

## Comprehensive analysis of somatic and parental gonosomal mosaicism using trio or singleton design.

### Step 1: Pre-requisites/ Software Installation

Packages and resources required to run the pipelines are listed in section 1.1 and 1.2. 

To enable testing of software installation and set up of the MosaiC-All pipelines, we have provided:

1. Toy data adapted from https://github.com/parklab/MosaicForecast/tree/master/demo
2. Commands and resources used to run the M3 pipeline on this data.
3. Output from all tools and resulting filtered variant list after post-processing.

These data and resources, along with further instructions, are provided in MosaiC-All/TestRun/
  

#### 1.1 Software Installation

Packages should be installed according to developer instructions. You can find these, along with documentation for each tool at:

1. MosaicHunter: https://github.com/zzhang526/MosaicHunter<br>
2. MosaicForecast: https://github.com/parklab/MosaicForecast<br>
3. GATK: https://github.com/broadinstitute/gatk/releases

We recommend installation via the docker images for each tool, which are available in the following locations:<br>

MosaicHunter: https://hub.docker.com/r/rborgesm/mosaichunter<br>

MosaicForecast: https://hub.docker.com/r/yanmei/mosaicforecast<br>

GATK: https://hub.docker.com/r/broadinstitute/gatk/tags<br>

#### 1.2 General Resources

The following resources are required and should be downloaded or generated:

|  Resources                    |     Example/Sources/Notes          | 
|-------------------------------|------------------------------------|  
|  Reference genome             |     Fasta file for the genome build to which your data was mapped e.g hs37d5.fa                  |
|  Variant databases            |     Variant population frequencies: dbSNP (e.g b37_dbsnp_138.b37.vcf), gnomAD (e.g somatic-b37_af-only-gnomad.raw.sites.vcf). Obtain these from the [GATK resource bundle](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle)     |
|  Repeat regions                    |     Repeat regions for filtering <br> e.g all_repeats.b37.bed; can be found in the [MosaicHunter repository](https://github.com/zzhang526/MosaicHunter/tree/master/resources) |
|  Exome error prone regions       |     Regions known to be error prone for filtering <br> e.g WES_Agilent_71M.error_prone.b37.bed; can be found in the [MosaicHunter repository](https://github.com/zzhang526/MosaicHunter/tree/master/resources)                |
|  Panel of Normals (PON)       |     PON should be prepared using samples that are not part of the analysis.As a suggestion for large cohort analysis, samples can be divided into two batches to create two PONs (PON_A and PON_B).<br> PONs are prepared based on GATK option-CreateSomaticPanelOfNormals (i.e https://gatk.broadinstitute.org/hc/en-us/articles/4405451431963-CreateSomaticPanelOfNormals-BETA)          |

### Step 2: Config-file
The config-file (MosaiC-ALL/config/Mosaic-All.config) is used to specify locations of required software and resources (as prepared in Step 1). 
Specify the locations of resources, genomebuild etc using the template config file as per included instructions.

### Step 3: Mosaic variant calling using three tools (for M3 and pGoM pipelines)

#### 3.1 Summary: 
Variant detection from either singleton or trio WES data is performed using three mosaic variant callers (MosaicHunter, MosaicForecast, Mutect2). Germline variants are called using GATK4 best practices workflow which should be run separately (https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels). 

The wrapper-script (MasterScript_MosaiC-All.sh) will run all tools for mosaic variant calling using the following command:

#### 3.2 Command:

`MosaiC-All/MASTERSCRIPT_MosaiC-All.sh -s SampleID.list -o /path/to/output/directory -c MosaiC-All/config/Mosaic-All.config`

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
- To filter MF calls and
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

#### 5.1 Prefilter
- Identify inherited variants that are predicted to be mosaic in a parent based on AAF and GT, using GATKHC outputs

Command:

`sbatch /MosaiC-ALL/postprocessing/pGoM.sh -v /path/to/directory_of_vcf -s FamilyID.txt -o /path/to/output/directory`

Requirements:

1. input directory (where to find the family.vcf).

2. sampleID list (one header row and then tab-delimited columns \$BAMdir,\$ProbandID,\$Gender,\$Mother,\$Father).
   
|  Directory of Bam files  | ProbandID | Gender   | MotherID | FatherID | FamilyVCF | 
|--------------------------|-----------|----------|----------|----------|-----------|
|   ./path                 |   001P    |   F      |  001M    |   001F   | Trio001.vcf |
|   ./path                 |   004P    |   F      |  004M    |   004F   | 004.family.vcf |

3. /path/to/output/directory	(A location for the output files).

#### 5.2 Postfilter
- Mosaic variants are identified from prefiltered pGoM variants using one or more mosaic variant calling tools
- Example script: MosaiC-ALL/postprocessing/pGoMpipeline.R

Requirements:
1. Three Output files from MosaiC-ALL/postprocessing/M3_CombineCalls.sh
2. pGoM.sh output file
3. Amend the working directory and output_file prefix in the R.script


# MosaiC-All

## Comprehensive analysis of somatic and parental gonosomal mosaicism using trio and singleton design.

The analysis is divided into three sections (Part A,B and C).

### Part A: Variant calling using four existing tools

Summary: 
Variants are either detected in singletons and/or trios using three mosaic variant callers (MosaicHunter, MosaicForecast, Mutect2) and GATK-HC. 
The results of the entire batch (as listed in SampleID) are stored in three separate files (MH.calls.txt,MF.calls.txt,Mutect2.,calls.txt), and the number of variants of each sample in Counts.txt
A wrapper-script (MasterScript_Mosaic_All.sh) was developed for this purpose, and can be executable in University of Adelaide HPC environment using the following command.

Wrapper-script Progress (will be deleted after finalising)
- Mosaic variant callers are tested on HPC, but not GATK-HC. 
- still need to include Variant Counting and coverage analysis

Command:

> MosaiC-All/MasterScript_Mosaic_All.sh -s SampleID -o Outputs -c MosaiC-All/Mosaic-All.config

Requirements:-
1. SampleID file (tab-separated-file) should be like this:- 

|  Directory of Bam files  | ProbandID | Gender   | MotherID | FatherID | 
|--------------------------|-----------|----------|----------|----------|
|   ./path                 |   001P    |   F      |  001M    |   001F   |

2. Necessary details/pathways are specified in the config
   - Reference genome
   - Installation of MosaicHunter, MosaicForecast, GATK (Mutect2 and HC)
   - Panel of Normals (PON)
     -> for Mutect2;
     -> in our study, we used two PONs developed using parental data of respective cohorts, to ensure samples that are accessed by Mutect2 for mosaic variants are not in PON.

3. Output directory: To store all final outputs
   
### Part B: Analysis for somatic mosaicism (M3 pipeline)
Summary: Using Outputs in Part A, R is used to find the overlap between variant callers.

### Part C: Analysis for parental gonosomal mosaicims (pGoM)
SCript npt included yet.



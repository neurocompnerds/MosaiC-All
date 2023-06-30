# MosaiC-All

## Comprehensive analysis of somatic and parental gonosomal mosaicism using trio and singleton design.

### Part A: Variant calling using four existing tools

Summary: 
Variants are either detected in singletons and/or trios using three mosaic variant callers (MosaicHunter, MosaicForecast, Mutect2) and GATK-HC. 
The results of the entire batch (as listed in SampleID) are stored in three separate files (MH.calls.txt,MF.calls.txt,Mutect2.,calls.txt), and the number of variants of each sample in Counts.txt

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

2. Specify the necessary details/pathways in the config
   - Reference genome
   - Installation of MosaicHunter, MosaicForecast, GATK (Mutect2 and HC)
   - Panel of Of Normals (for Mutect2); see (need to include link to Mutect2 page)

3. Output directory: To store all final outputs
   
### Part B: Analysis for somatic mosaicism (M3 pipeline)
Summary: Using Outputs in Part A, R is used to find the overlap.

### Part C: Analysis for parental gonosomal mosaicims (pGoM)
SCript npt included yet.



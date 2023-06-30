# MosaiC-All

## Comprehensive analysis of somatic and parental gonosomal mosaicism using trio and singleton design.

### Part A: Variant calling using four existing tools

Summary: 
To detect variants either in singletons or trios using three mosaic variant callers (MosaicHunter, MosaicForecast, Mutect2) and GATK-HC

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
### Part C: Analysis for parental gonosomal mosaicims (pGoM)

#Part A

Output:

#Part B
1. First item
2. Second item
3. Third item
    1. Indented item
    2. Indented item
4. Fourth item

#Part C

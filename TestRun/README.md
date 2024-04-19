## Demo using a test data (./MosaiC-All/TestRun/)

To test software installation and variant calling process, we have provided test data was adapted from MosaicForecast (https://github.com/parklab/MosaicForecast). 

### 1.Config and Resources
The resources used to generated the outputs for testdata can be found in the ./MosaiC-All/TestRun/Resources. 

Please specify the following resources in the **config** file:-

- RESOURCES="./MosaiC-All/TestRun/Resources"

- DBSNP=${RESOURCES}/b37_dbsnp_138.b37.vcf

- REPEATS=${RESOURCES}/all_repeats.b37.bed
  _#This can be found in subfolder of MosaicHunter (./MosaicHunter-master/resources)_

- COMMONERROR=${RESOURCES}/WES_Agilent_71M.error_prone.b37.bed
  _#This can be found in subfolder of MosaicHunter (./MosaicHunter-master/resources)_

- PON_A=${RESOURCES}/somatic-b37_Mutect2-exome-panel.vcf
  _#This can be downloaded from https://storage.googleapis.com/gatk-best-practices/somatic-b37_

- PON_B=${RESOURCES}/somatic-b37_Mutect2-exome-panel.vcf
  _#For testdata purpose, we just define the same Panel of Normals as used for $PON_A_

- GERMLINE_RESOURCES=${RESOURCES}/somatic-b37_af-only-gnomad.raw.sites.vcf
  _#This can be downloaded from https://storage.googleapis.com/gatk-best-practices/somatic-b37/af-only-gnomad.raw.sites.vcf_

### 2.Command

Execute the following command, which will perform variant calling and store the outputs in the $GIT/TestRun/Output

GIT=./MosaiC-All
CONFIG=$GIT/config/Mosaic-All.config
$GIT/MASTERSCRIPT_MosaiC-All.sh -s $GIT/TestRun/SampleID_Test -o $GIT/TestRun/Output -c $CONFIG

### 3.Outputs

The expected outputs are stored in $GIT/TestRun/Output_Expected, which can be used to compare to your outputs. We have also provided the list of variants after performing the post-processing step (FinalList_M3).


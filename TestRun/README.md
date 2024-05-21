## Demo using a test data (./MosaiC-All/TestRun/)

To test the variant calling process using manual installation of softwares, we have provided a test data which was extracted from MosaicForecast publication and following instructions

### 1.Pre-requisite

Packages should be installed according to developer instructions. You can find these, along with documentation for each tool at:

1. MosaicHunter: https://github.com/zzhang526/MosaicHunter<br>
2. MosaicForecast: https://github.com/parklab/MosaicForecast<br>
3. GATK: https://github.com/broadinstitute/gatk/releases

### 2.Config

Specify the directories in which the MosaiC-All git and variant calling tools softwares were downloaded as well as where log_files should be stored in ./MosaiC-All/config/Mosaic-All.TestRun.config.

### 3.Command

Execute the following command, which will download required resources, perform variant calling and store the outputs in the $GIT/TestRun/OUTPUT

git=/path/to/MosaiC-All

```$git/TestRun/TestRun_manual.sh -d $git -s $git/TestRun/SampleID_Test -c $git/config/Mosaic-All.TestRun.config -o $git/TestRun/OUTPUT```

### 4.Outputs

The expected outputs are stored in $gitTestRun/Output_Expected, which can be used to compare with the given outputs. Also, we provide the M3 variants upon post-processing step (FinalList_M3).


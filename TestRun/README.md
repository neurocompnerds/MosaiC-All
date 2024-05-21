## Demo using test data (./MosaiC-All/TestRun/)

To test the variant calling process following manual installation of software, we have provided test data extracted from the MosaicForecast publication.

### 1.Pre-requisites

Packages should be installed according to developer instructions. You can find these, along with documentation for each tool at:

1. MosaicHunter: https://github.com/zzhang526/MosaicHunter<br>
2. MosaicForecast: https://github.com/parklab/MosaicForecast<br>
3. GATK: https://github.com/broadinstitute/gatk/releases

### 2.Config

Specify the directories in which the MosaiC-All git and variant calling tools were downloaded, as well as where log_files should be stored in ./MosaiC-All/config/Mosaic-All.TestRun.config.

### 3.Command

Execute the following command, which will download required resources, perform variant calling and store the outputs in $git/TestRun/OUTPUT where git=/path/to/MosaiC-All

```$git/TestRun/TestRun_manual.sh -d $git -s $git/TestRun/SampleID_Test -c $git/config/Mosaic-All.TestRun.config -o $git/TestRun/OUTPUT```

### 4.Outputs

The expected outputs can be found in $git/TestRun/Output_Expected, which can be used to compare with your results. We have also provided the list of M3 variants after the post-processing step (FinalList_M3).


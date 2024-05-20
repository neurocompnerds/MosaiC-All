## Demo using docker (./MosaiC-All/TestRun_docker/)

To test the docker-installations and variant calling process, we have provided a test data which was extracted from MosaicForecast publication. 

### 1.Pre-requisites: Make sure to have these softwares installed prior to running the command-line
- samtools v1.20 (using htslib 1.20)
- bcftools v1.20 (using htslib 1.20)
- java (v16.0.1)
- wget
- gunzip
  
### 2.Command

Execute the following command, which will perform installation of dockers and resources required by each variant calling tool.

```./MosaiC-All/TestRun_docker/TestRun_docker.sh -d /path/to/MosaiC-All -s 1465_1024-pfc-bulk.rehead.bam -i 1465_1024-pfc-bulk -g M ```

### 3.Outputs

The expected outputs are stored in $GIT/TestRun/Output_Expected, which can be used to compare with the given outputs. Also, we provide the M3 variants upon post-processing step (FinalList_M3).


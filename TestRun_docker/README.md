## Demo using docker (./MosaiC-All/TestRun_docker/)

To test the docker-installations and variant calling process, we have provided test data which was extracted from the MosaicForecast publication. 

### 1.Command

Execute the following command, which will perform installation of dockers and download resources required by each variant calling tool.

```./MosaiC-All/TestRun_docker/TestRun_docker.sh -d /path/to/MosaiC-All -s 1465_1024-pfc-bulk.rehead.bam -i 1465_1024-pfc-bulk -g M ```

### 2.Outputs

The expected outputs can be found in $git/TestRun/Output_Expected, which can be used to compare with your results. We have also provided the list of M3 variants after the post-processing step (FinalList_M3).

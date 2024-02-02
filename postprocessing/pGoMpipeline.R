#README
#SCRIPT: pGoM pipeline.R 
#AIM: To filter variants that are mosaic in parents and germline in children
#REQUIREMENT: 
#(i) Three sets of variant lists (Mutect2, MF and MH) of parents that were identified and categorised based on each tool
#(ii) The output of pGoM.sh (single.file: $ProbandID.pGoM.vcf)
#AMEND/EDIT: STEP2 (directory) and STEP 3 (list of variants)

# 1. LOAD LIBRARY
library(tidyverse)

# 2. FILL UP THE VARIABLES HERE
setwd("YOUR_DIRECTORY_PATH_HERE")
Mut2_List<-"List_Mutect2variants.txt" 
MH_List<-"List_MHvariants.txt"
MF_List<-"List_MFvariants.txt"
prefilter_pGoM_List<-"List_prefilter_pGoMvariants.txt" ## Output of pGoM.sh
pGoM_finalOutput<-"pGoM_finalOutput.txt"

# 3. TIDY (Preparation before pGoM step)
# load three list of variants detected by each variant caller 
# and standardised some of the headers (so, check first)
# important columns required: "CHROM","POS","REF","ALT","SampleID","VC"

## 3.1 Mosaicvariant callers
Mutect2<-read.delim(file=Mut2_List)%>%
            rename("CHROM"="X.CHROM",
                   "SampleID"="Header_SAMPLEID", #amend this based on how the first column header was defined
                   "VC"="Mut",
                   "Mutect2_FORMAT"="SAMPLEID")%>%
            filter(CHROM!="#CHROM")

MH<-read.delim(file=MH_List, header=FALSE)%>%
            rename("SampleID"="V1",
                   "VC"="V2",
                   "CHROM"="V3",
                   "POS"="V4",
                   "REF"="V5",
                   "ALT"="V6")

MF<-read.delim(file="MF_List",header=FALSE)%>%
            separate(V3, c("file", "CHROM", "POS", "REF", "ALT"), sep = "~")

# 3.2 Prefiltered pGoM from GATKHC
GATKHC<-read.delim(file=prefilter_pGoM_List, header=TRUE)
# Remove "X" prefix and "." prefix from column names
names(GATKHC)[-1] <- gsub("^X", "", names(GATKHC)[-1])
names(GATKHC)[-1] <- gsub("^\\.", "", names(GATKHC)[-1])
#function to tidy the data to ensure there the FORMAT columns are separated based on SampleID specified in first column
separate_columns <- function(data) {
  unique_sample_ids <- unique(data$SampleID) 
  
  result <- data %>%
    mutate(across(starts_with(unique_sample_ids), ~ifelse(SampleID == sub("\\s.*", "", SampleID), ., NA_character_))) %>%
    separate(starts_with(unique_sample_ids), into = c("GT", "AD", "DP", "GQ", "PL"), sep = ":", remove=FALSE) %>%
    mutate(across(starts_with(unique_sample_ids), ~ifelse(is.na(.), "", .)))%>%
    separate(AD, into = c("AD1","AD2"))%>%
    mutate(AD1 = as.numeric(AD1), AD2 = as.numeric(AD2), DP = as.numeric(DP))%>%
    mutate(AAF = AD2 / DP)
  
  return(result)
}

#Execute the function of tidy-up the data
GATKHC_Tidy<- separate_columns(GATKHC)

# 4. MINOR FORMATTING
# amend the format of "POS" columns of Mutect2 and MF files
Mutect2$POS<-as.numeric(Mutect2$POS)
MF$POS<-as.numeric(MF$POS)
GATKHC$POS<-as.numeric(GATKHC$POS)

# 5.FLAG THE OVERLAP DETECTION (each tool and prefilter_pGoM) and TIDY UP
##5.1 Single Variant caller tool
Mut_pGoM <- inner_join(Mutect2, GATKHC, by = c("CHROM", "POS", "SampleID", "REF", "ALT"))%>%mutate(VC_MUT="MUT")
  MH_pGoM <- inner_join(MH, GATKHC, by = c("CHROM", "POS", "SampleID", "REF", "ALT"))%>%mutate(VC_MH="MH")
  MF_pGoM <- inner_join(MF, GATKHC, by = c("CHROM", "POS", "SampleID", "REF", "ALT"))%>%mutate(VC_MF="MF")

##5.2 More than 2 Variant callers
#merge three VC outputs (from 5.1), which will give three columns of VC
MF_MH<-full_join(MH_pGoM,MF_pGoM,by=c("CHROM", "POS", "SampleID", "REF", "ALT"))
  MF_MH_MUT<-full_join(MF_MH,Mut_pGoM,by=c("CHROM", "POS", "SampleID", "REF", "ALT"))
#flag VC that were not identified by respective tool
MF_MH_MUT$VC_MUT <- ifelse(is.na(MF_MH_MUT$VC_MUT), "NO", MF_MH_MUT$VC_MUT)
  MF_MH_MUT$VC_MH <- ifelse(is.na(MF_MH_MUT$VC_MH), "NO", MF_MH_MUT$VC_MH)
  MF_MH_MUT$VC_MF <- ifelse(is.na(MF_MH_MUT$VC_MF), "NO", MF_MH_MUT$VC_MF)
#Flag variants based on which tool or tools identified
ALL<-filter(MF_MH_MUT,VC_MF=="MF" & VC_MH=="MH" & VC_MUT=="MUT")%>%mutate(Group="ALL")
  MUT_MH<-filter(MF_MH_MUT,VC_MF=="NO" & VC_MH=="MH" & VC_MUT=="MUT")%>%mutate(Group="MUT_MH")
  MUT_MF<-filter(MF_MH_MUT,VC_MF=="MF" & VC_MH=="NO" & VC_MUT=="MUT")%>%mutate(Group="MUT_MF")
  MH_MF<-filter(MF_MH_MUT,VC_MF=="MF" & VC_MH=="MH" & VC_MUT=="NO")%>%mutate(Group="MH_MF")
  MH_only<-filter(MF_MH_MUT,VC_MF=="NO" & VC_MH=="MH" & VC_MUT=="NO")%>%mutate(Group="MH_ONLY")
  MUT_only<-filter(MF_MH_MUT,VC_MF=="NO" & VC_MH=="NO" & VC_MUT=="MUT")%>%mutate(Group="MUT_ONLY")
  MF_only<-filter(MF_MH_MUT,VC_MF=="MF" & VC_MH=="NO" & VC_MUT=="NO")%>%mutate(Group="MF_ONLY")
  NONE<-filter(MF_MH_MUT,VC_MF=="NO" & VC_MH=="NO" & VC_MUT=="NO")%>%mutate(Group="NONE")

#Compile all into one file   
pGoM<-rbind(ALL,MUT_MH,MUT_MF,MH_MF,MH_only,MUT_only,MF_only,NONE)%>%unique()

# 6. SAVE FILE
write.csv(pGoM, file = paste0(pGoM_finalOutput, ".csv"))

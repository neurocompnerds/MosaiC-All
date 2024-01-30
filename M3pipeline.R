#README
#SCRIPT: M3 pipeline.R 
#AIM: To flag variants that are identified by more than one variant calling tool
#REQUIREMENT: Three sets of variant lists (Mutect2, MF and MH) that were identified and categorised based on each tool
#AMEND/EDIT: STEP2 (directory) and STEP 3 (list of variants)

# 1. LOAD LIBRARY
library(tidyverse)

# 2. SET DIRECTORY
setwd("YOUR_DIRECTORY_PATH_HERE")

# 3. LOAD RAW FILE & TIDY
# load three list of variants detected by each variant caller 
# and standardised some of the headers (so, check first)
# important columns required: "CHROM","POS","REF","ALT","SampleID","VC"
Mutect2<-read.delim(file="Variants.Mutect2.txt")%>%
            rename("CHROM"="X.CHROM",
                   "SampleID"="Header_SAMPLEID",
                   "VC"="Mut",
                   "Proband"="SAMPLEID")%>%
            filter(CHROM!="#CHROM")
MH<-read.delim(file="Variants.MH.txt", header=FALSE)%>%
            rename("SampleID"="V1",
                   "VC"="V2",
                   "CHROM"="V3",
                   "POS"="V4",
                   "REF"="V5",
                   "ALT"="V6")
MF<-read.delim(file="Variants.MF.txt",header=FALSE)%>%
            separate(V3, c("file", "CHROM", "POS", "REF", "ALT"), sep = "~")

# 4. MINOR FORMATTING
# amend the format of "POS" columns of Mutect2 and MF files
Mutect2$POS<-as.numeric(Mutect2$POS)
MF$POS<-as.numeric(MF$POS)

# 5.FLAG THE OVERLAP DETECTION and TIDY UP

## 5.1. Mutect2 and MH
### Find the overlap between two tools and calculate the AAF, and extract only certain columns
Mut_MH <- inner_join(Mutect2, MH, by = c("CHROM", "POS", "SampleID", "REF", "ALT"))%>%
  separate(Proband, into = c("GT", "AD", "AF", "DP", "F1R2", "F2R1", "FAD", "SB"), sep = ":", remove=FALSE)%>%
  separate(AD, into = c("AD1", "AD2"), sep = ",")%>%
  mutate(AD1 = as.numeric(AD1), AD2 = as.numeric(AD2), DP = as.numeric(DP))%>%
  mutate(AAF = (AD2 / DP) * 100)%>%mutate(Category="MH_Mut")%>%
  select("CHROM","POS","ID", "REF","ALT","QUAL","FILTER","INFO","FORMAT","Proband","AAF","SampleID","Category")

## 5.2. Mutect2 and MF
### Find the overlap between two tools and calculate the AAF and extract only certain columns
Mut_MF <- inner_join(Mutect2, MF, by = c("CHROM", "POS", "SampleID", "REF", "ALT"))%>%
  separate(Proband, into = c("GT", "AD", "AF", "DP", "F1R2", "F2R1", "FAD", "SB"), sep = ":", remove=FALSE)%>%
  separate(AD, into = c("AD1", "AD2"), sep = ",")%>%
  mutate(AD1 = as.numeric(AD1), AD2 = as.numeric(AD2), DP = as.numeric(DP))%>%
  mutate(AAF = (AD2 / DP) * 100)%>%mutate(Category="MF_Mut")%>%
  select("CHROM","POS","ID", "REF","ALT","QUAL","FILTER","INFO","FORMAT","Proband","AAF","SampleID","Category")

## 5.3. MF and MH, and ALL
### Find the overlap between MF and MH, extract certain columns, followed by flagging for variants that are also detected by Mutect2
MH_MF_raw <- inner_join(MF, MH, by = c("CHROM", "POS", "SampleID", "REF", "ALT"))%>%
            mutate(Category="MF_MH")%>%
            select("CHROM","POS", "REF","ALT","AF","SampleID","Category")%>%
            left_join(Mut_MF,by = c("CHROM", "POS", "REF", "ALT", "SampleID"))

### filter MH_MF list for variants that are also flagged as MF_Mut
All<-filter(MH_MF_raw,Category.y=="MF_Mut")

### Remove duplicates by removing only flagging variants that are detected by ALL.
MH_MF<-anti_join(MH_MF_raw,All,by = c("CHROM", "POS", "REF", "ALT", "SampleID"))

# 6. SAVE FILE
write_delim(Mut_MH,file="Mut_MH.txt")
write_delim(Mut_MF,file="Mut_MF.txt")
write_delim(All,file="All.txt")
write_delim(MH_MF,file="MH_MF.txt")

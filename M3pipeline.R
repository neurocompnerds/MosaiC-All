#M3 pipeline: 
#Overlap between two or more mosaic variant caller, flagging each group

#set the directory
setwd("D:/Genomic Autopsy")

library(tidyverse)

#loadfiles which has list of variants detected by each variant caller and tidy up accordingly.
Mutect2<-read.delim(file="Variants.Mutect2.Batch1and2.txt")%>%
            rename("CHROM"="X.CHROM","SampleID"="PED004D","VC"="Mut","Proband"="SCO_PED004D_GA0014_D1")%>%
            filter(CHROM!="#CHROM")
MH<-read.delim(file="Variants.MH.Batch1and2.txt", header=FALSE)%>%
            rename("SampleID"="V1","VC"="V2","CHROM"="V3","POS"="V4","REF"="V5","ALT"="V6")
MF<-read.delim(file="Variants.MF.Batch1and2.txt",header=FALSE)%>%
            separate(V3, c("file", "CHROM", "POS", "REF", "ALT"), sep = "~")

#amend the types
Mutect2$POS<-as.numeric(Mutect2$POS)
MF$POS<-as.numeric(MF$POS)

#find the overlap
Mut_MH <- inner_join(Mutect2, MH, by = c("CHROM", "POS", "SampleID", "REF", "ALT"))%>%
  separate(Proband, into = c("GT", "AD", "AF", "DP", "F1R2", "F2R1", "FAD", "SB"), sep = ":", remove=FALSE)%>%
  separate(AD, into = c("AD1", "AD2"), sep = ",")%>%
  mutate(AD1 = as.numeric(AD1), AD2 = as.numeric(AD2), DP = as.numeric(DP))%>%
  mutate(AAF = (AD2 / DP) * 100)%>%mutate(Category="MH_Mut")%>%
  select("CHROM","POS","ID", "REF","ALT","QUAL","FILTER","INFO","FORMAT","Proband","AAF","SampleID","Category")

Mut_MF <- inner_join(Mutect2, MF, by = c("CHROM", "POS", "SampleID", "REF", "ALT"))%>%
  separate(Proband, into = c("GT", "AD", "AF", "DP", "F1R2", "F2R1", "FAD", "SB"), sep = ":", remove=FALSE)%>%
  separate(AD, into = c("AD1", "AD2"), sep = ",")%>%
  mutate(AD1 = as.numeric(AD1), AD2 = as.numeric(AD2), DP = as.numeric(DP))%>%
  mutate(AAF = (AD2 / DP) * 100)%>%mutate(Category="MF_Mut")%>%
  select("CHROM","POS","ID", "REF","ALT","QUAL","FILTER","INFO","FORMAT","Proband","AAF","SampleID","Category")

MH_MF_raw <- inner_join(MF, MH, by = c("CHROM", "POS", "SampleID", "REF", "ALT"))%>%
            mutate(Category="MF_MH")%>%
            select("CHROM","POS", "REF","ALT","AF","SampleID","Category")%>%
            left_join(Mut_MF,by = c("CHROM", "POS", "REF", "ALT", "SampleID"))

All<-filter(MH_MF_raw,Category.y=="MF_Mut")

MH_MF<-anti_join(MH_MF_raw,All,by = c("CHROM", "POS", "REF", "ALT", "SampleID"))

#save each file for annotations
write_delim(Mut_MH,file="GA_Mut_MH_batch1and2_2.txt")
write_delim(Mut_MF,file="GA_Mut_MF_batch1and2.txt")
write_delim(All,file="GA_All_batch1and2.txt")
write_delim(MH_MF,file="GA_MH_MF_batch1and2.txt")

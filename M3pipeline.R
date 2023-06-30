#M3 pipeline: Overlap between two or more mosaic variant caller

setwd("D:/Genomic Autopsy")

library(tidyverse)

#loadfiles
GA_Mutect2<-read.delim(file="Variants.Mutect2.Batch1and2.txt")%>%
            rename("CHROM"="X.CHROM","SampleID"="PED004D","VC"="Mut","Proband"="SCO_PED004D_GA0014_D1")%>%filter(CHROM!="#CHROM")
GA_MH<-read.delim(file="Variants.MH.Batch1and2.txt", header=FALSE)%>%
          rename("SampleID"="V1","VC"="V2","CHROM"="V3","POS"="V4","REF"="V5","ALT"="V6")
GA_MF<-read.delim(file="Variants.MF.Batch1and2.txt")%>%
        separate(filename, c("file", "CHROM", "POS", "REF", "ALT"), sep = "~")

#amend the types
GA_Mutect2$POS<-as.numeric(GA_Mutect2$POS)
GA_MF$POS<-as.numeric(GA_MF$POS)


#find the overlap
GA_Mut_MH <- inner_join(GA_Mutect2, GA_MH, by = c("CHROM", "POS", "SampleID", "REF", "ALT"))%>%
  separate(Proband, into = c("GT", "AD", "AF", "DP", "F1R2", "F2R1", "FAD", "SB"), sep = ":", remove=FALSE)%>%
  separate(AD, into = c("AD1", "AD2"), sep = ",")%>%
  mutate(AD1 = as.numeric(AD1), AD2 = as.numeric(AD2), DP = as.numeric(DP))%>%
  mutate(AAF = (AD2 / DP) * 100)%>%mutate(Category="MH_Mut")%>%
  select("CHROM","POS","ID", "REF","ALT","QUAL","FILTER","INFO","FORMAT","Proband","AAF","SampleID","Category")
write_delim(GA_Mut_MH,file="GA_Mut_MH_batch1and2_2.txt")

GA_Mut_MF <- inner_join(GA_Mutect2, GA_MF, by = c("CHROM", "POS", "SampleID", "REF", "ALT"))%>%
  separate(Proband, into = c("GT", "AD", "AF", "DP", "F1R2", "F2R1", "FAD", "SB"), sep = ":", remove=FALSE)%>%
  separate(AD, into = c("AD1", "AD2"), sep = ",")%>%
  mutate(AD1 = as.numeric(AD1), AD2 = as.numeric(AD2), DP = as.numeric(DP))%>%
  mutate(AAF = (AD2 / DP) * 100)%>%mutate(Category="MF_Mut")%>%
  select("CHROM","POS","ID", "REF","ALT","QUAL","FILTER","INFO","FORMAT","Proband","AAF","SampleID","Category")
write_delim(GA_Mut_MF,file="GA_Mut_MF_batch1and2.txt")

GA_MH_MF_raw <- inner_join(GA_MF, GA_MH, by = c("CHROM", "POS", "SampleID", "REF", "ALT"))%>%
            mutate(Category="MF_MH")%>%
            select("CHROM","POS", "REF","ALT","AF","SampleID","Category")%>%
            left_join(GA_Mut_MF,by = c("CHROM", "POS", "REF", "ALT", "SampleID"))

GA_All<-filter(GA_MH_MF_raw,Category.y=="MF_Mut")
write_delim(GA_All,file="GA_All_batch1and2.txt")

GA_MH_MF<-anti_join(GA_MH_MF_raw,GA_All,by = c("CHROM", "POS", "REF", "ALT", "SampleID"))
write_delim(GA_MH_MF,file="GA_MH_MF_batch1and2.txt")

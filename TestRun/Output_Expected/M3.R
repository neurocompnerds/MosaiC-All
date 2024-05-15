setwd("D:/WRITE UPS/manuscript for checking/ToyData")

#load tables
Mutect2_Test<-read.delim(file="1465_1024-pfc-bulk.mutect2.singlemode.PASS.aaf.vcf")%>%
  mutate(VC="Mutect2")%>%
  rename("CHROM"="X.CHROM")

  filter(CHROM!="#CHROM")

MH_Test<-read.delim(file="1465_1024-pfc-bulk.final.passed.tsv", header=FALSE)%>%
  rename("CHROM"="V1",
         "POS"="V2",
         "REF"="V3",
         "ALT"="V9")%>%mutate(VC="MH")

MF_Test<-read.delim(file="1465_1024-pfc-bulk.mosaicforecast.genotype.predictions.refined.bed")%>%
  separate(id, c("file", "CHROM", "POS", "REF", "ALT"), sep = "~")%>%
  mutate(VC="MF")

Mutect2_Test$POS<-as.numeric(Mutect2_Test$POS)
MF_Test$POS<-as.numeric(MF_Test$POS)
MF_Test$CHROM<-as.numeric(MF_Test$CHROM)
MH_Test$POS<-as.numeric(MH_Test$POS)

Mut_MH_Test <- inner_join(Mutect2_Test, MH_Test, by = c("CHROM", "POS", "REF", "ALT"))%>%
  separate(X1465_1024.pfc.bulk, into = c("GT", "AD", "AF", "DP", "F1R2", "F2R1", "FAD", "SB"), sep = ":", remove=FALSE)%>%
  separate(AD, into = c("AD1", "AD2"), sep = ",")%>%
  mutate(AD1 = as.numeric(AD1), AD2 = as.numeric(AD2), DP = as.numeric(DP))%>%
  mutate(AAF = (AD2 / DP) * 100)%>%mutate(Category="MH_Mut")%>%
  select("CHROM","POS","ID", "REF","ALT","QUAL","FILTER","INFO","FORMAT","AAF","Category")

Mut_MF_Test <- inner_join(Mutect2_Test, MF_Test, by = c("CHROM", "POS", "REF", "ALT"))%>%
  separate(X1465_1024.pfc.bulk, into = c("GT", "AD", "AF", "DP", "F1R2", "F2R1", "FAD", "SB"), sep = ":", remove=FALSE)%>%
  separate(AD, into = c("AD1", "AD2"), sep = ",")%>%
  mutate(AD1 = as.numeric(AD1), AD2 = as.numeric(AD2), DP = as.numeric(DP))%>%
  mutate(AAF = (AD2 / DP) * 100)%>%mutate(Category="MF_Mut")%>%
  select("CHROM","POS","ID", "REF","ALT","QUAL","FILTER","INFO","FORMAT","AAF","Category")

MH_MF_raw_Test <- inner_join(MF_Test, MH_Test, by = c("CHROM", "POS", "REF", "ALT"))%>%
  mutate(Category="MF_MH")%>%
  select("CHROM","POS", "REF","ALT","Category")%>%
  left_join(Mut_MF_Test,by = c("CHROM", "POS", "REF", "ALT"))

All_Test<-filter(MH_MF_raw_Test,Category.y=="MF_Mut")%>%
  mutate(Category="All")

MH_MF_Test<-anti_join(MH_MF_raw_Test,All_Test,by = c("CHROM", "POS", "REF", "ALT"))%>%
  mutate(Category="MH_MF")

Test_Results<-full_join(Mut_MH_Test,Mut_MF_Test,by = c("CHROM", "POS", "REF", "ALT"))%>%
              full_join(All_Test,by = c("CHROM", "POS", "REF", "ALT"))%>%
              full_join(MH_MF_Test,by = c("CHROM", "POS", "REF", "ALT"))
write_delim(Test_Results,file="M3_variants_Test.txt")
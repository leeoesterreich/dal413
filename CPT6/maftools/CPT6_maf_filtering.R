## ------------------------------------------------------------------ ##
##
## Script name: 
##
## Purpose of script:
##
## Author: Osama Shiraz Shah
##
## Date Created: 2021-05-08
##
## Copyright (c) Osama Shiraz Shah, 2021
## Email: shaho@upmc.edu
## Github: https://github.com/osamashiraz
## 
## ------------------------------------------------------------------ ##
##
## Notes: (YYYY-MM-DD -> Update)
##    
##
## ------------------------------------------------------------------ ##




## ****************************************************************** ##
##                            load essentials
## __________________________________________________________________ ##

`%notin%`=Negate(`%in%`)

library(magrittr)
library(ggplot2)
library(cowplot)
library(dplyr)
library(EnvStats)
library(ggpubr)
library(gridExtra)
library(ggthemes)
library(reshape2)


## ................................................................... ##





##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##
#......................1.CHUNK START                                ####
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##
#{



### LOAD VEP MAFs - FILTERED FOR DP > 15, AF > 15 AND VARIANTS IN CHR1-22, CHRX
mafFiles = list.files("/bgfs/alee/LO_LAB/Personal/Daisong/CPT6/maftools/filtered", pattern = "*.maf", full.names = T)

maf_list = lapply(mafFiles, function(x){
  read.delim(x, skip = 1)
})
names(maf_list) = gsub(x=mafFiles, pattern = "/bgfs/alee/LO_LAB/Personal/Daisong/CPT6/maftools/filtered/CPT6.maf/|.maf", replacement = "")

CPT6_MAF = do.call(rbind, maf_list)
dim(CPT6_MAF)
tumorIDs = names(table(CPT6_MAF$Tumor_Sample_Barcode))
newIDs = tumorIDs
newIDs = gsub(x=newIDs, pattern = "B6ILC-cells", replacement = "CPT6")
names(newIDs) = tumorIDs

CPT6_MAF$Tumor_Sample_Barcode = newIDs[CPT6_MAF$Tumor_Sample_Barcode]
table(CPT6_MAF$Tumor_Sample_Barcode) %>% sort() # basic filtering SNV count
#  CPT6
#  19201




### FILTERING PLAN

## 1. PON - GATK
# will use GATK PON HG38 for further filtering
# GATK PON https://gatk.broadinstitute.org/hc/en-us/artCPT6s/360035890811-Resource-bundle
#GATK_PON = as.data.frame(vroom::vroom("D:/TEMP_BOX/02_Datasets/Resources/GATK/1000g_pon.hg38.vcf"))
#head(GATK_PON)
#GATK_PON$uniqueID = paste0(GATK_PON$v1,GATK_PON$v2, GATK_PON$v4, GATK_PON$v5)

#CPT6_MAF$uniqueID = paste0(CPT6_MAF$Chromosome, CPT6_MAF$Start_Position, CPT6_MAF$Reference_Allele, CPT6_MAF$Tumor_Seq_Allele2)
#dim(CPT6_MAF) # 330443    123

#CPT6_MAF_filter1 = subset(CPT6_MAF, uniqueID %notin% GATK_PON$uniqueID)
#dim(CPT6_MAF_filter1) # 224515    123




## 2. gnomAD - remove common polymorphisms
#summary(subset(CPT6_MAF_filter1, FILTER %in% "common_variant")$gnomAD_AF)
#cutoff = min(subset(CPT6_MAF_filter1, FILTER %in% "common_variant")$gnomAD_AF) # 1.091e-05

#CPT6_MAF_filter1$gnomAD_AF[is.na(CPT6_MAF_filter1$gnomAD_AF)]=0

#CPT6_MAF_filter2 = subset(CPT6_MAF_filter1, gnomAD_AF < 0.0001)

#dim(CPT6_MAF_filter2) #167164    
#table(CPT6_MAF_filter2$Tumor_Sample_Barcode) %>% sort()


## 3. Supporting reads
CPT6_MAF_filter3 = subset(CPT6_MAF, t_alt_count > 10)

dim(CPT6_MAF_filter3) # 5513

# subset(CPT6_MAF_filter3[,c("Tumor_Sample_Barcode", "Hugo_Symbol", "Variant_Classification", "Variant_Type", "Exon_Number", "BIOTYPE","HGVSp_Short", "PolyPhen", "SIFT", "AF", "t_depth", "t_alt_count")], Hugo_Symbol %in% "COL14A1")
subset(CPT6_MAF_filter3[,c("Hugo_Symbol", "Variant_Classification", "Variant_Type", "Exon_Number", "BIOTYPE","HGVSp_Short", "PolyPhen", "SIFT", "AF", "t_depth", "t_alt_count","IMPACT")], Hugo_Symbol %in% "COL14A1")


table(CPT6_MAF_filter3$Tumor_Sample_Barcode)
# CPT6
# 5513



## 4. Remove non-protein coding genes

table(CPT6_MAF_filter3$BIOTYPE) %>% sort()
table(subset(CPT6_MAF_filter3, BIOTYPE %in% "processed_transcript")$Hugo_Symbol) %>% sort() #nonsense_mediated_decay
head(subset(CPT6_MAF_filter3, BIOTYPE %in% "processed_transcript"))
CPT6_MAF_filter4 = subset(CPT6_MAF_filter3, BIOTYPE %in% "protein_coding")
dim(CPT6_MAF_filter4)

table((CPT6_MAF_filter4$HGVSc==""))

#head(subset(CPT6_MAF_filter3, HGVSc == "")) "HGVSC" is codon change, it will be "" for non-protein coding genes
table(CPT6_MAF_filter3$BIOTYPE, CPT6_MAF_filter3$Tumor_Sample_Barcode)
head(subset(CPT6_MAF_filter4, HGVSc == "")) ## out gene region effects
table(subset(CPT6_MAF_filter4, HGVSc == "")$Hugo_Symbol) %>% sort() %>% tail(20) ## genes with mutation in regulatory regions




# 5. Remove mutations outside CDS
CPT6_MAF_filter5 = subset(CPT6_MAF_filter4, HGVSc %notin% "")

dim(CPT6_MAF_filter5)

table(CPT6_MAF_filter5$Tumor_Sample_Barcode)


CPT6_MAF_filter5$AF = CPT6_MAF_filter5$t_alt_count/CPT6_MAF_filter5$t_depth


plot(density(na.omit(CPT6_MAF_filter5$t_depth[CPT6_MAF_filter5$t_depth<100] %>% as.numeric())))
min(na.omit(CPT6_MAF_filter5$AF %>% as.numeric()))
plot(density(na.omit(CPT6_MAF_filter5$AF %>% as.numeric())))
boxplot((na.omit(CPT6_MAF_filter5$AF %>% as.numeric())))



## 6. REMOVE SAME EVENTS OCCURING IN MANY SAMPLES
## if same mutation observed in many samples remove these
#head(CPT6_MAF_filter5$uniqueID)
#sort(table(CPT6_MAF_filter5$uniqueID)) %>% tail()
#freq_in_pool = as.data.frame(table(CPT6_MAF_filter5$uniqueID))
#freq_in_pool$Var1 = as.character(freq_in_pool$Var1)
#rownames(freq_in_pool) = freq_in_pool$Var1
#head(freq_in_pool)

### ADD FREQ IN POOL
#for(i in 1:nrow(freq_in_pool)){
#   uniqueID = freq_in_pool[i, "Var1"]
#   CPT6_MAF_filter5[CPT6_MAF_filter5$uniqueID %in% uniqueID, "freq_in_pool"] = freq_in_pool[uniqueID, "Freq"]
#}

#table(CPT6_MAF_filter5$freq_in_pool)
#plot(density(CPT6_MAF_filter5$freq_in_pool))
#freq_mutations = unique(subset(CPT6_MAF_filter5, freq_in_pool %in% 5)$uniqueID)
#subset(CPT6_MAF_filter5, uniqueID %in% freq_mutations[10])[,c("Tumor_Sample_Barcode", "Hugo_Symbol","Reference_Allele","t_GT","Tumor_Seq_Allele2", "Variant_Classification", "Variant_Type", "Exon_Number", "BIOTYPE","HGVSp_Short", "PolyPhen", "SIFT", "AF", "t_depth", "t_alt_count")]

#(15*20)/100

#CPT6_MAF_filter6 = subset(CPT6_MAF_filter5, freq_in_pool <= 3)
#dim(CPT6_MAF_filter6)
#table(CPT6_MAF_filter6$Tumor_Sample_Barcode)

# BCK4      CAMA1    HCC2185    HCC2218     IPH926 MDAMB134VI   MDAMB330   MDAMB468     MPE600      OCUBM      SKBR3    SUM44PE   UACC3133 
# 2809       3085       4277       5521       6653       2970       4921       5875       3800       5922       2700       3613       3377 
# WCRC25     ZR7530 
# 4347       7507

#table(CPT6_MAF_filter6$Hugo_Symbol) %>% sort() %>% tail(20)


# 6. Remove mutations with AF < 0.15 and DP < 10
CPT6_MAF_filter7 = subset(CPT6_MAF_filter5, AF >= 0.1 & t_depth >= 10)
dim(CPT6_MAF_filter7)

table(CPT6_MAF_filter7$Tumor_Sample_Barcode) %>% sort()

table(CPT6_MAF_filter7$Chromosome) %>% sort()
table(CPT6_MAF_filter7$Variant_Classification) %>% sort()
table(CPT6_MAF_filter7$Variant_Type) %>% sort()
table(CPT6_MAF_filter7$VARIANT_CLASS) %>% sort()
table(CPT6_MAF_filter7$HGVSp) %>% sort() %>% tail()
table(CPT6_MAF_filter7$Hugo_Symbol) %>% sort() %>% tail(30)
table(subset(CPT6_MAF_filter7, Variant_Classification %notin% c("Silent","Intron","Missense_Mutation"))$Hugo_Symbol) %>% sort() %>% tail(50)
head(CPT6_MAF_filter7)


# VARIANTS WITH EFFECT ON PROTEIN CHANGE
CPT6_MAF_filter8 = subset(CPT6_MAF_filter7, HGVSp_Short %notin% "") # & t_ref_count != t_alt_count - also not heterozygous(intogen used this?)
dim(CPT6_MAF_filter8)

table(CPT6_MAF_filter8$Tumor_Sample_Barcode) %>% sort()
table(CPT6_MAF_filter8$Hugo_Symbol) %>% sort() %>% tail(40)
table(CPT6_MAF_filter8$EXON) %>% sort() %>% tail()

desired_columns <- c("Tumor_Sample_Barcode", "Hugo_Symbol", "Chromosome", "Variant_Classification", "Variant_Type", "Exon_Number", "BIOTYPE","HGVSp_Short", "PolyPhen", "SIFT", "AF", "t_depth", "t_alt_count","VARIANT_CLASS","IMPACT")
CPT6_MAF_filter8 <- CPT6_MAF_filter8[, intersect(desired_columns, colnames(CPT6_MAF_filter8))]
colnames(CPT6_MAF_filter8)

subset(CPT6_MAF_filter8[,c("Tumor_Sample_Barcode", "Hugo_Symbol", "Chromosome", "Variant_Classification", "Variant_Type", "Exon_Number", "BIOTYPE","HGVSp_Short", "PolyPhen", "SIFT", "AF", "t_depth", "t_alt_count","VARIANT_CLASS","IMPACT")], Hugo_Symbol %in% "Kras")

subset(CPT6_MAF_filter8[,c("Tumor_Sample_Barcode", "Hugo_Symbol", "Chromosome", "Variant_Classification", "Variant_Type", "Exon_Number", "BIOTYPE","HGVSp_Short", "PolyPhen", "SIFT", "AF", "t_depth", "t_alt_count","VARIANT_CLASS")], Hugo_Symbol %in% "Muc3a")
MUC3A_muts = subset(CPT6_MAF_filter8, Hugo_Symbol %in% "MUC3A")
table(MUC3A_muts$Tumor_Sample_Barcode, MUC3A_muts$HGVSp_Short)
table(MUC3A_muts$HGVSp_Short) %>% sort()

subset(CPT6_MAF_filter8[,c("Tumor_Sample_Barcode", "Hugo_Symbol", "Variant_Classification", "Chromosome", "Variant_Type", "Exon_Number", "BIOTYPE","HGVSp_Short", "PolyPhen", "SIFT", "AF", "t_depth", "t_alt_count","VARIANT_CLASS")], Hugo_Symbol %in% "Tp53")

subset(CPT6_MAF_filter8[,c("Tumor_Sample_Barcode", "Hugo_Symbol", "Variant_Classification", "Chromosome", "Variant_Type", "Exon_Number", "BIOTYPE","HGVSp_Short", "PolyPhen", "SIFT", "AF", "t_depth", "t_alt_count","VARIANT_CLASS")], Hugo_Symbol %in% "Pik3ca")

subset(CPT6_MAF_filter8[,c("Tumor_Sample_Barcode", "Hugo_Symbol", "Variant_Classification", "Chromosome", "Variant_Type", "Exon_Number", "BIOTYPE","HGVSp_Short", "PolyPhen", "SIFT", "AF", "t_depth", "t_alt_count","VARIANT_CLASS")], Hugo_Symbol %in% "Pten")

subset(CPT6_MAF_filter8[,c("Tumor_Sample_Barcode", "Hugo_Symbol", "Variant_Classification", "Chromosome", "Variant_Type", "Exon_Number", "BIOTYPE","HGVSp_Short", "PolyPhen", "SIFT", "AF", "t_depth", "t_alt_count","VARIANT_CLASS")], Hugo_Symbol %in% "Tbx3")



ClassDF = CPT6_MAF_filter8 %>% group_by(Tumor_Sample_Barcode, Variant_Classification) %>% summarise(n=n()) %>% as.data.frame()

ggplot(ClassDF, aes(x=reorder(Tumor_Sample_Barcode,n), y=n, fill=Variant_Classification)) + geom_bar(stat = "identity", position = "fill") + 
  ggsci::scale_fill_igv() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))


TypeDF = CPT6_MAF_filter8 %>% group_by(Tumor_Sample_Barcode, Variant_Type) %>% summarise(n=n()) %>% as.data.frame()

ggplot(TypeDF, aes(x=reorder(Tumor_Sample_Barcode,n), y=n, fill=Variant_Type)) + geom_bar(stat = "identity", position = "fill") + 
  ggsci::scale_fill_igv() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))


effectDF = CPT6_MAF_filter8 %>% group_by(Tumor_Sample_Barcode, gsub(x=PolyPhen,pattern = "0.|[[:digit:]]+|[()]", replacement = "")) %>% summarise(n=n()) %>% as.data.frame()
names(effectDF) = c("Tumor_Sample_Barcode", "PolyPhen","n")
head(effectDF)

ggplot(effectDF, aes(x=reorder(Tumor_Sample_Barcode,n), y=n, fill=PolyPhen)) + geom_bar(stat = "identity", position = "fill") + 
  ggsci::scale_fill_igv() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))


effectDF = CPT6_MAF_filter8 %>% group_by(Tumor_Sample_Barcode, IMPACT) %>% summarise(n=n()) %>% as.data.frame()
head(effectDF)

ggplot(effectDF, aes(x=reorder(Tumor_Sample_Barcode,n), y=n, fill=IMPACT)) + geom_bar(stat = "identity", position = "fill") + 
  ggsci::scale_fill_igv() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

table(CPT6_MAF_filter8$Tumor_Sample_Barcode)


save(CPT6_MAF_filter8, file = "/bgfs/alee/LO_LAB/Personal/Daisong/CPT6/maftools/filtered/basci_filtering.Rdata")
vroom::vroom_write(CPT6_MAF_filter8, path = "/bgfs/alee/LO_LAB/Personal/Daisong/CPT6/maftools/filtered/basci_filtering.tsv", delim = "\t", quote = "none")


library(maftools)

maf_file <- CPT6_MAF_filter8
laml.maf <- read.maf(maf = maf_file)

plotmafSummary(maf = laml.maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)

oncoplot(maf = laml.maf, top = 10)

lollipopPlot(
  maf = laml.maf,
  gene = 'KRAS',
  AACol = 'Amino_acids',
  showMutationRate = TRUE,
  labelPos = 882
)

library(EDASeq)
getGeneLengthAndGCContent(ensembl_list, "hsa")

head(CPT6_MAF_filter8)

gene_map = CPT6_MAF_filter8[,c("Hugo_Symbol")]
gene_map = na.omit(gene_map[!duplicated(gene_map$Gene)])
rownames(gene_map) = gene_map$Gene

geneMutation_rate = as.data.frame(table(CPT6_MAF_filter8$Gene))
rownames(geneMutation_rate) = geneMutation_rate$Var1
geneMutation_rate$Var1 = as.character(geneMutation_rate$Var1)
geneMutation_rate$gene = gene_map[geneMutation_rate$Var1, "Hugo_Symbol"]
geneMutation_rate = geneMutation_rate[!duplicated(geneMutation_rate$Var1),]

ensembl_list <- geneMutation_rate$Var1 %>% as.character()
geneLength = EDASeq::getGeneLengthAndGCContent(ensembl_list, "hsa", mode="org.db")
head(geneLength)


geneMutation_rate_df = cbind(geneMutation_rate, geneLength)
head(geneMutation_rate_df)

geneMutation_rate_df$log2Length = log2(geneMutation_rate_df$length+1)
ggplot(geneMutation_rate_df, aes(Freq, log2Length, color = gc)) + geom_point(size=5) + 
       ggrepel::geom_text_repel(data=subset(geneMutation_rate_df, Freq > 10),color="red", aes(label=gene))
ggplot(geneMutation_rate_df, aes(Freq, gc, color = log2Length)) + geom_point(size=5)









## my chunk 1 -------------------------------------------------------

## ..................................................................



## my chunk 1 -------------------------------------------------------

## ..................................................................



#} 
#......................# CHUNK END 
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##

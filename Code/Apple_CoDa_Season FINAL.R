#Two-year monitoring of Lm in apple packing houses
#Analysis of microbiomes and mycobiomes using CoDa approach at ASV level
#Laura Rolon
#Last updated: 07/26/22

#Load packages
library(ape)
library(vegan)
library(ggplot2)
library(phyloseq)
library(cowplot)
library(tidyr)
library(dplyr)
library(compositions)
library(zCompositions)
library(viridis)
library(readxl)
library(pairwiseAdonis)
library(psych)
library(svglite)

set.seed(336)

#Set working directory to where files are located
setwd("/storage/home/m/mlr355/work/Apple/Downstream")

#### IMPORT DATA ####

#Import data - 16s 
asvs_16s<-read.csv('ASV_16s.csv', header = TRUE, row.names = 1)
taxon_16s<-read.csv('Taxon_16s.csv', header = TRUE, row.names = 1)
metadata_16s<-read.csv('metadata_apple_16s.csv', header=TRUE, row.names = 1)


#Import ASV table - ITS
asvs_ITS<-read.csv('ASV_ITS.csv', header = TRUE, row.names = 1)
taxon_ITS<-read.csv('Taxon_ITS.csv', header = TRUE, row.names = 1)
metadata_ITS<-read.csv('metadata_apple_ITS.csv', header=TRUE, row.names = 1)

#Clean up Taxon tables
#Add '_unclassified' marker to NAs in taxon table
taxon_16s$Phylum<-ifelse(is.na(taxon_16s$Phylum), paste(taxon_16s$Kingdom, "unclassified", sep = '_'), taxon_16s$Phylum)
taxon_16s$Class<-ifelse(is.na(taxon_16s$Class), paste(taxon_16s$Phylum, "unclassified", sep = '_'), taxon_16s$Class)
taxon_16s$Order<-ifelse(is.na(taxon_16s$Order), paste(taxon_16s$Class, "unclassified", sep = '_'), taxon_16s$Order)
taxon_16s$Family<-ifelse(is.na(taxon_16s$Family), paste(taxon_16s$Order, "unclassified", sep = '_'), taxon_16s$Family)
taxon_16s$Genus<-ifelse(is.na(taxon_16s$Genus), paste(taxon_16s$Family, "unclassified", sep = '_'), taxon_16s$Genus)
taxon_16s$Species<-ifelse(is.na(taxon_16s$Species), paste(taxon_16s$Genus, "unclassified", sep = '_'), taxon_16s$Species)

taxon_ITS$Phylum<-ifelse(is.na(taxon_ITS$Phylum), paste(taxon_ITS$Kingdom, "unclassified", sep = '_'), taxon_ITS$Phylum)
taxon_ITS$Class<-ifelse(is.na(taxon_ITS$Class), paste(taxon_ITS$Phylum, "unclassified", sep = '_'), taxon_ITS$Class)
taxon_ITS$Order<-ifelse(is.na(taxon_ITS$Order), paste(taxon_ITS$Class, "unclassified", sep = '_'), taxon_ITS$Order)
taxon_ITS$Family<-ifelse(is.na(taxon_ITS$Family), paste(taxon_ITS$Order, "unclassified", sep = '_'), taxon_ITS$Family)
taxon_ITS$Genus<-ifelse(is.na(taxon_ITS$Genus), paste(taxon_ITS$Family, "unclassified", sep = '_'), taxon_ITS$Genus)
taxon_ITS$Species<-ifelse(is.na(taxon_ITS$Species), paste(taxon_ITS$Genus, "unclassified", sep = '_'), taxon_ITS$Species)

#Remove extra _unclassified
taxon_16s$Class<-gsub("_unclassified_unclassified", '_unclassified', taxon_16s$Class)
taxon_16s$Order<-gsub("_unclassified_unclassified_unclassified", '_unclassified', taxon_16s$Order)
taxon_16s$Order<-gsub("_unclassified_unclassified", '_unclassified', taxon_16s$Order)
taxon_16s$Family<-gsub("_unclassified_unclassified_unclassified_unclassified", '_unclassified', taxon_16s$Family)
taxon_16s$Family<-gsub("_unclassified_unclassified_unclassified", '_unclassified', taxon_16s$Family)
taxon_16s$Family<-gsub("_unclassified_unclassified", '_unclassified', taxon_16s$Family)
taxon_16s$Genus<-gsub("_unclassified_unclassified_unclassified_unclassified_unclassified", '_unclassified', taxon_16s$Genus)
taxon_16s$Genus<-gsub("_unclassified_unclassified_unclassified_unclassified", '_unclassified', taxon_16s$Genus)
taxon_16s$Genus<-gsub("_unclassified_unclassified_unclassified", '_unclassified', taxon_16s$Genus)
taxon_16s$Genus<-gsub("_unclassified_unclassified", '_unclassified', taxon_16s$Genus)
taxon_16s$Species<-gsub("_unclassified_unclassified_unclassified_unclassified_unclassified_unclassified", '_unclassified', taxon_16s$Species)
taxon_16s$Species<-gsub("_unclassified_unclassified_unclassified_unclassified_unclassified", '_unclassified', taxon_16s$Species)
taxon_16s$Species<-gsub("_unclassified_unclassified_unclassified_unclassified", '_unclassified', taxon_16s$Species)
taxon_16s$Species<-gsub("_unclassified_unclassified_unclassified", '_unclassified', taxon_16s$Species)
taxon_16s$Species<-gsub("_unclassified_unclassified", '_unclassified', taxon_16s$Species)

taxon_ITS$Class<-gsub("_unclassified_unclassified", '_unclassified', taxon_ITS$Class)
taxon_ITS$Order<-gsub("_unclassified_unclassified_unclassified", '_unclassified', taxon_ITS$Order)
taxon_ITS$Order<-gsub("_unclassified_unclassified", '_unclassified', taxon_ITS$Order)
taxon_ITS$Family<-gsub("_unclassified_unclassified_unclassified_unclassified", '_unclassified', taxon_ITS$Family)
taxon_ITS$Family<-gsub("_unclassified_unclassified_unclassified", '_unclassified', taxon_ITS$Family)
taxon_ITS$Family<-gsub("_unclassified_unclassified", '_unclassified', taxon_ITS$Family)
taxon_ITS$Genus<-gsub("_unclassified_unclassified_unclassified_unclassified_unclassified", '_unclassified', taxon_ITS$Genus)
taxon_ITS$Genus<-gsub("_unclassified_unclassified_unclassified_unclassified", '_unclassified', taxon_ITS$Genus)
taxon_ITS$Genus<-gsub("_unclassified_unclassified_unclassified", '_unclassified', taxon_ITS$Genus)
taxon_ITS$Genus<-gsub("_unclassified_unclassified", '_unclassified', taxon_ITS$Genus)
taxon_ITS$Species<-gsub("_unclassified_unclassified_unclassified_unclassified_unclassified_unclassified", '_unclassified', taxon_ITS$Species)
taxon_ITS$Species<-gsub("_unclassified_unclassified_unclassified_unclassified_unclassified", '_unclassified', taxon_ITS$Species)
taxon_ITS$Species<-gsub("_unclassified_unclassified_unclassified_unclassified", '_unclassified', taxon_ITS$Species)
taxon_ITS$Species<-gsub("_unclassified_unclassified_unclassified", '_unclassified', taxon_ITS$Species)
taxon_ITS$Species<-gsub("_unclassified_unclassified", '_unclassified', taxon_ITS$Species)

#Remove taxonomic rank from ITS taxon table
taxon_ITS$Kingdom<-gsub("k__","", taxon_ITS$Kingdom)
taxon_ITS$Phylum<-gsub("k__","", taxon_ITS$Phylum)
taxon_ITS$Phylum<-gsub("p__","", taxon_ITS$Phylum)
taxon_ITS$Class<-gsub("k__","", taxon_ITS$Class)
taxon_ITS$Class<-gsub("p__","", taxon_ITS$Class)
taxon_ITS$Class<-gsub("c__","", taxon_ITS$Class)
taxon_ITS$Order<-gsub("k__","", taxon_ITS$Order)
taxon_ITS$Order<-gsub("p__","", taxon_ITS$Order)
taxon_ITS$Order<-gsub("c__","", taxon_ITS$Order)
taxon_ITS$Order<-gsub("o__","", taxon_ITS$Order)
taxon_ITS$Family<-gsub("k__","", taxon_ITS$Family)
taxon_ITS$Family<-gsub("p__","", taxon_ITS$Family)
taxon_ITS$Family<-gsub("c__","", taxon_ITS$Family)
taxon_ITS$Family<-gsub("o__","", taxon_ITS$Family)
taxon_ITS$Family<-gsub("f__","", taxon_ITS$Family)
taxon_ITS$Genus<-gsub("k__","", taxon_ITS$Genus)
taxon_ITS$Genus<-gsub("p__","", taxon_ITS$Genus)
taxon_ITS$Genus<-gsub("c__","", taxon_ITS$Genus)
taxon_ITS$Genus<-gsub("o__","", taxon_ITS$Genus)
taxon_ITS$Genus<-gsub("f__","", taxon_ITS$Genus)
taxon_ITS$Genus<-gsub("g__","", taxon_ITS$Genus)
taxon_ITS$Species<-gsub("k__","", taxon_ITS$Species)
taxon_ITS$Species<-gsub("p__","", taxon_ITS$Species)
taxon_ITS$Species<-gsub("c__","", taxon_ITS$Species)
taxon_ITS$Species<-gsub("o__","", taxon_ITS$Species)
taxon_ITS$Species<-gsub("f__","", taxon_ITS$Species)
taxon_ITS$Species<-gsub("g__","", taxon_ITS$Species)
taxon_ITS$Species<-gsub("s__","", taxon_ITS$Species)


#Convert asv and taxon tables to matrix
asvs_16s<-as.matrix(asvs_16s)
taxon_16s<-as.matrix(taxon_16s)

asvs_ITS<-as.matrix(asvs_ITS)
taxon_ITS<-as.matrix(taxon_ITS)

#Make phyloseq object
ps_16s<-phyloseq(otu_table(asvs_16s, taxa_are_rows = FALSE), tax_table(taxon_16s), sample_data(metadata_16s))
ps_ITS<-phyloseq(otu_table(asvs_ITS, taxa_are_rows = FALSE), tax_table(taxon_ITS), sample_data(metadata_ITS))

#Remove Chloroplast and Mitochondria reads from ASV table
physeq_16s <- ps_16s %>%  subset_taxa( Order!="Chloroplast" | is.na(Order) )
physeq_16s <- physeq_16s %>% subset_taxa( Family!= "Mitochondria" | is.na("Family"))

#Filter samples that were extracted with Qiagen kit
ps.16s<-subset_samples(physeq_16s, Extraction == "Qiagen")
ps.ITS<-subset_samples(ps_ITS, Extraction == "Qiagen")

#Get ASV table from phyloseq object
asv_16s<-as.data.frame(t(otu_table(ps.16s)))
tail(rowSums(asv_16s))

asv_ITS<-as.data.frame(t(otu_table(ps.ITS)))
tail(rowSums(asv_ITS))

#Get Taxon table from phyloseq object
taxon.16s<-as.matrix(tax_table(ps.16s))
taxon.ITS<-as.matrix(tax_table(ps.ITS))

#Remove ASVs with zero counts in all samples
asv_16s<-asv_16s[ which(rowSums(asv_16s)>0),]
asv_16s<-t(asv_16s)

asv_ITS<-asv_ITS[ which(rowSums(asv_ITS)>0),]
asv_ITS<-t(asv_ITS)

#Get metadata 
metadata.16s<-subset(metadata_16s, Extraction == "Qiagen")
metadata.ITS<-subset(metadata_ITS, Extraction == "Qiagen")

#Save ASV, taxon, and metadata files to use in Core, Network, and Differential Abundance Analysis
write.csv(asv_16s, file="ASV_16s_clean.csv")
write.csv(asv_ITS, file="ASV_ITS_clean.csv")
write.csv(taxon.16s, file="Taxon_16s_clean.csv")
write.csv(taxon.ITS, file="Taxon_ITS_clean.csv")
write.csv(metadata.16s, file="metadata_16s_clean.csv")
write.csv(metadata.ITS, file="metadata_ITS_clean.csv")

#### COMPOSITIONAL ANALYSIS OF MICROBIOME AT ASV LEVEL ####
#Based on Microbiome Analysis in R. Chap 10.

#Step 1: Convert ASV table to appropriate format
#Following step requires samples on rows and ASVs in columns
head(asv_16s) 
head(asv_ITS)

#Step 2: Replace zero values before clr transformation. 
#Uses Count Zero Multiplicative method to replace zeros and outputs pseudo-counts
asv.n0_16s<-t(cmultRepl(asv_16s, label=0, method="CZM", output="p-counts")) #No. corrected values:  160752
asv.n0_ITS<-t(cmultRepl(asv_ITS, label=0, method="CZM", output="p-counts")) #No. corrected values:  24968 

#Note: Check the output to make sure there are no negative numbers. If samples or ASV are sparse, the CZM method will 
#add a negative number that interferes with the log normalization. If necessary use function below to convert negative values
#into positives
asv_n0_16s<-ifelse(asv.n0_16s < 0, asv.n0_16s*(-1), asv.n0_16s)
asv_n0_ITS<-ifelse(asv.n0_ITS < 0, asv.n0_ITS*(-1), asv.n0_ITS)

#output table needs to have samples in columns and ASVs in rows
head(asv_n0_16s) 
head(asv_n0_ITS)

#Step 3: Convert data to proportions
asv.n0_16s_prop<-apply(asv_n0_16s, 2, function(x) {x/sum(x)})
asv.n0_ITS_prop<-apply(asv_n0_ITS, 2, function(x) {x/sum(x)})

#Step 4: Perform abundance and sample filtering and deal sparsity
#Filter the data to remove all taxa with less than 0.00001% abundance in any sample
asv.n0_16s_prop_f<-asv_n0_16s[apply(asv.n0_16s_prop, 1, min) > 0.000001, ]
asv.n0_ITS_prop_f<-asv_n0_ITS[apply(asv.n0_ITS_prop, 1, min) > 0.000001, ]

#Check that samples are on columns and ASVs in rows
head(asv.n0_16s_prop_f) 
head(asv.n0_ITS_prop_f)

#Step 5: perform CLR transformation
asv.n0.clr_16s<-t(apply(asv.n0_16s_prop_f, 2, function(x){log(x)-mean(log(x))}))
asv.n0.clr_ITS<-t(apply(asv.n0_ITS_prop_f, 2, function(x){log(x)-mean(log(x))}))

#Check output table. Samples should be in rows and ASVs in columns
head(asv.n0.clr_16s) 
head(asv.n0.clr_ITS)

#Step 6: Perform Singular Value Decomposition (PCA)
pc.clr_16s<-prcomp(asv.n0.clr_16s) #Doesn't give error in R 4.1, but gives error in R 4.0
pc.clr_ITS<-prcomp(asv.n0.clr_ITS)

#library(ade4) #This library contains the PCA function that doesn't give an error. Alternatively use prcomp()
#pc.clr_16s<-dudi.pca(asv.n0.clr_16s, scannf= FALSE, nf=5)
#pc.clr_ITS<-dudi.pca(asv.n0.clr_ITS, scannf= FALSE, nf=5)

png("Screeplot - PCA2.png", width = 400, height = 300, units = 'px')
par(mar=c(2,2,2,2))
par(mfrow=c(1,2))
screeplot(pc.clr_16s, type='barplot', main="Bacteria")
screeplot(pc.clr_ITS, type='barplot', main="Fungi ")
dev.off()

#Calculate total variance of the data
mvar.clr_16s<-mvar(asv.n0.clr_16s)
mvar.clr_ITS<-mvar(asv.n0.clr_ITS)

#Display results - 16s 
row_16s<-rownames(asv.n0.clr_16s) #Make vector with sample names
pc_out_16s<-as.data.frame(pc.clr_16s$x[,1:2]) #Get PC1 and PC2 
pc_out_meta_16s<-as.data.frame(bind_cols(pc_out_16s,metadata.16s)) #Add metadata information
row.names(pc_out_meta_16s)<-row_16s #Add rownames to dataframe
pc_out_meta_16s$Facility<-as.factor(pc_out_meta_16s$Facility)
pc_out_meta_16s$L..monocytogenes<-as.factor(pc_out_meta_16s$L.monocytogenes)
pc_out_meta_16s$Year<-as.factor(pc_out_meta_16s$Year)


# Make PCA plot - First 2 axis- color by facility/Y and shape by L..mono
#Fig 2A
PCA_16s <- ggplot(pc_out_meta_16s, aes(x=PC1,y=PC2, color=Facility, shape=Year))+
  geom_point(size=3)+
  scale_shape_manual(values = c(15,12))+
  theme(legend.text=element_text(size=15,color='black'), legend.title= element_text(size=15, face='italic')) +
  theme(legend.position = 'right')+
  #theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  theme(axis.title = element_text(size=15,color='black')) +
  theme(panel.background = element_rect(fill='white', color = NA),
        plot.background = element_rect(fill = 'white',color = NA), 
        panel.border = element_rect(color='black',fill = NA,size=1))+
  scale_x_continuous(name = paste("PC1: ", round(pc.clr_16s$sdev[1]^2/mvar.clr_16s*100, digits=1), "%", sep="")) +
  scale_y_continuous(name = paste("PC2: ", round(pc.clr_16s$sdev[2]^2/mvar.clr_16s*100, digits=1), "%", sep="")) +
  ggtitle("Bacteria", subtitle = "PCA by Facility and Season")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=16))+
  scale_color_viridis_d(begin = 0.2, end = 0.8, option='inferno')
PCA_16s
ggsave("PCA_Bacteria_ASV2.png", plot =PCA_16s, device="png", width=6, height=5, units="in",dpi=600)
ggsave("PCA_Bacteria_ASV2.svg", plot =PCA_16s, device="svg", width=6, height=5, units="in",dpi=600)

#Display results -ITS
row_ITS<-rownames(asv.n0.clr_ITS) #Make vector with sample names
pc_out_ITS<-as.data.frame(pc.clr_ITS$x[,1:2]) #Get PC1 and PC2
pc_out_meta_ITS<-as.data.frame(bind_cols(pc_out_ITS,metadata.ITS)) #Add metadata information
row.names(pc_out_meta_ITS)<-row_ITS #Add rownames to dataframe
pc_out_meta_ITS$Facility<-as.factor(pc_out_meta_ITS$Facility)
pc_out_meta_ITS$L..monocytogenes<-as.factor(pc_out_meta_ITS$L.monocytogenes)
pc_out_meta_ITS$Year<-as.factor(pc_out_meta_ITS$Year)


# Make PCA plot - First 2 axis- color by facility/Y and shape by L..mono
#Fig 2D
PCA_ITS <- ggplot(pc_out_meta_ITS, aes(x=PC1,y=PC2, color=Facility, shape=Year))+
  geom_point(size=3)+
  scale_shape_manual(values = c(15,12))+
  theme(legend.text=element_text(size=15,color='black'), legend.title= element_text(size=13, face='italic')) +
  theme(legend.position = 'right')+
  #theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  theme(axis.title = element_text(size=15,color='black')) +
  theme(panel.background = element_rect(fill='white', color = NA),
        plot.background = element_rect(fill = 'white',color = NA),
        panel.border = element_rect(color='black',fill = NA,size=1))+
  scale_x_continuous(name = paste("PC1: ", round(pc.clr_ITS$sdev[1]^2/mvar.clr_ITS*100, digits=1), "%", sep="")) +
  scale_y_continuous(name = paste("PC2: ", round(pc.clr_ITS$sdev[2]^2/mvar.clr_ITS*100, digits=1), "%", sep="")) +
    ggtitle("Fungi", subtitle = "PCA by Facility and Season")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=16))+
  scale_color_viridis_d(begin = 0.2, end = 0.8, option='inferno')
PCA_ITS
ggsave("PCA_Fungi_ASV.png", plot =PCA_ITS, device="png", width=6, height=5, units="in",dpi=600)
ggsave("PCA_Fungi_ASV.svg", plot =PCA_ITS, device="svg", width=6, height=5, units="in",dpi=600)

# PERMANOVA #
#Calculate Aitchinson distance
dist_16s<-dist(asv.n0.clr_16s, method='euclidean')
dist_ITS<-dist(asv.n0.clr_ITS, method='euclidean')

#16s-
permanova_16s<-pairwise.adonis2(dist_16s~Year:Facility+Facility+Year, data=metadata.16s, perm = 999, p.adjust.m = 'bonferroni')
permanova_16s

#ITS
permanova_ITS<-pairwise.adonis2(dist_ITS~Year:Facility+Facility+Year, data=metadata.ITS, perm = 999, p.adjust.m = 'bonferroni')
permanova_ITS



##NOTE: Permanova showed significant interaction effect for Year:Facility for fungi and
##Continue by (1) splitting ASV table by Year to see Facility effect within each year and
##(2) splitting ASV table by Facility to see Seasonal effect within each Facility

#### Composition of microbiota by year - Bubble plots by facility and year #### 
#ASV level plots

#Transform sample counts into compositions
asv_n0.acomp_16s<-as.data.frame(acomp(t(asv_n0_16s)), total=1)
asv_n0.acomp_ITS<-as.data.frame(acomp(t(asv_n0_ITS)), total=1)

write.csv(asv_n0.acomp_16s, file = 'ASVtable_RA_16s.csv')
write.csv(asv_n0.acomp_ITS, file = 'ASVtable_RA_ITS.csv')


#Make Phyloseq object for each year
phyloseq_16s<-phyloseq(otu_table(asv_n0.acomp_16s, taxa_are_rows = FALSE), tax_table(taxon.16s), sample_data(metadata.16s))
phyloseq_ITS<-phyloseq(otu_table(asv_n0.acomp_ITS, taxa_are_rows = FALSE), tax_table(taxon.ITS), sample_data(metadata.ITS))


#Make long format table from Phyloseq object
asv_16s_long <- phyloseq_16s %>%  
  transform_sample_counts(function(x) {x * 100} ) %>% #Recalculates the relative abundance 
  psmelt() %>%  #Melts to long format
  arrange(desc(Abundance))

asv_ITS_long <- phyloseq_ITS %>%  
  transform_sample_counts(function(x) {x * 100} ) %>% #Recalculates the relative abundance 
  psmelt() %>%  #Melts to long format
  arrange(desc(Abundance))

#Save ASV composition as .csv
write.csv(asv_16s_long, "asv_16s.csv")
write.csv(asv_ITS_long, "asv_ITS.csv")

#Make vector with ASV names above 1 or 10% relative abundance in at least one sample
asv_16s_over5<-unique(c(asv_16s_long$OTU[which(asv_16s_long$Abundance >=0.5)]))
asv_ITS_over5<-unique(c(asv_ITS_long$OTU[which(asv_ITS_long$Abundance >=10)])) 

#Filter table to obtain only ASVs with over 1 or 10% in at least one sample
asv_16s_over5abund <- filter(asv_16s_long, OTU %in% asv_16s_over5)
asv_ITS_over5abund <- filter(asv_ITS_long, OTU %in% asv_ITS_over5)

### Calculate mean relative abundance by Facility each ASV
asv_16s_over5abund_mean<-asv_16s_over5abund%>%
  group_by(OTU, Facility, Genus, Year,SampleID, SampleOrder, Week)%>%
  summarize(Mean=mean(Abundance))%>%
  arrange(desc(Mean))

asv_ITS_over5abund_mean<-asv_ITS_over5abund%>%
  group_by(OTU, Facility, Genus, Year,SampleID, SampleOrder, Week)%>%
  summarize(Mean=mean(Abundance))%>%
  arrange(desc(Mean))

# #Bubble plot
# bubbleplot_ASV16s<-ggplot(asv_16s_over5abund_mean, aes(x=Facility, y=reorder(OTU,desc(OTU)), size=Mean, color=Facility))+
#   geom_point()+facet_grid(Genus~Year, scales = "free", space = 'free')+
#   scale_size(limits=c(0,1),name = "Relative abundance (%)")+
#   theme(axis.title = element_blank(), axis.text.x=element_text(color='black', size=10), 
#         axis.ticks=element_line(color='black'),
#         axis.text.y = element_text(color='black', size=10)) + 
#   guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, nrow = 1)) +
#   theme(legend.position = "bottom")+
#   theme(panel.background = element_rect(fill="grey99", color =NA),
#         plot.background = element_rect(fill="white", color =NA)) +
#   theme(strip.background= element_rect(fill="white", color = "black"), 
#         strip.text.y = element_text(size=10, angle=0, face = 'italic'), 
#         panel.border = element_rect(color="black", fill=NA))+
#   theme(legend.text=element_text(size=10,color='black'), legend.title= element_text(size=10)) +
#   ggtitle("Bacteria")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
#   scale_color_viridis_d(begin = 0.2, end = 0.8, option='inferno')+
#   theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))
# bubbleplot_otu16s
# ggsave("Bubbleplot_otu16s.png", plot=bubbleplot_otu16s, device="png", width=8, height=10, units="in", dpi=600)
# 
# bubbleplot_otuITS<-ggplot(asv_ITS_over5abund_mean, aes(x=Facility, y=reorder(OTU,desc(OTU)), size=Mean, color=Facility))+
#   geom_point()+facet_grid(Genus~Year, scales = "free", space = 'free')+
#   scale_size(limits=c(0,1),name = "Relative abundance (%)")+
#   theme(axis.title = element_blank(), axis.text.x=element_text(color='black', size=10), 
#         axis.ticks=element_line(color='black'),
#         axis.text.y = element_text(color='black', size=10)) + 
#   guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, nrow = 1)) +
#   theme(legend.position = "bottom")+
#   theme(panel.background = element_rect(fill="grey99", color =NA),
#         plot.background = element_rect(fill="white", color =NA)) +
#   theme(strip.background= element_rect(fill="white", color = "black"), 
#         strip.text.y = element_text(size=10, angle=0, face = 'italic'), 
#         panel.border = element_rect(color="black", fill=NA))+
#   theme(legend.text=element_text(size=10,color='black'), legend.title= element_text(size=10)) +
#   ggtitle("Bacteria")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
#   scale_color_viridis_d(begin = 0.2, end = 0.8, option='inferno')+
#   theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))
# bubbleplot_otuITS
# ggsave("Bubbleplot_otuITS.png", plot=bubbleplot_otuITS, device="png", width=8, height=10, units="in", dpi=600)


#Barplots
barplot_16s<-ggplot(asv_16s_over5abund_mean, aes(x=SampleOrder, y=Mean, fill=Genus))+
  geom_bar(stat='identity')+facet_grid(Year~Facility, scales = "free", space = 'free')+
  ylab("Relative abundance (%)")+
  theme(axis.title.x = element_blank(), axis.text.x=element_text(color='black', size=10),
       axis.ticks=element_line(color='black'),
       axis.text.y = element_text(color='black', size=10)) +
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, nrow = 10)) +
  theme(legend.position = "bottom")+
  theme(panel.background = element_rect(fill="grey99", color =NA),
        plot.background = element_rect(fill="white", color =NA)) +
  theme(strip.background= element_rect(fill="white", color = "black"), 
        strip.text.y = element_text(size=10, angle=0, face = 'italic'), 
        panel.border = element_rect(color="black", fill=NA))+
  theme(legend.text=element_text(size=10,color='black'), legend.title= element_text(size=10)) +
  ggtitle("Bacteria")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))
  #scale_color_viridis_d(begin = 0.2, end = 0.8, option='inferno')+
  #theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))
barplot_16s
ggsave("barplot_16s2.png", plot=barplot_16s, device="png", width=8, height=10, units="in", dpi=600)
ggsave("barplot_16s2.svg", plot=barplot_16s, device="svg", width=12, height=10, units="in", dpi=600)

barplot_ITS<-ggplot(asv_ITS_over5abund_mean, aes(x=SampleOrder, y=Mean, fill=Genus))+
  geom_bar(stat='identity')+facet_grid(Year~Facility, scales = "free", space = 'free')+
  ylab("Relative abundance (%)")+
  theme(axis.title.x = element_blank(), axis.text.x=element_text(color='black', size=10),
       axis.ticks=element_line(color='black'),
       axis.text.y = element_text(color='black', size=10)) +
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, nrow = 10)) +
  theme(legend.position = "bottom")+
  theme(panel.background = element_rect(fill="grey99", color =NA),
        plot.background = element_rect(fill="white", color =NA)) +
  theme(strip.background= element_rect(fill="white", color = "black"), 
        strip.text.y = element_text(size=10, angle=0, face = 'italic'), 
        panel.border = element_rect(color="black", fill=NA))+
  theme(legend.text=element_text(size=10,color='black'), legend.title= element_text(size=10)) +
  ggtitle("Bacteria")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
  #scale_color_viridis_d(begin = 0.2, end = 0.8, option='inferno')+
  theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))
barplot_ITS
ggsave("barplot_ITS.png", plot=barplot_ITS, device="png", width=8, height=10, units="in", dpi=600)

# #ITS
# #Fig S3-B
# bubbleplot_otuITS<-ggplot(otu_ITS_over5abund_mean, aes(x=Facility, y=reorder(OTU,desc(OTU)), size=Mean, color=Facility))+
#   geom_point()+facet_grid(Family_clean~Year, scales = "free", space = 'free')+
#   scale_size(limits=c(0,50), name = "Relative abundance (%)")+
#   theme(axis.title = element_blank(), axis.text.x=element_text(color='black', size=16), 
#         axis.ticks=element_line(color='black'),
#         axis.text.y = element_text(color='black', size=16)) + 
#   guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
#   theme(panel.background = element_rect(fill="grey99", color =NA),
#         plot.background = element_rect(fill="transparent", color =NA)) +
#   theme(strip.background= element_rect(fill="white", color = "black"), 
#         strip.text.y = element_text(size=16, angle=0, face = 'italic'), 
#         panel.border = element_rect(color="black", fill=NA))+
#   theme(legend.text=element_text(size=15,color='black'), legend.title= element_text(size=16)) +
#   ggtitle("Bacteria")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=18))+
#   scale_color_viridis_d(begin = 0.2, end = 0.8, option='inferno')+
#   theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))
# bubbleplot_otuITS
# ggsave("Bubbleplot_otuITS.svg", plot=bubbleplot_otuITS, device="svg", width=8, height=10, units="in", dpi=600)


#### Heatmap of microbiota by presence of Lm at ASV level by year ####


### Calculate mean relative abundance by Facility each ASV
asv_16s_over5abund_mean_Lm<-asv_16s_over5abund%>%
  group_by(OTU, Facility, Genus, Year, L.monocytogenes)%>%
  summarize(Mean=mean(Abundance))%>%
  arrange(desc(Mean))

asv_ITS_over5abund_mean_Lm<-asv_ITS_over5abund%>%
  group_by(OTU, Facility, Genus, Year, L.monocytogenes)%>%
  summarize(Mean=mean(Abundance))%>%
  arrange(desc(Mean))


#Heatmaps by facility at the ASV level
#Fig 5A
#16s
heatmap_asv16s<-ggplot(asv_16s_over5abund_mean_Lm, aes(x=L.monocytogenes, y=reorder(OTU,desc(OTU)), fill=Mean))+
  geom_tile(color='black')+facet_grid(Genus~Facility+Year, scales = "free", space = 'free')+
  scale_size(limits=c(0,50),name = "Relative abundance (%)")+
  theme(axis.title = element_blank(), axis.text.x=element_text(color='black', size=10), 
        axis.ticks=element_line(color='black'),
        axis.text.y = element_text(color='black', size=10)) + 
  theme(legend.position = "bottom")+
  theme(panel.background = element_rect(fill="grey99", color =NA),
        plot.background = element_rect(fill="white", color =NA)) +
  theme(strip.background= element_rect(fill="white", color = "black"), 
        strip.text.y = element_text(size=10, angle=0, face = 'italic'), 
        panel.border = element_rect(color="black", fill=NA))+
  theme(legend.text=element_text(size=10,color='black'), legend.title= element_text(size=10)) +
  ggtitle("Bacteria")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
  scale_fill_viridis_c(begin = 1, end = 0, option='inferno')
  #theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))
heatmap_asv16s
ggsave("Heatmap_asv16s2.svg", plot=heatmap_asv16s, device="svg", width=8, height=10, units="in", dpi=600)
ggsave("Heatmap_asv16s2.png", plot=heatmap_asv16s, device="png", width=8, height=10, units="in", dpi=600)

#ITS
#Fig 6A
heatmap_asvITS<-ggplot(asv_ITS_over5abund_mean_Lm, aes(x=L.monocytogenes, y=reorder(OTU,desc(OTU)), fill=Mean))+
  geom_tile(color='black')+facet_grid(Genus~Facility+Year, scales = "free", space = 'free')+
  scale_size(limits=c(0,50),name = "Relative abundance (%)")+
  theme(axis.title = element_blank(), axis.text.x=element_text(color='black', size=10), 
        axis.ticks=element_line(color='black'),
        axis.text.y = element_text(color='black', size=10)) + 
  theme(legend.position = "bottom")+
  theme(panel.background = element_rect(fill="grey99", color =NA),
        plot.background = element_rect(fill="white", color =NA)) +
  theme(strip.background= element_rect(fill="white", color = "black"), 
        strip.text.y = element_text(size=10, angle=0, face = 'italic'), 
        panel.border = element_rect(color="black", fill=NA))+
  theme(legend.text=element_text(size=10,color='black'), legend.title= element_text(size=10)) +
  ggtitle("Fungi")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=18))+
  scale_fill_viridis_c(begin = 1, end = 0, option='inferno')+
  theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))
heatmap_asvITS
ggsave("Heatmap_asvITS.svg", plot=heatmap_asvITS, device="svg", width=8, height=10, units="in", dpi=600)
ggsave("Heatmap_asvITS.png", plot=heatmap_asvITS, device="png", width=8, height=10, units="in", dpi=600)



#### Compositional analysis by year at the ASV level ####
#Filter samples that were extracted with Qiagen kit
ps.16sY1<-subset_samples(ps.16s, Year == "Y1")
ps.16sY2<-subset_samples(ps.16s, Year == "Y2")

ps.ITSY1<-subset_samples(ps.ITS, Year == "Y1")
ps.ITSY2<-subset_samples(ps.ITS, Year == "Y2")

#Get ASV table from phyloseq object
asv_16sY1<-as.data.frame(t(otu_table(ps.16sY1)))
asv_16sY2<-as.data.frame(t(otu_table(ps.16sY2)))

asv_ITSY1<-as.data.frame(t(otu_table(ps.ITSY1)))
asv_ITSY2<-as.data.frame(t(otu_table(ps.ITSY2)))

#Remove ASVs with zero counts in all samples
asv_16sY1<-asv_16sY1[ which(rowSums(asv_16sY1)>0),]
asv_16sY1<-t(asv_16sY1)

asv_16sY2<-asv_16sY2[ which(rowSums(asv_16sY2)>0),]
asv_16sY2<-t(asv_16sY2)

asv_ITSY1<-asv_ITSY1[ which(rowSums(asv_ITSY1)>0),]
asv_ITSY1<-t(asv_ITSY1)

asv_ITSY2<-asv_ITSY2[ which(rowSums(asv_ITSY2)>0),]
asv_ITSY2<-t(asv_ITSY2)

#Get metadata 
metadata.16sY1<-subset(metadata.16s, Year == "Y1")
metadata.16sY2<-subset(metadata.16s, Year == "Y2")

metadata.ITSY1<-subset(metadata.ITS, Year == "Y1")
metadata.ITSY2<-subset(metadata.ITS, Year == "Y2")

#Step 1: Convert ASV table to appropriate format
#Following step requires samples on rows and ASVs in columns
head(asv_16sY1)
head(asv_ITSY1)

head(asv_16sY2)
head(asv_ITSY2)

#Step 2: Replace zero values with a small non- zero value before clr transformation. 
#Uses Count Zero Multiplicative method to replace zeros and outputs pseudo-counts
asv.n0_16sY1<-t(cmultRepl(asv_16sY1, label=0, method="CZM", output="p-counts")) #No. adjusted imputations:  160752 
asv.n0_ITSY1<-t(cmultRepl(asv_ITSY1, label=0, method="CZM", output="p-counts")) #No. adjusted imputations:  12658 

asv.n0_16sY2<-t(cmultRepl(asv_16sY2, label=0, method="CZM", output="p-counts")) #No. adjusted imputations:  68530 
asv.n0_ITSY2<-t(cmultRepl(asv_ITSY2, label=0, method="CZM", output="p-counts")) #No. adjusted imputations:  873 

#Note: Check the output to make sure there are no negative numbers. If samples or asv are sparce, the CZM method will 
#add a negative number that interferes with the log normalization. If necessary use function elow to convert negative values
#into positives
asv_n0_16sY1<-ifelse(asv.n0_16sY1 < 0, asv.n0_16sY1*(-1), asv.n0_16sY1)
asv_n0_ITSY1<-ifelse(asv.n0_ITSY1 < 0, asv.n0_ITSY1*(-1), asv.n0_ITSY1)

asv_n0_16sY2<-ifelse(asv.n0_16sY2 < 0, asv.n0_16sY2*(-1), asv.n0_16sY2)
asv_n0_ITSY2<-ifelse(asv.n0_ITSY2 < 0, asv.n0_ITSY2*(-1), asv.n0_ITSY2)

#output table needs to have samples in columns and ASVs in rows
head(asv_n0_16sY1) 
head(asv_n0_ITSY1)

head(asv_n0_16sY2) 
head(asv_n0_ITSY2)

#Step 3: Convert data to proportions
asv.n0_16sY1_prop<-apply(asv_n0_16sY1, 2, function(x) {x/sum(x)})
asv.n0_ITSY1_prop<-apply(asv_n0_ITSY1, 2, function(x) {x/sum(x)})

asv.n0_16sY2_prop<-apply(asv_n0_16sY2, 2, function(x) {x/sum(x)})
asv.n0_ITSY2_prop<-apply(asv_n0_ITSY2, 2, function(x) {x/sum(x)})

head(asv.n0_16sY1_prop)
head(asv.n0_ITSY1_prop)

head(asv.n0_16sY2_prop)
head(asv.n0_ITSY2_prop)

#Step 4: Perform abundance and sample filtering and deal sparsity
#Filter the data to remove all taxa with less than 0.00001% abundance in any sample
asv.n0_16sY1_prop_f<-asv_n0_16sY1[apply(asv.n0_16sY1_prop, 1, min) > 0.0000001, ]
asv.n0_ITSY1_prop_f<-asv_n0_ITSY1[apply(asv.n0_ITSY1_prop, 1, min) > 0.0000001, ]

asv.n0_16sY2_prop_f<-asv_n0_16sY2[apply(asv.n0_16sY2_prop, 1, min) > 0.0000001, ]
asv.n0_ITSY2_prop_f<-asv_n0_ITSY2[apply(asv.n0_ITSY2_prop, 1, min) > 0.0000001, ]

#Check that samples are on columns and ASVs in rows
head(asv.n0_16sY1_prop_f) 
head(asv.n0_ITSY1_prop_f)

head(asv.n0_16sY2_prop_f)
head(asv.n0_ITSY2_prop_f)

#Step 5: perform CLR transformation
asv.n0.clr_16sY1<-t(apply(asv.n0_16sY1_prop_f, 2, function(x){log(x)-mean(log(x))}))
asv.n0.clr_ITSY1<-t(apply(asv.n0_ITSY1_prop_f, 2, function(x){log(x)-mean(log(x))}))

asv.n0.clr_16sY2<-t(apply(asv.n0_16sY2_prop_f, 2, function(x){log(x)-mean(log(x))}))
asv.n0.clr_ITSY2<-t(apply(asv.n0_ITSY2_prop_f, 2, function(x){log(x)-mean(log(x))}))

#Check output table. Samples should be in rows and ASVs in columns
head(asv.n0.clr_16sY1) 
head(asv.n0.clr_ITSY1)

head(asv.n0.clr_16sY2) 
head(asv.n0.clr_ITSY2)

#Step 6: Perform Singular Value Decomposition (PCA)
pc.clr_16sY1<-prcomp(asv.n0.clr_16sY1)
pc.clr_ITSY1<-prcomp(asv.n0.clr_ITSY1)

pc.clr_16sY2<-prcomp(asv.n0.clr_16sY2)
pc.clr_ITSY2<-prcomp(asv.n0.clr_ITSY2)

png("Screeplot - PCA by year .png", width = 1000, height = 500, units = 'px')
par(mar=c(2,2,2,2))
par(mfrow=c(2,2))
screeplot(pc.clr_16sY1, type='lines', main="Bacteria Y1")
screeplot(pc.clr_16sY2, type='lines', main="Bacteria Y2")
screeplot(pc.clr_ITSY1, type='lines', main="Fungi Y1")
screeplot(pc.clr_ITSY2, type='lines', main="Fungi Y2")
dev.off()

#Calculate total variance of the data
mvar.clr_16sY1<-mvar(asv.n0.clr_16sY1)
mvar.clr_ITSY1<-mvar(asv.n0.clr_ITSY1)

mvar.clr_16sY2<-mvar(asv.n0.clr_16sY2)
mvar.clr_ITSY2<-mvar(asv.n0.clr_ITSY2)

#Display results - 16sY1
row_16sY1<-rownames(asv.n0.clr_16sY1) #Make vector with sample names
pc_out_16sY1<-as.data.frame(pc.clr_16sY1$x[,1:2]) #Get PC1 and PC2 
pc_out_meta_16sY1<-as.data.frame(bind_cols(pc_out_16sY1,metadata.16sY1)) #Add metadata information
row.names(pc_out_meta_16sY1)<-row_16sY1 #Add rownames to dataframe
pc_out_meta_16sY1$Facility<-as.factor(pc_out_meta_16sY1$Facility)
pc_out_meta_16sY1$L..monocytogenes<-as.factor(pc_out_meta_16sY1$L.monocytogenes)

# Make PCA plot - First 2 axis- color by facility and shape by L..mono
#Fig 2B
PCA_16sY1<- ggplot(pc_out_meta_16sY1, aes(x=PC1,y=PC2, color=Facility, shape=L..monocytogenes))+
  geom_point(size=3)+
  theme(legend.text=element_text(size=15,color='black'), legend.title= element_text(size=15, face='italic')) +
  theme(legend.position = 'right')+
  #theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  theme(axis.title = element_text(size=15,color='black')) +
  theme(panel.background = element_rect(fill='white', color = NA),
        plot.background = element_rect(fill = 'white',color = NA), 
        panel.border = element_rect(color='black',fill = NA,size=1))+
  scale_x_continuous(name = paste("PC1: ", round(pc.clr_16sY1$sdev[1]^2/mvar.clr_16sY1*100, digits=1), "%", sep="")) +
  scale_y_continuous(name = paste("PC2: ", round(pc.clr_16sY1$sdev[2]^2/mvar.clr_16sY1*100, digits=1), "%", sep="")) +
  ggtitle("Bacteria - Y1", subtitle = "PCA by Facility and Lm presence")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=16))+
  scale_color_viridis_d(begin = 0.2, end = 0.8, option='inferno')
PCA_16sY1
ggsave("PCA_16sY1_2.png", plot =PCA_16sY1, device="png", width=6, height=5, units="in",dpi=600)
ggsave("PCA_16sY1_2.svg", plot =PCA_16sY1, device="svg", width=6, height=5, units="in",dpi=600)


#Display results - 16sY2 
row_16sY2<-rownames(asv.n0.clr_16sY2) #Make vector with sample names
pc_out_16sY2<-as.data.frame(pc.clr_16sY2$x[,1:2]) #Get PC1 and PC2 
pc_out_meta_16sY2<-as.data.frame(bind_cols(pc_out_16sY2,metadata.16sY2)) #Add metadata information
row.names(pc_out_meta_16sY2)<-row_16sY2 #Add rownames to dataframe
pc_out_meta_16sY2$Facility<-as.factor(pc_out_meta_16sY2$Facility)
pc_out_meta_16sY2$L..monocytogenes<-as.factor(pc_out_meta_16sY2$L.monocytogenes)

# Make PCA plot - First 2 axis- color by facility and shape by L..mono
#Fig 2C
PCA_16sY2 <- ggplot(pc_out_meta_16sY2, aes(x=PC1,y=PC2, color=Facility, shape=L..monocytogenes))+
  geom_point(size=3)+
  theme(legend.text=element_text(size=15,color='black'), legend.title= element_text(size=15, face='italic')) +
  theme(legend.position = 'right')+
  #theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  theme(axis.title = element_text(size=15,color='black')) +
  theme(panel.background = element_rect(fill='white', color = NA),
        plot.background = element_rect(fill = 'white',color = NA), 
        panel.border = element_rect(color='black',fill = NA,size=1))+
  scale_x_continuous(name = paste("PC1: ", round(pc.clr_16sY2$sdev[1]^2/mvar.clr_16sY2*100, digits=1), "%", sep="")) +
  scale_y_continuous(name = paste("PC2: ", round(pc.clr_16sY2$sdev[2]^2/mvar.clr_16sY2*100, digits=1), "%", sep="")) +
  ggtitle("Bacteria - Y2", subtitle = "PCA by Facility and Lm presence")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=16))+
  scale_color_viridis_d(begin = 0.2, end = 0.8, option='inferno')
PCA_16sY2
ggsave("PCA_16sY2_2.png", plot =PCA_16sY2, device="png", width=6, height=5, units="in",dpi=600)
ggsave("PCA_16sY2_2.svg", plot =PCA_16sY2, device="svg", width=6, height=5, units="in",dpi=600)

#Display results - ITSY1
row_ITSY1<-rownames(asv.n0.clr_ITSY1) #Make vector with sample names
pc_out_ITSY1<-as.data.frame(pc.clr_ITSY1$x[,1:2]) #Get PC1 and PC2 
pc_out_meta_ITSY1<-as.data.frame(bind_cols(pc_out_ITSY1,metadata.ITSY1)) #Add metadata information
row.names(pc_out_meta_ITSY1)<-row_ITSY1 #Add rownames to dataframe
pc_out_meta_ITSY1$Facility<-as.factor(pc_out_meta_ITSY1$Facility)
pc_out_meta_ITSY1$L..monocytogenes<-as.factor(pc_out_meta_ITSY1$L.monocytogenes)

# Make PCA plot - First 2 axis- color by facility and shape by L..mono
#Fig 2E
PCA_ITSY1 <- ggplot(pc_out_meta_ITSY1, aes(x=PC1,y=PC2, color=Facility, shape=L..monocytogenes))+
  geom_point(size=3)+
  theme(legend.text=element_text(size=15,color='black'), legend.title= element_text(size=15, face='italic')) +
  theme(legend.position = 'right')+
  theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  theme(axis.title = element_text(size=15,color='black')) +
  theme(panel.background = element_rect(fill='white', color = NA),
        plot.background = element_rect(fill = 'white',color = NA), 
        panel.border = element_rect(color='black',fill = NA,size=1))+
  scale_x_continuous(name = paste("PC1: ", round(pc.clr_ITSY1$sdev[1]^2/mvar.clr_ITSY1*100, digits=1), "%", sep="")) +
  scale_y_continuous(name = paste("PC2: ", round(pc.clr_ITSY1$sdev[2]^2/mvar.clr_ITSY1*100, digits=1), "%", sep="")) +
  ggtitle("Fungi - Y1", subtitle = "PCA by Facility and Lm presence")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=16))+
  scale_color_viridis_d(begin = 0.2, end = 0.8, option='inferno')
PCA_ITSY1
ggsave("PCA_ITSY1.png", plot =PCA_ITSY1, device="png", width=6, height=5, units="in",dpi=600)
ggsave("PCA_ITSY1.svg", plot =PCA_ITSY1, device="svg", width=6, height=5, units="in",dpi=600)


#Display results - ITSY2
row_ITSY2<-rownames(asv.n0.clr_ITSY2) #Make vector with sample names
pc_out_ITSY2<-as.data.frame(pc.clr_ITSY2$x[,1:2]) #Get PC1 and PC2 
pc_out_meta_ITSY2<-as.data.frame(bind_cols(pc_out_ITSY2,metadata.ITSY2)) #Add metadata information
row.names(pc_out_meta_ITSY2)<-row_ITSY2 #Add rownames to dataframe
pc_out_meta_ITSY2$Facility<-as.factor(pc_out_meta_ITSY2$Facility)
pc_out_meta_ITSY2$L..monocytogenes<-as.factor(pc_out_meta_ITSY2$L.monocytogenes)

# Make PCA plot - First 2 axis- color by facility and shape by L..mono
#Fig 2E
PCA_ITSY2 <- ggplot(pc_out_meta_ITSY2, aes(x=PC1,y=PC2, color=Facility, shape=L..monocytogenes))+
  geom_point(size=3)+
  theme(legend.text=element_text(size=15,color='black'), legend.title= element_text(size=15, face='italic')) +
  theme(legend.position = 'right')+
  theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  theme(axis.title = element_text(size=15,color='black')) +
  theme(panel.background = element_rect(fill='white', color = NA),
        plot.background = element_rect(fill = 'white',color = NA), 
        panel.border = element_rect(color='black',fill = NA,size=1))+
  scale_x_continuous(name = paste("PC1: ", round(pc.clr_ITSY2$sdev[1]^2/mvar.clr_ITSY2*100, digits=1), "%", sep="")) +
  scale_y_continuous(name = paste("PC2: ", round(pc.clr_ITSY2$sdev[2]^2/mvar.clr_ITSY2*100, digits=1), "%", sep="")) +
  ggtitle("Fungi - Y2", subtitle = "PCA by Facility and Lm presence")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=16))+
  scale_color_viridis_d(begin = 0.2, end = 0.8, option='inferno')
PCA_ITSY2
ggsave("PCA_ITSY2.png", plot =PCA_ITSY2, device="png", width=6, height=5, units="in",dpi=600)
ggsave("PCA_ITSY2.svg", plot =PCA_ITSY2, device="svg", width=6, height=5, units="in",dpi=600)


# PERMANOVA #
#Calculate Aitchinson distance
dist_16sY1<-dist(asv.n0.clr_16sY1, method='euclidean')
dist_16sY2<-dist(asv.n0.clr_16sY2, method='euclidean')
dist_ITSY1<-dist(asv.n0.clr_ITSY1, method='euclidean')
dist_ITSY2<-dist(asv.n0.clr_ITSY2, method='euclidean')

#Two-way anova by facility and Lmono presence
#16s
permanova_16sY1<-pairwise.adonis2(dist_16sY1~L.monocytogenes:Facility+Facility+L.monocytogenes, data=metadata.16sY1, perm = 999, p.adjust.m = 'bonferroni')
permanova_16sY1

permanova_16sY2<-pairwise.adonis2(dist_16sY2~L.monocytogenes:Facility+Facility+L.monocytogenes, data=metadata.16sY2, perm = 999, p.adjust.m = 'bonferroni')
permanova_16sY2

#ITS
permanova_ITSY1<-pairwise.adonis2(dist_ITSY1~L.monocytogenes:Facility+Facility+L.monocytogenes, data=metadata.ITSY1, perm = 999, p.adjust.m = 'bonferroni')
permanova_ITSY1

permanova_ITSY2<-pairwise.adonis2(dist_ITSY2~L.monocytogenes:Facility+Facility+L.monocytogenes, data=metadata.ITSY2, perm = 999, p.adjust.m = 'bonferroni')
permanova_ITSY2

#NOTE: For the 4 datasets, there was a significant difference by facilities. Run one-way anova to see where the differences are.

#Pairwize one-way anova by facility
pairwise.adonis(dist_16sY1, factors=metadata.16sY1$Facility, perm = 999, p.adjust.m = 'bonferroni')
pairwise.adonis(dist_16sY2, factors=metadata.16sY2$Facility, perm = 999, p.adjust.m = 'bonferroni')

pairwise.adonis(dist_ITSY1, factors=metadata.ITSY1$Facility, perm = 999, p.adjust.m = 'bonferroni')
pairwise.adonis(dist_ITSY2, factors=metadata.ITSY2$Facility, perm = 999, p.adjust.m = 'bonferroni')

#### Composition of microbiota for each year - Barplots facility and year #### 
#ASV level plots

#Transform sample counts into compositions
asv_n0.acomp_16sY1<-as.data.frame(acomp(t(asv_n0_16sY1)), total=1)
asv_n0.acomp_ITSY1<-as.data.frame(acomp(t(asv_n0_ITSY1)), total=1)
asv_n0.acomp_16sY2<-as.data.frame(acomp(t(asv_n0_16sY2)), total=1)
asv_n0.acomp_ITSY2<-as.data.frame(acomp(t(asv_n0_ITSY2)), total=1)

#Make Phyloseq object for each year
phyloseq_16sY1<-phyloseq(otu_table(asv_n0.acomp_16sY1, taxa_are_rows = FALSE), tax_table(taxon_16s), sample_data(metadata.16sY1))
phyloseq_ITSY1<-phyloseq(otu_table(asv_n0.acomp_ITSY1, taxa_are_rows = FALSE), tax_table(taxon_ITS), sample_data(metadata.ITSY1))
phyloseq_16sY2<-phyloseq(otu_table(asv_n0.acomp_16sY2, taxa_are_rows = FALSE), tax_table(taxon_16s), sample_data(metadata.16sY2))
phyloseq_ITSY2<-phyloseq(otu_table(asv_n0.acomp_ITSY2, taxa_are_rows = FALSE), tax_table(taxon_ITS), sample_data(metadata.ITSY2))

#Make long format table from Phyloseq object
asv_16sY1_long <- phyloseq_16sY1 %>%  
  transform_sample_counts(function(x) {x * 100} ) %>% #Recalculates the relative abundance 
  psmelt() %>%  #Melts to long format
  arrange(desc(Abundance))

asv_ITSY1_long <- phyloseq_ITSY1 %>%  
  transform_sample_counts(function(x) {x * 100} ) %>% #Recalculates the relative abundance 
  psmelt() %>%  #Melts to long format
  arrange(desc(Abundance))

asv_16sY2_long <- phyloseq_16sY2 %>%  
  transform_sample_counts(function(x) {x * 100} ) %>% #Recalculates the relative abundance 
  psmelt() %>%  #Melts to long format
  arrange(desc(Abundance))

asv_ITSY2_long <- phyloseq_ITSY2 %>%  
  transform_sample_counts(function(x) {x * 100} ) %>% #Recalculates the relative abundance 
  psmelt() %>%  #Melts to long format
  arrange(desc(Abundance))


#Make vector with ASV names above 5% relative abundance in at least one sample
asv_16sY1_over5<-unique(c(asv_16sY1_long$OTU[which(asv_16sY1_long$Abundance >=1)]))
asv_ITSY1_over5<-unique(c(asv_ITSY1_long$OTU[which(asv_ITSY1_long$Abundance >=10)])) 

asv_16sY2_over5<-unique(c(asv_16sY2_long$OTU[which(asv_16sY2_long$Abundance >=0.5)]))
asv_ITSY2_over5<-unique(c(asv_ITSY2_long$OTU[which(asv_ITSY2_long$Abundance >=10)])) 

#Filter table to obtain only ASVs with over 1% in at least one sample
asv_16sY1_over5abund <- filter(asv_16sY1_long, OTU %in% asv_16sY1_over5)
asv_ITSY1_over5abund <- filter(asv_ITSY1_long, OTU %in% asv_ITSY1_over5)

asv_16sY2_over5abund <- filter(asv_16sY2_long, OTU %in% asv_16sY2_over5)
asv_ITSY2_over5abund <- filter(asv_ITSY2_long, OTU %in% asv_ITSY2_over5)

### Calculate mean relative abundance by Facility each ASV
asv_16sY1_over5abund_mean<-asv_16sY1_over5abund%>%
  group_by(OTU, Facility, Genus, Year,SampleID, SampleOrder, Week)%>%
  summarize(Mean=mean(Abundance))%>%
  arrange(desc(Mean))

asv_ITSY1_over5abund_mean<-asv_ITSY1_over5abund%>%
  group_by(OTU, Facility, Genus, Year,SampleID, SampleOrder, Week)%>%
  summarize(Mean=mean(Abundance))%>%
  arrange(desc(Mean))

asv_16sY2_over5abund_mean<-asv_16sY2_over5abund%>%
  group_by(OTU, Facility, Genus, Year,SampleID, SampleOrder, Week)%>%
  summarize(Mean=mean(Abundance))%>%
  arrange(desc(Mean))

asv_ITSY2_over5abund_mean<-asv_ITSY2_over5abund%>%
  group_by(OTU, Facility, Genus, Year,SampleID, SampleOrder, Week)%>%
  summarize(Mean=mean(Abundance))%>%
  arrange(desc(Mean))

#Merge data sets
asv_mean_16s<-rbind(asv_16sY1_over5abund_mean,asv_16sY2_over5abund_mean)
asv_mean_ITS<-rbind(asv_ITSY1_over5abund_mean,asv_ITSY2_over5abund_mean)

#Barplots
barplot_16s2<-ggplot(asv_mean_16s, aes(x=SampleOrder, y=Mean, fill=Genus))+
  geom_bar(stat='identity')+facet_grid(Year~Facility, scales = "free", space = 'free')+
  ylab("Relative abundance (%)")+
  theme(axis.title.x = element_blank(), axis.text.x=element_text(color='black', size=10),
        axis.ticks=element_line(color='black'),
        axis.text.y = element_text(color='black', size=10)) +
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol = 4)) +
  theme(legend.position = "bottom")+
  theme(panel.background = element_rect(fill="grey99", color =NA),
        plot.background = element_rect(fill="white", color =NA)) +
  theme(strip.background= element_rect(fill="white", color = "black"), 
        strip.text.y = element_text(size=10, angle=0, face = 'italic'), 
        panel.border = element_rect(color="black", fill=NA))+
  theme(legend.text=element_text(size=10,color='black'), legend.title= element_text(size=10)) +
  ggtitle("Bacteria")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
  #scale_color_viridis_d(begin = 0.2, end = 0.8, option='inferno')+
  theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))
barplot_16s2
ggsave("barplot_16s2.png", plot=barplot_16s2, device="png", width=8, height=10, units="in", dpi=600)
ggsave("barplot_16s2.svg", plot=barplot_16s2, device="svg", width=8, height=10, units="in", dpi=600)

barplot_ITS2<-ggplot(asv_mean_ITS, aes(x=SampleOrder, y=Mean, fill=Genus))+
  geom_bar(stat='identity')+facet_grid(Year~Facility, scales = "free", space = 'free')+
  ylab("Relative abundance (%)")+
  theme(axis.title.x = element_blank(), axis.text.x=element_text(color='black', size=10),
        axis.ticks=element_line(color='black'),
        axis.text.y = element_text(color='black', size=10)) +
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol = 4)) +
  theme(legend.position = "bottom")+
  theme(panel.background = element_rect(fill="grey99", color =NA),
        plot.background = element_rect(fill="white", color =NA)) +
  theme(strip.background= element_rect(fill="white", color = "black"), 
        strip.text.y = element_text(size=10, angle=0, face = 'italic'), 
        panel.border = element_rect(color="black", fill=NA))+
  theme(legend.text=element_text(size=10,color='black'), legend.title= element_text(size=10)) +
  ggtitle("Bacteria")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
  #scale_color_viridis_d(begin = 0.2, end = 0.8, option='inferno')+
  theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))
barplot_ITS2
ggsave("barplot_ITS2.png", plot=barplot_ITS2, device="png", width=8, height=10, units="in", dpi=600)
ggsave("barplot_ITS2.svg", plot=barplot_ITS2, device="svg", width=8, height=10, units="in", dpi=600)


#Calculate the mean relative abundance of genera by year and facility
genus_16sY1_mean_fac<-asv_16sY1_long%>%
  group_by(Genus,Facility,SampleID)%>%
  summarize(Sum=sum(Abundance))%>%
  arrange(desc(Sum))

genus_16sY1_fac<-genus_16sY1_mean_fac%>%
  group_by(Genus,Facility)%>%
  summarize(Mean=mean(Sum))%>%
  arrange(desc(Mean))

genus_ITSY1_mean_fac<-asv_ITSY1_long%>%
  group_by(Genus,Facility,SampleID)%>%
  summarize(Sum=sum(Abundance))%>%
  arrange(desc(Sum))

genus_ITSY1_fac<-genus_ITSY1_mean_fac%>%
  group_by(Genus,Facility)%>%
  summarize(Mean=mean(Sum))%>%
  arrange(desc(Mean))

genus_16sY2_mean_fac<-asv_16sY2_long%>%
  group_by(Genus,Facility,SampleID)%>%
  summarize(Sum=sum(Abundance))%>%
  arrange(desc(Sum))

genus_16sY2_fac<-genus_16sY2_mean_fac%>%
  group_by(Genus,Facility)%>%
  summarize(Mean=mean(Sum))%>%
  arrange(desc(Mean))

genus_ITSY2_mean_fac<-asv_ITSY2_long%>%
  group_by(Genus,Facility,SampleID)%>%
  summarize(Sum=sum(Abundance))%>%
  arrange(desc(Sum))

genus_ITSY2_fac<-genus_ITSY2_mean_fac%>%
  group_by(Genus,Facility)%>%
  summarize(Mean=mean(Sum))%>%
  arrange(desc(Mean))

#Calculate the mean relative abundance of asv by facility throughout each season
asv_16sY1_mean_fac<-asv_16sY1_long%>%
  group_by(OTU,Facility,SampleID)%>%
  summarize(Sum=sum(Abundance))%>%
  arrange(desc(Sum))

asv_16sY1_fac<-asv_16sY1_mean_fac%>%
  group_by(OTU,Facility)%>%
  summarize(Mean=mean(Sum))%>%
  arrange(desc(Mean))

asv_ITSY1_mean_fac<-asv_ITSY1_long%>%
  group_by(OTU,Facility,SampleID)%>%
  summarize(Sum=sum(Abundance))%>%
  arrange(desc(Sum))

asv_ITSY1_fac<-asv_ITSY1_mean_fac%>%
  group_by(OTU,Facility)%>%
  summarize(Mean=mean(Sum))%>%
  arrange(desc(Mean))

asv_16sY2_mean_fac<-asv_16sY2_long%>%
  group_by(OTU,Facility,SampleID)%>%
  summarize(Sum=sum(Abundance))%>%
  arrange(desc(Sum))

asv_16sY2_fac<-asv_16sY2_mean_fac%>%
  group_by(OTU,Facility)%>%
  summarize(Mean=mean(Sum))%>%
  arrange(desc(Mean))

asv_ITSY2_mean_fac<-asv_ITSY2_long%>%
  group_by(OTU,Facility,SampleID)%>%
  summarize(Sum=sum(Abundance))%>%
  arrange(desc(Sum))

asv_ITSY2_fac<-asv_ITSY2_mean_fac%>%
  group_by(OTU,Facility)%>%
  summarize(Mean=mean(Sum))%>%
  arrange(desc(Mean))


#### Compositional analysis by facility at the ASV level ####
#Filter samples that were extracted with Qiagen kit
ps.16sF1<-subset_samples(ps.16s, Facility == "F1")
ps.16sF2<-subset_samples(ps.16s, Facility == "F2")
ps.16sF3<-subset_samples(ps.16s, Facility == "F3")

ps.ITSF1<-subset_samples(ps.ITS, Facility == "F1")
ps.ITSF2<-subset_samples(ps.ITS, Facility == "F2")
ps.ITSF3<-subset_samples(ps.ITS, Facility == "F3")

#Get ASV table from phyloseq object
asv_16sF1<-as.data.frame(t(otu_table(ps.16sF1)))
asv_16sF2<-as.data.frame(t(otu_table(ps.16sF2)))
asv_16sF3<-as.data.frame(t(otu_table(ps.16sF3)))

asv_ITSF1<-as.data.frame(t(otu_table(ps.ITSF1)))
asv_ITSF2<-as.data.frame(t(otu_table(ps.ITSF2)))
asv_ITSF3<-as.data.frame(t(otu_table(ps.ITSF3)))

#Remove ASVs with zero counts in all samples
asv_16sF1<-asv_16sF1[ which(rowSums(asv_16sF1)>0),]
asv_16sF2<-asv_16sF2[ which(rowSums(asv_16sF2)>0),]
asv_16sF3<-asv_16sF3[ which(rowSums(asv_16sF3)>0),]

asv_ITSF1<-asv_ITSF1[ which(rowSums(asv_ITSF1)>0),]
asv_ITSF2<-asv_ITSF2[ which(rowSums(asv_ITSF2)>0),]
asv_ITSF3<-asv_ITSF3[ which(rowSums(asv_ITSF3)>0),]

#Get metadata 
metadata.16sF1<-subset(metadata_16s, Facility == "F1")
metadata.16sF2<-subset(metadata_16s, Facility == "F2")
metadata.16sF3<-subset(metadata_16s, Facility == "F3")

metadata.ITSF1<-subset(metadata_ITS, Facility == "F1")
metadata.ITSF2<-subset(metadata_ITS, Facility == "F2")
metadata.ITSF3<-subset(metadata_ITS, Facility == "F3")

#Step 1: Convert ASV table to appropriate format
#Following step requires samples on rows and ASVs in columns
head(t(asv_16sF1))
head(t(asv_16sF2))
head(t(asv_16sF2))

head(t(asv_ITSF1))
head(t(asv_ITSF2))
head(t(asv_ITSF2))

#Step 2: Replace zero values with a small non- zero value before clr transformation. 
#Uses Count Zero Multiplicative method to replace zeros and outputs pseudo-counts
asv.n0_16sF1<-t(cmultRepl(t(asv_16sF1), label=0, method="CZM", output="p-counts")) #No. adjusted imputations:  160752 
asv.n0_16sF2<-t(cmultRepl(t(asv_16sF2), label=0, method="CZM", output="p-counts")) #No. adjusted imputations:  160752 
asv.n0_16sF3<-t(cmultRepl(t(asv_16sF3), label=0, method="CZM", output="p-counts")) #No. adjusted imputations:  160752 

asv.n0_ITSF1<-t(cmultRepl(t(asv_ITSF1), label=0, method="CZM", output="p-counts")) #No. adjusted imputations:  160752 
asv.n0_ITSF2<-t(cmultRepl(t(asv_ITSF2), label=0, method="CZM", output="p-counts")) #No. adjusted imputations:  160752 
asv.n0_ITSF3<-t(cmultRepl(t(asv_ITSF3), label=0, method="CZM", output="p-counts")) #No. adjusted imputations:  160752 

#Note: Check the output to make sure there are no negative numbers. If samples or asv are sparce, the CZM method will 
#add a negative number that interferes with the log normalization. If necessary use function elow to convert negative values
#into positives
asv_n0_16sF1<-ifelse(asv.n0_16sF1 < 0, asv.n0_16sF1*(-1), asv.n0_16sF1)
asv_n0_16sF2<-ifelse(asv.n0_16sF2 < 0, asv.n0_16sF2*(-1), asv.n0_16sF2)
asv_n0_16sF3<-ifelse(asv.n0_16sF3 < 0, asv.n0_16sF3*(-1), asv.n0_16sF3)

asv_n0_ITSF1<-ifelse(asv.n0_ITSF1 < 0, asv.n0_ITSF1*(-1), asv.n0_ITSF1)
asv_n0_ITSF2<-ifelse(asv.n0_ITSF2 < 0, asv.n0_ITSF2*(-1), asv.n0_ITSF2)
asv_n0_ITSF3<-ifelse(asv.n0_ITSF3 < 0, asv.n0_ITSF3*(-1), asv.n0_ITSF3)


#output table needs to have samples in columns and ASVs in rows
head(asv_n0_16sF1) 
head(asv_n0_16sF2) 
head(asv_n0_16sF3) 

head(asv_n0_ITSF1) 
head(asv_n0_ITSF2) 
head(asv_n0_ITSF3) 

#Step 3: Convert data to proportions
asv.n0_16sF1_prop<-apply(asv_n0_16sF1, 2, function(x) {x/sum(x)})
asv.n0_16sF2_prop<-apply(asv_n0_16sF2, 2, function(x) {x/sum(x)})
asv.n0_16sF3_prop<-apply(asv_n0_16sF3, 2, function(x) {x/sum(x)})

asv.n0_ITSF1_prop<-apply(asv_n0_ITSF1, 2, function(x) {x/sum(x)})
asv.n0_ITSF2_prop<-apply(asv_n0_ITSF2, 2, function(x) {x/sum(x)})
asv.n0_ITSF3_prop<-apply(asv_n0_ITSF3, 2, function(x) {x/sum(x)})

#Step 4: Perform abundance and sample filtering and deal sparsity
#Filter the data to remove all taxa with less than 0.00001% abundance in any sample
asv.n0_16sF1_prop_f<-asv_n0_16sF1[apply(asv.n0_16sF1_prop, 1, min) > 0.0000001, ]
asv.n0_16sF2_prop_f<-asv_n0_16sF2[apply(asv.n0_16sF2_prop, 1, min) > 0.0000001, ]
asv.n0_16sF3_prop_f<-asv_n0_16sF3[apply(asv.n0_16sF3_prop, 1, min) > 0.0000001, ]

asv.n0_ITSF1_prop_f<-asv_n0_ITSF1[apply(asv.n0_ITSF1_prop, 1, min) > 0.0000001, ]
asv.n0_ITSF2_prop_f<-asv_n0_ITSF2[apply(asv.n0_ITSF2_prop, 1, min) > 0.0000001, ]
asv.n0_ITSF3_prop_f<-asv_n0_ITSF3[apply(asv.n0_ITSF3_prop, 1, min) > 0.0000001, ]

#Check that samples are on columns and ASVs in rows
head(asv.n0_16sF1_prop_f) 
head(asv.n0_16sF2_prop_f) 
head(asv.n0_16sF3_prop_f) 

head(asv.n0_ITSF1_prop_f) 
head(asv.n0_ITSF2_prop_f) 
head(asv.n0_ITSF3_prop_f) 

#Step 5: perform CLR transformation
asv.n0.clr_16sF1<-t(apply(asv.n0_16sF1_prop_f, 2, function(x){log(x)-mean(log(x))}))
asv.n0.clr_16sF2<-t(apply(asv.n0_16sF2_prop_f, 2, function(x){log(x)-mean(log(x))}))
asv.n0.clr_16sF3<-t(apply(asv.n0_16sF3_prop_f, 2, function(x){log(x)-mean(log(x))}))

asv.n0.clr_ITSF1<-t(apply(asv.n0_ITSF1_prop_f, 2, function(x){log(x)-mean(log(x))}))
asv.n0.clr_ITSF2<-t(apply(asv.n0_ITSF2_prop_f, 2, function(x){log(x)-mean(log(x))}))
asv.n0.clr_ITSF3<-t(apply(asv.n0_ITSF3_prop_f, 2, function(x){log(x)-mean(log(x))}))


#Check output table. Samples should be in rows and ASVs in columns
head(asv.n0.clr_16sF1) 
head(asv.n0.clr_16sF2) 
head(asv.n0.clr_16sF3) 

head(asv.n0.clr_ITSF1) 
head(asv.n0.clr_ITSF2) 
head(asv.n0.clr_ITSF3) 

#Step 6: Perform Singular Value Decomposition (PCA)
pc.clr_16sF1<-prcomp(asv.n0.clr_16sF1)
pc.clr_16sF2<-prcomp(asv.n0.clr_16sF2)
pc.clr_16sF3<-prcomp(asv.n0.clr_16sF3)

pc.clr_ITSF1<-prcomp(asv.n0.clr_ITSF1)
pc.clr_ITSF2<-prcomp(asv.n0.clr_ITSF2)
pc.clr_ITSF3<-prcomp(asv.n0.clr_ITSF3)

png("Screeplot - PCA by facility .png", width = 1000, height = 500, units = 'px')
par(mar=c(2,2,2,2))
par(mfrow=c(3,2))
screeplot(pc.clr_16sF1, type='lines', main="Bacteria F1")
screeplot(pc.clr_16sF2, type='lines', main="Bacteria F2")
screeplot(pc.clr_16sF3, type='lines', main="Bacteria F3")
screeplot(pc.clr_ITSF1, type='lines', main="Fungi F1")
screeplot(pc.clr_ITSF2, type='lines', main="Fungi F2")
screeplot(pc.clr_ITSF3, type='lines', main="Fungi F3")
dev.off()

#Calculate total variance of the data
mvar.clr_16sF1<-mvar(asv.n0.clr_16sF1)
mvar.clr_16sF2<-mvar(asv.n0.clr_16sF2)
mvar.clr_16sF3<-mvar(asv.n0.clr_16sF3)

mvar.clr_ITSF1<-mvar(asv.n0.clr_ITSF1)
mvar.clr_ITSF2<-mvar(asv.n0.clr_ITSF2)
mvar.clr_ITSF3<-mvar(asv.n0.clr_ITSF3)

#Display results - 16s
row_16sF1<-rownames(asv.n0.clr_16sF1) #Make vector with sample names
pc_out_16sF1<-as.data.frame(pc.clr_16sF1$x[,1:2]) #Get PC1 and PC2 
pc_out_meta_16sF1<-as.data.frame(bind_cols(pc_out_16sF1,metadata.16sF1)) #Add metadata information
row.names(pc_out_meta_16sF1)<-row_16sF1 #Add rownames to dataframe
pc_out_meta_16sF1$Facility<-as.factor(pc_out_meta_16sF1$Facility)
pc_out_meta_16sF1$Year<-as.factor(pc_out_meta_16sF1$Year)

row_16sF2<-rownames(asv.n0.clr_16sF2) #Make vector with sample names
pc_out_16sF2<-as.data.frame(pc.clr_16sF2$x[,1:2]) #Get PC1 and PC2 
pc_out_meta_16sF2<-as.data.frame(bind_cols(pc_out_16sF2,metadata.16sF2)) #Add metadata information
row.names(pc_out_meta_16sF2)<-row_16sF2 #Add rownames to dataframe
pc_out_meta_16sF2$Facility<-as.factor(pc_out_meta_16sF2$Facility)
pc_out_meta_16sF2$Year<-as.factor(pc_out_meta_16sF2$Year)

row_16sF3<-rownames(asv.n0.clr_16sF3) #Make vector with sample names
pc_out_16sF3<-as.data.frame(pc.clr_16sF3$x[,1:2]) #Get PC1 and PC2 
pc_out_meta_16sF3<-as.data.frame(bind_cols(pc_out_16sF3,metadata.16sF3)) #Add metadata information
row.names(pc_out_meta_16sF3)<-row_16sF3 #Add rownames to dataframe
pc_out_meta_16sF3$Facility<-as.factor(pc_out_meta_16sF3$Facility)
pc_out_meta_16sF3$Year<-as.factor(pc_out_meta_16sF3$Year)

# Make PCA plot - First 2 axis- color by facility and shape by L..mono

PCA_16sF1<- ggplot(pc_out_meta_16sF1, aes(x=PC1,y=PC2, color=Facility, shape=Year))+
  geom_point(size=3)+
  scale_shape_manual(values = c(15,12))+
  theme(legend.text=element_text(size=15,color='black'), legend.title= element_text(size=15, face='italic')) +
  theme(legend.position = 'right')+
  theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  theme(axis.title = element_text(size=15,color='black')) +
  theme(panel.background = element_rect(fill='white', color = NA),
        plot.background = element_rect(fill = 'white',color = NA), 
        panel.border = element_rect(color='black',fill = NA,size=1))+
  scale_x_continuous(name = paste("PC1: ", round(pc.clr_16sF1$sdev[1]^2/mvar.clr_16sF1*100, digits=1), "%", sep="")) +
  scale_y_continuous(name = paste("PC2: ", round(pc.clr_16sF1$sdev[2]^2/mvar.clr_16sF1*100, digits=1), "%", sep="")) +
  ggtitle("Bacteria - F1", subtitle = "PCA by Year")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=16))+
  scale_color_manual(values = "#420A68FF")
PCA_16sF1

PCA_16sF2<- ggplot(pc_out_meta_16sF2, aes(x=PC1,y=PC2, color=Facility, shape=Year))+
  geom_point(size=3)+
  scale_shape_manual(values = c(15,12))+
  theme(legend.text=element_text(size=15,color='black'), legend.title= element_text(size=15, face='italic')) +
  theme(legend.position = 'right')+
  theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  theme(axis.title = element_text(size=15,color='black')) +
  theme(panel.background = element_rect(fill='white', color = NA),
        plot.background = element_rect(fill = 'white',color = NA), 
        panel.border = element_rect(color='black',fill = NA,size=1))+
  scale_x_continuous(name = paste("PC1: ", round(pc.clr_16sF2$sdev[1]^2/mvar.clr_16sF2*100, digits=1), "%", sep="")) +
  scale_y_continuous(name = paste("PC2: ", round(pc.clr_16sF2$sdev[2]^2/mvar.clr_16sF2*100, digits=1), "%", sep="")) +
  ggtitle("Bacteria - F2", subtitle = "PCA by Year")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=16))+
  scale_color_manual(values = "#BB3754FF")
PCA_16sF2

PCA_16sF3<- ggplot(pc_out_meta_16sF3, aes(x=PC1,y=PC2, color=Facility, shape=Year))+
  geom_point(size=3)+
  scale_shape_manual(values = c(15,12))+
  theme(legend.text=element_text(size=15,color='black'), legend.title= element_text(size=15, face='italic')) +
  theme(legend.position = 'right')+
  theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  theme(axis.title = element_text(size=15,color='black')) +
  theme(panel.background = element_rect(fill='white', color = NA),
        plot.background = element_rect(fill = 'white',color = NA), 
        panel.border = element_rect(color='black',fill = NA,size=1))+
  scale_x_continuous(name = paste("PC1: ", round(pc.clr_16sF3$sdev[1]^2/mvar.clr_16sF3*100, digits=1), "%", sep="")) +
  scale_y_continuous(name = paste("PC2: ", round(pc.clr_16sF3$sdev[2]^2/mvar.clr_16sF3*100, digits=1), "%", sep="")) +
  ggtitle("Bacteria - F3", subtitle = "PCA by Year")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=16))+
  scale_color_manual(values = "#FCA50AFF")
PCA_16sF3

#Display results - ITS
row_ITSF1<-rownames(asv.n0.clr_ITSF1) #Make vector with sample names
pc_out_ITSF1<-as.data.frame(pc.clr_ITSF1$x[,1:2]) #Get PC1 and PC2 
pc_out_meta_ITSF1<-as.data.frame(bind_cols(pc_out_ITSF1,metadata.ITSF1)) #Add metadata information
row.names(pc_out_meta_ITSF1)<-row_ITSF1 #Add rownames to dataframe
pc_out_meta_ITSF1$Facility<-as.factor(pc_out_meta_ITSF1$Facility)
pc_out_meta_ITSF1$Year<-as.factor(pc_out_meta_ITSF1$Year)

row_ITSF2<-rownames(asv.n0.clr_ITSF2) #Make vector with sample names
pc_out_ITSF2<-as.data.frame(pc.clr_ITSF2$x[,1:2]) #Get PC1 and PC2 
pc_out_meta_ITSF2<-as.data.frame(bind_cols(pc_out_ITSF2,metadata.ITSF2)) #Add metadata information
row.names(pc_out_meta_ITSF2)<-row_ITSF2 #Add rownames to dataframe
pc_out_meta_ITSF2$Facility<-as.factor(pc_out_meta_ITSF2$Facility)
pc_out_meta_ITSF2$Year<-as.factor(pc_out_meta_ITSF2$Year)

row_ITSF3<-rownames(asv.n0.clr_ITSF3) #Make vector with sample names
pc_out_ITSF3<-as.data.frame(pc.clr_ITSF3$x[,1:2]) #Get PC1 and PC2 
pc_out_meta_ITSF3<-as.data.frame(bind_cols(pc_out_ITSF3,metadata.ITSF3)) #Add metadata information
row.names(pc_out_meta_ITSF3)<-row_ITSF3 #Add rownames to dataframe
pc_out_meta_ITSF3$Facility<-as.factor(pc_out_meta_ITSF3$Facility)
pc_out_meta_ITSF3$Year<-as.factor(pc_out_meta_ITSF3$Year)

# Make PCA plot - First 2 axis- color by facility and shape by L..mono

PCA_ITSF1<- ggplot(pc_out_meta_ITSF1, aes(x=PC1,y=PC2, color=Facility, shape=Year))+
  geom_point(size=3)+
  scale_shape_manual(values = c(15,12))+
  theme(legend.text=element_text(size=15,color='black'), legend.title= element_text(size=15, face='italic')) +
  theme(legend.position = 'right')+
  theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  theme(axis.title = element_text(size=15,color='black')) +
  theme(panel.background = element_rect(fill='white', color = NA),
        plot.background = element_rect(fill = 'white',color = NA), 
        panel.border = element_rect(color='black',fill = NA,size=1))+
  scale_x_continuous(name = paste("PC1: ", round(pc.clr_ITSF1$sdev[1]^2/mvar.clr_ITSF1*100, digits=1), "%", sep="")) +
  scale_y_continuous(name = paste("PC2: ", round(pc.clr_ITSF1$sdev[2]^2/mvar.clr_ITSF1*100, digits=1), "%", sep="")) +
  ggtitle("Fungi - F1", subtitle = "PCA by Year")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=16))+
  scale_color_manual(values = "#420A68FF")
PCA_ITSF1

PCA_ITSF2<- ggplot(pc_out_meta_ITSF2, aes(x=PC1,y=PC2, color=Facility, shape=Year))+
  geom_point(size=3)+
  scale_shape_manual(values = c(15,12))+
  theme(legend.text=element_text(size=15,color='black'), legend.title= element_text(size=15, face='italic')) +
  theme(legend.position = 'right')+
  theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  theme(axis.title = element_text(size=15,color='black')) +
  theme(panel.background = element_rect(fill='white', color = NA),
        plot.background = element_rect(fill = 'white',color = NA), 
        panel.border = element_rect(color='black',fill = NA,size=1))+
  scale_x_continuous(name = paste("PC1: ", round(pc.clr_ITSF2$sdev[1]^2/mvar.clr_ITSF2*100, digits=1), "%", sep="")) +
  scale_y_continuous(name = paste("PC2: ", round(pc.clr_ITSF2$sdev[2]^2/mvar.clr_ITSF2*100, digits=1), "%", sep="")) +
  ggtitle("Fungi - F2", subtitle = "PCA by Year")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=16))+
  scale_color_manual(values = "#BB3754FF")
PCA_ITSF2

PCA_ITSF3<- ggplot(pc_out_meta_ITSF3, aes(x=PC1,y=PC2, color=Facility, shape=Year))+
  geom_point(size=3)+
  scale_shape_manual(values = c(15,12))+
  theme(legend.text=element_text(size=15,color='black'), legend.title= element_text(size=15, face='italic')) +
  theme(legend.position = 'right')+
  theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  theme(axis.title = element_text(size=15,color='black')) +
  theme(panel.background = element_rect(fill='white', color = NA),
        plot.background = element_rect(fill = 'white',color = NA), 
        panel.border = element_rect(color='black',fill = NA,size=1))+
  scale_x_continuous(name = paste("PC1: ", round(pc.clr_ITSF3$sdev[1]^2/mvar.clr_ITSF3*100, digits=1), "%", sep="")) +
  scale_y_continuous(name = paste("PC2: ", round(pc.clr_ITSF3$sdev[2]^2/mvar.clr_ITSF3*100, digits=1), "%", sep="")) +
  ggtitle("Fungi - F3", subtitle = "PCA by Year")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=16))+
  scale_color_manual(values = "#FCA50AFF")
PCA_ITSF3

#Combine plots
PCA_Fac = plot_grid(PCA_16sF1, PCA_16sF2, PCA_16sF3,
                    PCA_ITSF1, PCA_ITSF2, PCA_ITSF3,
                    ncol=3, nrow=3, labels = c("A","B","C","D","E", "F"), label_size = 20, vjust = 2, hjust = -1.5)


ggsave("PCA_Facility.png", plot =PCA_Fac, device="png", width=10, height=12, units="in",dpi=600)
ggsave("PCA_ITSY2.svg", plot =PCA_ITSY2, device="svg", width=6, height=5, units="in",dpi=600)


# PERMANOVA #
#Calculate Aitchinson distance
dist_16sF1<-dist(asv.n0.clr_16sF1, method='euclidean')
dist_16sF2<-dist(asv.n0.clr_16sF2, method='euclidean')
dist_16sF3<-dist(asv.n0.clr_16sF3, method='euclidean')

dist_ITSF1<-dist(asv.n0.clr_ITSF1, method='euclidean')
dist_ITSF2<-dist(asv.n0.clr_ITSF2, method='euclidean')
dist_ITSF3<-dist(asv.n0.clr_ITSF3, method='euclidean')

#One-way anova by facility and Lmono presence
#16s

#Pairwize one-way anova by year
pairwise.adonis(dist_16sF1, factors=metadata.16sF1$Year, perm = 999, p.adjust.m = 'bonferroni')
pairwise.adonis(dist_16sF2, factors=metadata.16sF2$Year, perm = 999, p.adjust.m = 'bonferroni')
pairwise.adonis(dist_16sF3, factors=metadata.16sF3$Year, perm = 999, p.adjust.m = 'bonferroni')

pairwise.adonis(dist_ITSF1, factors=metadata.ITSF1$Year, perm = 999, p.adjust.m = 'bonferroni')
pairwise.adonis(dist_ITSF2, factors=metadata.ITSF2$Year, perm = 999, p.adjust.m = 'bonferroni')
pairwise.adonis(dist_ITSF3, factors=metadata.ITSF3$Year, perm = 999, p.adjust.m = 'bonferroni')


#Check # of Listeria spp reads in samples that were used for Nanopore sequencing
asv_Nanopore_Lm<- asvs_16s[rownames(asvs_16s) %in% c("SRR12559261","SRR12559222","SRR12559271"), colnames(asv_Nanopore) %in% c("ASV11509", "ASV34463", "ASV34514")] 

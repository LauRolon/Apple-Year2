#Two-year monitoring of Lm in apple packing houses
#Analysis of core microbiomes and networks
#Laura Rolon
#Last updated: 08/10/21

#Attach libraries
library(phyloseq)
library(ape)
library(zCompositions)
library(compositions)
library(dplyr)
library(ggplot2)
library(cowplot)
library(svglite)

#Set working directory to where files are located
setwd()

#### IMPORT DATA ####

#Import OTU table - 16s (SILVA version 132)
otus_16s<-as.data.frame(import_mothur(mothur_shared_file = 'apple16s_132.shared'))

#Import taxonomy table -16s (SILVA version 132)
taxon_16s <- as.data.frame(import_mothur(mothur_constaxonomy_file = 'apple16s_132.cons.taxonomy'))
colnames(taxon_16s) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus")
taxon_16s$otu_id<-rownames(taxon_16s)
TAX_16s =tax_table(as.matrix(taxon_16s))

#Import metadata -16s
metadata_16s <-read.csv("metadata_apple_16s.csv", header=TRUE, row.names=1)

META_16s = sample_data(metadata_16s)


#Import OTU table - ITS
otus_ITS<-as.data.frame(import_mothur(mothur_shared_file = 'appleITS.shared'))

#Import taxonomy table -ITSY1
taxon_ITS <- as.data.frame(import_mothur(mothur_constaxonomy_file = 'appleITS.cons.taxonomy'))
colnames(taxon_ITS) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus")
TAX_ITS =tax_table(as.matrix(taxon_ITS))

#Import metadata -ITS
metadata_ITS <-read.csv("metadata_apple_ITS.csv", header=TRUE, row.names=1)
META_ITS = sample_data(metadata_ITS)

#### Prepare data for analyses ####
#Make phyloseq
OTU_16s <- otu_table(otus_16s, taxa_are_rows = TRUE)
phyloseq_16s = phyloseq(OTU_16s, TAX_16s, META_16s)
TREE_16s = rtree(ntaxa(phyloseq_16s), rooted=TRUE, tip.label = taxa_names(phyloseq_16s))  
phyloseq16s <- phyloseq(OTU_16s, TAX_16s, TREE_16s, META_16s)

OTU_ITS <- otu_table(otus_ITS, taxa_are_rows = TRUE)
phyloseq_ITS = phyloseq(OTU_ITS, TAX_ITS, META_ITS)
TREE_ITS = rtree(ntaxa(phyloseq_ITS), rooted=TRUE, tip.label = taxa_names(phyloseq_ITS))
phyloseqITS <- phyloseq(OTU_ITS, TAX_ITS, TREE_ITS, META_ITS)


# Subset Phyloseq for each facility
physeq_16sF1 <- subset_samples(phyloseq16s, Facility == "F1") 
physeq_16sF2 <- subset_samples(phyloseq16s, Facility == "F2") 
physeq_16sF3 <- subset_samples(phyloseq16s, Facility == "F3") 

physeq_ITSF1 <- subset_samples(phyloseqITS, Facility == "F1") 
physeq_ITSF2 <- subset_samples(phyloseqITS, Facility == "F2") 
physeq_ITSF3 <- subset_samples(phyloseqITS, Facility == "F3") 

#Export OTU table for each facility from Phyloseq 
otus_16sF1<-as.data.frame(as(otu_table(physeq_16sF1), "matrix"))
otus_16sF2<-as.data.frame(as(otu_table(physeq_16sF2), "matrix"))
otus_16sF3<-as.data.frame(as(otu_table(physeq_16sF3), "matrix"))

otus_ITSF1<-as.data.frame(as(otu_table(physeq_ITSF1), "matrix"))
otus_ITSF2<-as.data.frame(as(otu_table(physeq_ITSF2), "matrix"))
otus_ITSF3<-as.data.frame(as(otu_table(physeq_ITSF3), "matrix"))


#Subset metadata
metadata_16sF1<-subset(metadata_16s, Facility=="F1")
metadata_16sF2<-subset(metadata_16s, Facility=="F2")
metadata_16sF3<-subset(metadata_16s, Facility=="F3")

metadata_ITSF1<-subset(metadata_ITS, Facility=="F1")
metadata_ITSF2<-subset(metadata_ITS, Facility=="F2")
metadata_ITSF3<-subset(metadata_ITS, Facility=="F3")

#Remove OTUs that have count zero in all samples - Necessary step for the zero imputation function
otus_16sF1<-otus_16sF1[ which(rowSums(otus_16sF1)>0),]
otus_16sF2<-otus_16sF2[ which(rowSums(otus_16sF2)>0),]
otus_16sF3<-otus_16sF3[ which(rowSums(otus_16sF3)>0),]

otus_ITSF1<-otus_ITSF1[ which(rowSums(otus_ITSF1)>0),]
otus_ITSF2<-otus_ITSF2[ which(rowSums(otus_ITSF2)>0),]
otus_ITSF3<-otus_ITSF3[ which(rowSums(otus_ITSF3)>0),]

#### Occupancy abundance curves for whole data set ####
#Individual sampling sites are considered independent.
#Calculate presence absence
otus_16sF1_PA <- 1*((otus_16sF1>0)==1) 
otus_16sF2_PA <- 1*((otus_16sF2>0)==1) 
otus_16sF3_PA <- 1*((otus_16sF3>0)==1) 

otus_ITSF1_PA <- 1*((otus_ITSF1>0)==1) 
otus_ITSF2_PA <- 1*((otus_ITSF2>0)==1) 
otus_ITSF3_PA <- 1*((otus_ITSF3>0)==1) 

#Calculate mean occupancy for each OTU
otus_16sF1_occ <- as.data.frame(rowSums(otus_16sF1_PA)/ncol(otus_16sF1_PA))
otus_16sF2_occ <- as.data.frame(rowSums(otus_16sF2_PA)/ncol(otus_16sF2_PA))
otus_16sF3_occ <- as.data.frame(rowSums(otus_16sF3_PA)/ncol(otus_16sF3_PA))

otus_ITSF1_occ <- as.data.frame(rowSums(otus_ITSF1_PA)/ncol(otus_ITSF1_PA))
otus_ITSF2_occ <- as.data.frame(rowSums(otus_ITSF2_PA)/ncol(otus_ITSF2_PA))
otus_ITSF3_occ <- as.data.frame(rowSums(otus_ITSF3_PA)/ncol(otus_ITSF3_PA))

#Rename column
otus_16sF1_occ <-otus_16sF1_occ %>% rename(MeanOcc=`rowSums(otus_16sF1_PA)/ncol(otus_16sF1_PA)`)
otus_16sF2_occ <-otus_16sF2_occ %>% rename(MeanOcc=`rowSums(otus_16sF2_PA)/ncol(otus_16sF2_PA)`)
otus_16sF3_occ <-otus_16sF3_occ %>% rename(MeanOcc=`rowSums(otus_16sF3_PA)/ncol(otus_16sF3_PA)`)

otus_ITSF1_occ <-otus_ITSF1_occ %>% rename(MeanOcc=`rowSums(otus_ITSF1_PA)/ncol(otus_ITSF1_PA)`)
otus_ITSF2_occ <-otus_ITSF2_occ %>% rename(MeanOcc=`rowSums(otus_ITSF2_PA)/ncol(otus_ITSF2_PA)`)
otus_ITSF3_occ <-otus_ITSF3_occ %>% rename(MeanOcc=`rowSums(otus_ITSF3_PA)/ncol(otus_ITSF3_PA)`)

#Calculate relative abundance for each OTU by Facility
#First replace zero values with a small non- zero value before clr transformation. 
#Uses Count Zero Multiplicative method to replace zeros and outputs pseudo-counts
otu.n0_16sF1<-t(cmultRepl(t(otus_16sF1), label=0, method="CZM", output="p-counts")) 
otu.n0_16sF2<-t(cmultRepl(t(otus_16sF2), label=0, method="CZM", output="p-counts")) 
otu.n0_16sF3<-t(cmultRepl(t(otus_16sF3), label=0, method="CZM", output="p-counts")) 

otu.n0_ITSF1<-t(cmultRepl(t(otus_ITSF1), label=0, method="CZM", output="p-counts")) 
otu.n0_ITSF2<-t(cmultRepl(t(otus_ITSF2), label=0, method="CZM", output="p-counts")) 
otu.n0_ITSF3<-t(cmultRepl(t(otus_ITSF3), label=0, method="CZM", output="p-counts")) 

#Then convert to compositions with Aitchinson 
otu.n0.acomp_16sF1<-as.data.frame(acomp(t(otu.n0_16sF1)), total=1)
otu.n0.acomp_16sF2<-as.data.frame(acomp(t(otu.n0_16sF2)), total=1)
otu.n0.acomp_16sF3<-as.data.frame(acomp(t(otu.n0_16sF3)), total=1)

otu.n0.acomp_ITSF1<-as.data.frame(acomp(t(otu.n0_ITSF1)), total=1)
otu.n0.acomp_ITSF2<-as.data.frame(acomp(t(otu.n0_ITSF2)), total=1)
otu.n0.acomp_ITSF3<-as.data.frame(acomp(t(otu.n0_ITSF3)), total=1)


#Calculate the mean relative abundance for each OTU in each facility
otus_16sF1_abun <- as.data.frame(rowSums(t(otu.n0.acomp_16sF1))/ncol(t(otu.n0.acomp_16sF1)))
otus_16sF2_abun <- as.data.frame(rowSums(t(otu.n0.acomp_16sF2))/ncol(t(otu.n0.acomp_16sF2)))
otus_16sF3_abun <- as.data.frame(rowSums(t(otu.n0.acomp_16sF3))/ncol(t(otu.n0.acomp_16sF3)))

otus_ITSF1_abun <- as.data.frame(rowSums(t(otu.n0.acomp_ITSF1))/ncol(t(otu.n0.acomp_ITSF1)))
otus_ITSF2_abun <- as.data.frame(rowSums(t(otu.n0.acomp_ITSF2))/ncol(t(otu.n0.acomp_ITSF2)))
otus_ITSF3_abun <- as.data.frame(rowSums(t(otu.n0.acomp_ITSF3))/ncol(t(otu.n0.acomp_ITSF3)))

#Rename column
otus_16sF1_abun <-otus_16sF1_abun %>% rename(MeanAbun=`rowSums(t(otu.n0.acomp_16sF1))/ncol(t(otu.n0.acomp_16sF1))`)
otus_16sF2_abun <-otus_16sF2_abun %>% rename(MeanAbun=`rowSums(t(otu.n0.acomp_16sF2))/ncol(t(otu.n0.acomp_16sF2))`)
otus_16sF3_abun <-otus_16sF3_abun %>% rename(MeanAbun=`rowSums(t(otu.n0.acomp_16sF3))/ncol(t(otu.n0.acomp_16sF3))`)

otus_ITSF1_abun <-otus_ITSF1_abun %>% rename(MeanAbun=`rowSums(t(otu.n0.acomp_ITSF1))/ncol(t(otu.n0.acomp_ITSF1))`)
otus_ITSF2_abun <-otus_ITSF2_abun %>% rename(MeanAbun=`rowSums(t(otu.n0.acomp_ITSF2))/ncol(t(otu.n0.acomp_ITSF2))`)
otus_ITSF3_abun <-otus_ITSF3_abun %>% rename(MeanAbun=`rowSums(t(otu.n0.acomp_ITSF3))/ncol(t(otu.n0.acomp_ITSF3))`)


#Bind columns
OccAbun_16sF1<-bind_cols(otus_16sF1_occ,otus_16sF1_abun)
OccAbun_16sF2<-bind_cols(otus_16sF2_occ,otus_16sF2_abun)
OccAbun_16sF3<-bind_cols(otus_16sF3_occ,otus_16sF3_abun)

OccAbun_ITSF1<-bind_cols(otus_ITSF1_occ,otus_ITSF1_abun)
OccAbun_ITSF2<-bind_cols(otus_ITSF2_occ,otus_ITSF2_abun)
OccAbun_ITSF3<-bind_cols(otus_ITSF3_occ,otus_ITSF3_abun)


#Plot abundance-occupancy 
Plot_OccAbun_16sF1<-ggplot(OccAbun_16sF1, aes(x=MeanAbun, y=MeanOcc))+geom_point()+scale_x_log10()+
  scale_size(limits=c(0,50), name = "Relative abundance (%)")+
  ylab('Occupancy')+xlab("log10 Mean Relative Abundance")+
  theme(axis.title = element_text('black', size=17), axis.text=element_text(color='black', size=16), 
        axis.ticks=element_line(color='black'), panel.grid.major.y = element_line(color='grey50'))+
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  theme(panel.background = element_rect(fill="grey99", color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_rect(fill="white", color = "black"), 
        strip.text.y = element_text(size=16, angle=0, face = 'italic'), 
        panel.border = element_rect(color="black", fill=NA))+
  theme(legend.text=element_text(size=15,color='black'), legend.title= element_text(size=16)) +
  ggtitle("Bacteria F1")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=18))+
  scale_color_viridis_d(begin = 0.2, end = 0.8, option='inferno')+
  theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))


Plot_OccAbun_16sF2<-ggplot(OccAbun_16sF2, aes(x=MeanAbun, y=MeanOcc))+geom_point()+scale_x_log10()+
  scale_size(limits=c(0,50), name = "Relative abundance (%)")+
  ylab('Occupancy')+xlab("log10 Mean Relative Abundance")+
  theme(axis.title = element_text('black', size=17), axis.text=element_text(color='black', size=16), 
        axis.ticks=element_line(color='black'), panel.grid.major.y =element_line(color='grey50'))+
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  theme(panel.background = element_rect(fill="grey99", color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_rect(fill="white", color = "black"), 
        strip.text.y = element_text(size=16, angle=0, face = 'italic'), 
        panel.border = element_rect(color="black", fill=NA))+
  theme(legend.text=element_text(size=15,color='black'), legend.title= element_text(size=16)) +
  ggtitle("Bacteria F2")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=18))+
  scale_color_viridis_d(begin = 0.2, end = 0.8, option='inferno')+
  theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))

Plot_OccAbun_16sF3<-ggplot(OccAbun_16sF3, aes(x=MeanAbun, y=MeanOcc))+geom_point()+scale_x_log10()+
  scale_size(limits=c(0,50), name = "Relative abundance (%)")+
  ylab('Occupancy')+xlab("log10 Mean Relative Abundance")+
  theme(axis.title = element_text('black', size=17), axis.text=element_text(color='black', size=16), 
        axis.ticks=element_line(color='black'), panel.grid.major.y =element_line(color='grey50'))+
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  theme(panel.background = element_rect(fill="grey99", color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_rect(fill="white", color = "black"), 
        strip.text.y = element_text(size=16, angle=0, face = 'italic'), 
        panel.border = element_rect(color="black", fill=NA))+
  theme(legend.text=element_text(size=15,color='black'), legend.title= element_text(size=16)) +
  ggtitle("Bacteria F3")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=18))+
  scale_color_viridis_d(begin = 0.2, end = 0.8, option='inferno')+
  theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))

Plot_OccAbun_ITSF1<-ggplot(OccAbun_ITSF1, aes(x=MeanAbun, y=MeanOcc))+geom_point()+scale_x_log10()+
  scale_size(limits=c(0,50), name = "Relative abundance (%)")+
  ylab('Occupancy')+xlab("log10 Mean Relative Abundance")+
  theme(axis.title = element_text('black', size=17), axis.text=element_text(color='black', size=16), 
        axis.ticks=element_line(color='black'), panel.grid.major.y =element_line(color='grey50'))+
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  theme(panel.background = element_rect(fill="grey99", color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_rect(fill="white", color = "black"), 
        strip.text.y = element_text(size=16, angle=0, face = 'italic'), 
        panel.border = element_rect(color="black", fill=NA))+
  theme(legend.text=element_text(size=15,color='black'), legend.title= element_text(size=16)) +
  ggtitle("Fungi F1")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=18))+
  scale_color_viridis_d(begin = 0.2, end = 0.8, option='inferno')+
  theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))


Plot_OccAbun_ITSF2<-ggplot(OccAbun_ITSF2, aes(x=MeanAbun, y=MeanOcc))+geom_point()+scale_x_log10()+
  scale_size(limits=c(0,50), name = "Relative abundance (%)")+
  ylab('Occupancy')+xlab("log10 Mean Relative Abundance")+
  theme(axis.title = element_text('black', size=17), axis.text=element_text(color='black', size=16), 
        axis.ticks=element_line(color='black'), panel.grid.major.y =element_line(color='grey50'))+
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  theme(panel.background = element_rect(fill="grey99", color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_rect(fill="white", color = "black"), 
        strip.text.y = element_text(size=16, angle=0, face = 'italic'), 
        panel.border = element_rect(color="black", fill=NA))+
  theme(legend.text=element_text(size=15,color='black'), legend.title= element_text(size=16)) +
  ggtitle("Fungi F2")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=18))+
  scale_color_viridis_d(begin = 0.2, end = 0.8, option='inferno')+
  theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))

Plot_OccAbun_ITSF3<-ggplot(OccAbun_ITSF3, aes(x=MeanAbun, y=MeanOcc))+geom_point()+scale_x_log10()+
  scale_size(limits=c(0,50), name = "Relative abundance (%)")+
  ylab('Occupancy')+xlab("log10 Mean Relative Abundance")+
  theme(axis.title = element_text('black', size=17), axis.text=element_text(color='black', size=16), 
        axis.ticks=element_line(color='black'), panel.grid.major.y =element_line(color='grey50'))+
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  theme(panel.background = element_rect(fill="grey99", color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_rect(fill="white", color = "black"), 
        strip.text.y = element_text(size=16, angle=0, face = 'italic'), 
        panel.border = element_rect(color="black", fill=NA))+
  theme(legend.text=element_text(size=15,color='black'), legend.title= element_text(size=16)) +
  ggtitle("Fungi F3")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=18))+
  scale_color_viridis_d(begin = 0.2, end = 0.8, option='inferno')+
  theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))


#combine plots
OccAbun_plots = plot_grid(Plot_OccAbun_16sF1,Plot_OccAbun_16sF2, Plot_OccAbun_16sF3,
                          Plot_OccAbun_ITSF1, Plot_OccAbun_ITSF2, Plot_OccAbun_ITSF3,
                               ncol=3, nrow=2, labels = c("A","B", "C", "D", "E", "F"), label_size = 20)
OccAbun_plots


#Subset OTUs with Occupancy=1
otus_16sF1_occ1<-subset(OccAbun_16sF1, MeanOcc==1)%>%
  arrange(rownames(otus_16sF1_occ1))
otus_16sF2_occ1<-subset(OccAbun_16sF2, MeanOcc==1)%>%
  arrange(rownames(otus_16sF2_occ1))
otus_16sF3_occ1<-subset(OccAbun_16sF3, MeanOcc==1)%>%
  arrange(rownames(otus_16sF3_occ1))

otus_ITSF1_occ1<-subset(OccAbun_ITSF1, MeanOcc==1)%>%
  arrange(rownames(otus_ITSF1_occ1))
otus_ITSF2_occ1<-subset(OccAbun_ITSF2, MeanOcc==1)%>%
  arrange(rownames(otus_ITSF2_occ1))
otus_ITSF3_occ1<-subset(OccAbun_ITSF3, MeanOcc==1)%>%
  arrange(rownames(otus_ITSF3_occ1))

#Add taxonomic information
otus_16sF1_occ1_tax<-subset(taxon_16s, rownames(taxon_16s) %in% rownames(otus_16sF1_occ1))
otus_16sF2_occ1_tax<-subset(taxon_16s, rownames(taxon_16s) %in% rownames(otus_16sF2_occ1))
otus_16sF3_occ1_tax<-subset(taxon_16s, rownames(taxon_16s) %in% rownames(otus_16sF3_occ1))

otus_ITSF1_occ1_tax<-subset(taxon_ITS, rownames(taxon_ITS) %in% rownames(otus_ITSF1_occ1))
otus_ITSF2_occ1_tax<-subset(taxon_ITS, rownames(taxon_ITS) %in% rownames(otus_ITSF2_occ1))
otus_ITSF3_occ1_tax<-subset(taxon_ITS, rownames(taxon_ITS) %in% rownames(otus_ITSF3_occ1))

#Bind columns
otus_16sF1_occ1<-bind_cols(otus_16sF1_occ1,otus_16sF1_occ1_tax)
otus_16sF2_occ1<-bind_cols(otus_16sF2_occ1,otus_16sF2_occ1_tax)
otus_16sF3_occ1<-bind_cols(otus_16sF3_occ1,otus_16sF3_occ1_tax)

otus_ITSF1_occ1<-bind_cols(otus_ITSF1_occ1,otus_ITSF1_occ1_tax)
otus_ITSF2_occ1<-bind_cols(otus_ITSF2_occ1,otus_ITSF2_occ1_tax)
otus_ITSF3_occ1<-bind_cols(otus_ITSF3_occ1,otus_ITSF3_occ1_tax)




#### Occupancy - abundance curves by collection time point ####
#Calculate mean occupancy by collection time. 
#We consider that the the samples collected at each facility each time are technical replicates.
#We want to identify those OTUs that were present in each facility at every timepoint throughout two years
phyloseq_16sF1<-phyloseq(otu_table(otus_16sF1, taxa_are_rows = TRUE), TAX_16s, sample_data(metadata_16sF1))
phyloseq_16sF2<-phyloseq(otu_table(otus_16sF2, taxa_are_rows = TRUE), TAX_16s, sample_data(metadata_16sF2))
phyloseq_16sF3<-phyloseq(otu_table(otus_16sF3, taxa_are_rows = TRUE), TAX_16s, sample_data(metadata_16sF3))

phyloseq_ITSF1<-phyloseq(otu_table(otus_ITSF1, taxa_are_rows = TRUE), TAX_ITS, sample_data(metadata_ITSF1))
phyloseq_ITSF2<-phyloseq(otu_table(otus_ITSF2, taxa_are_rows = TRUE), TAX_ITS, sample_data(metadata_ITSF2))
phyloseq_ITSF3<-phyloseq(otu_table(otus_ITSF3, taxa_are_rows = TRUE), TAX_ITS, sample_data(metadata_ITSF3))

#Merge samples taken in the same collection date
phyloseq_16sF1_time<- merge_samples(phyloseq_16sF1, group = "Date")
phyloseq_16sF2_time<- merge_samples(phyloseq_16sF2, group = "Date")
phyloseq_16sF3_time<- merge_samples(phyloseq_16sF3, group = "Date")

phyloseq_ITSF1_time<- merge_samples(phyloseq_ITSF1, group = "Date")
phyloseq_ITSF2_time<- merge_samples(phyloseq_ITSF2, group = "Date")
phyloseq_ITSF3_time<- merge_samples(phyloseq_ITSF3, group = "Date")

#Export OTU table for each facility from Phyloseq 
otus_16sF1_time<-as.data.frame(as(otu_table(phyloseq_16sF1_time), "matrix"))
otus_16sF2_time<-as.data.frame(as(otu_table(phyloseq_16sF2_time), "matrix"))
otus_16sF3_time<-as.data.frame(as(otu_table(phyloseq_16sF3_time), "matrix"))

otus_ITSF1_time<-as.data.frame(as(otu_table(phyloseq_ITSF1_time), "matrix"))
otus_ITSF2_time<-as.data.frame(as(otu_table(phyloseq_ITSF2_time), "matrix"))
otus_ITSF3_time<-as.data.frame(as(otu_table(phyloseq_ITSF3_time), "matrix"))

#Transpose otu table to have OTUs in rows
otus_16sF1_time<-t(otus_16sF1_time)
otus_16sF2_time<-t(otus_16sF2_time)
otus_16sF3_time<-t(otus_16sF3_time)

otus_ITSF1_time<-t(otus_ITSF1_time)
otus_ITSF2_time<-t(otus_ITSF2_time)
otus_ITSF3_time<-t(otus_ITSF3_time)

#Calculate presence absence
otus_16sF1_time_PA <- 1*((otus_16sF1_time>0)==1) 
otus_16sF2_time_PA <- 1*((otus_16sF2_time>0)==1) 
otus_16sF3_time_PA <- 1*((otus_16sF3_time>0)==1) 

otus_ITSF1_time_PA <- 1*((otus_ITSF1_time>0)==1) 
otus_ITSF2_time_PA <- 1*((otus_ITSF2_time>0)==1) 
otus_ITSF3_time_PA <- 1*((otus_ITSF3_time>0)==1) 

#Calculate mean occupancy for each OTU
otus_16sF1_time_occ <- as.data.frame(rowSums(otus_16sF1_time_PA)/ncol(otus_16sF1_time_PA))
otus_16sF2_time_occ <- as.data.frame(rowSums(otus_16sF2_time_PA)/ncol(otus_16sF2_time_PA))
otus_16sF3_time_occ <- as.data.frame(rowSums(otus_16sF3_time_PA)/ncol(otus_16sF3_time_PA))

otus_ITSF1_time_occ <- as.data.frame(rowSums(otus_ITSF1_time_PA)/ncol(otus_ITSF1_time_PA))
otus_ITSF2_time_occ <- as.data.frame(rowSums(otus_ITSF2_time_PA)/ncol(otus_ITSF2_time_PA))
otus_ITSF3_time_occ <- as.data.frame(rowSums(otus_ITSF3_time_PA)/ncol(otus_ITSF3_time_PA))

#Rename column
otus_16sF1_time_occ <-otus_16sF1_time_occ %>% rename(MeanOcc=`rowSums(otus_16sF1_time_PA)/ncol(otus_16sF1_time_PA)`)
otus_16sF2_time_occ <-otus_16sF2_time_occ %>% rename(MeanOcc=`rowSums(otus_16sF2_time_PA)/ncol(otus_16sF2_time_PA)`)
otus_16sF3_time_occ <-otus_16sF3_time_occ %>% rename(MeanOcc=`rowSums(otus_16sF3_time_PA)/ncol(otus_16sF3_time_PA)`)

otus_ITSF1_time_occ <-otus_ITSF1_time_occ %>% rename(MeanOcc=`rowSums(otus_ITSF1_time_PA)/ncol(otus_ITSF1_time_PA)`)
otus_ITSF2_time_occ <-otus_ITSF2_time_occ %>% rename(MeanOcc=`rowSums(otus_ITSF2_time_PA)/ncol(otus_ITSF2_time_PA)`)
otus_ITSF3_time_occ <-otus_ITSF3_time_occ %>% rename(MeanOcc=`rowSums(otus_ITSF3_time_PA)/ncol(otus_ITSF3_time_PA)`)


#Calculate relative abundance for each OTU by Facility
#First replace zero values with a small non- zero value before clr transformation. 
#Uses Count Zero Multiplicative method to replace zeros and outputs pseudo-counts
otu.n0_time_16sF1<-t(cmultRepl(t(otus_16sF1_time), label=0, method="CZM", output="p-counts")) 
otu.n0_time_16sF2<-t(cmultRepl(t(otus_16sF2_time), label=0, method="CZM", output="p-counts")) 
otu.n0_time_16sF3<-t(cmultRepl(t(otus_16sF3_time), label=0, method="CZM", output="p-counts")) 

otu.n0_time_ITSF1<-t(cmultRepl(t(otus_ITSF1_time), label=0, method="CZM", output="p-counts")) 
otu.n0_time_ITSF2<-t(cmultRepl(t(otus_ITSF2_time), label=0, method="CZM", output="p-counts")) 
otu.n0_time_ITSF3<-t(cmultRepl(t(otus_ITSF3_time), label=0, method="CZM", output="p-counts")) 


#Then convert to compositions with Aitchinson 
otu.n0.acomp_16sF1_time<-as.data.frame(acomp(t(otu.n0_time_16sF1)), total=1)
otu.n0.acomp_16sF2_time<-as.data.frame(acomp(t(otu.n0_time_16sF2)), total=1)
otu.n0.acomp_16sF3_time<-as.data.frame(acomp(t(otu.n0_time_16sF3)), total=1)

otu.n0.acomp_ITSF1_time<-as.data.frame(acomp(t(otu.n0_time_ITSF1)), total=1)
otu.n0.acomp_ITSF2_time<-as.data.frame(acomp(t(otu.n0_time_ITSF2)), total=1)
otu.n0.acomp_ITSF3_time<-as.data.frame(acomp(t(otu.n0_time_ITSF3)), total=1)

#Calculate the mean relative abundance for each OTU in each facility
otus_16sF1_abun_time <- as.data.frame(rowSums(t(otu.n0.acomp_16sF1_time))/ncol(t(otu.n0.acomp_16sF1_time)))
otus_16sF2_abun_time <- as.data.frame(rowSums(t(otu.n0.acomp_16sF2_time))/ncol(t(otu.n0.acomp_16sF2_time)))
otus_16sF3_abun_time <- as.data.frame(rowSums(t(otu.n0.acomp_16sF3_time))/ncol(t(otu.n0.acomp_16sF3_time)))

otus_ITSF1_abun_time <- as.data.frame(rowSums(t(otu.n0.acomp_ITSF1_time))/ncol(t(otu.n0.acomp_ITSF1_time)))
otus_ITSF2_abun_time <- as.data.frame(rowSums(t(otu.n0.acomp_ITSF2_time))/ncol(t(otu.n0.acomp_ITSF2_time)))
otus_ITSF3_abun_time <- as.data.frame(rowSums(t(otu.n0.acomp_ITSF3_time))/ncol(t(otu.n0.acomp_ITSF3_time)))

#Rename column
otus_16sF1_abun_time <-otus_16sF1_abun_time %>% rename(MeanAbun=`rowSums(t(otu.n0.acomp_16sF1_time))/ncol(t(otu.n0.acomp_16sF1_time))`)
otus_16sF2_abun_time <-otus_16sF2_abun_time %>% rename(MeanAbun=`rowSums(t(otu.n0.acomp_16sF2_time))/ncol(t(otu.n0.acomp_16sF2_time))`)
otus_16sF3_abun_time <-otus_16sF3_abun_time %>% rename(MeanAbun=`rowSums(t(otu.n0.acomp_16sF3_time))/ncol(t(otu.n0.acomp_16sF3_time))`)

otus_ITSF1_abun_time <-otus_ITSF1_abun_time %>% rename(MeanAbun=`rowSums(t(otu.n0.acomp_ITSF1_time))/ncol(t(otu.n0.acomp_ITSF1_time))`)
otus_ITSF2_abun_time <-otus_ITSF2_abun_time %>% rename(MeanAbun=`rowSums(t(otu.n0.acomp_ITSF2_time))/ncol(t(otu.n0.acomp_ITSF2_time))`)
otus_ITSF3_abun_time <-otus_ITSF3_abun_time %>% rename(MeanAbun=`rowSums(t(otu.n0.acomp_ITSF3_time))/ncol(t(otu.n0.acomp_ITSF3_time))`)

#Bind columns
OccAbun_16sF1_time<-bind_cols(otus_16sF1_time_occ,otus_16sF1_abun_time)
OccAbun_16sF2_time<-bind_cols(otus_16sF2_time_occ,otus_16sF2_abun_time)
OccAbun_16sF3_time<-bind_cols(otus_16sF3_time_occ,otus_16sF3_abun_time)

OccAbun_ITSF1_time<-bind_cols(otus_ITSF1_time_occ,otus_ITSF1_abun_time)
OccAbun_ITSF2_time<-bind_cols(otus_ITSF2_time_occ,otus_ITSF2_abun_time)
OccAbun_ITSF3_time<-bind_cols(otus_ITSF3_time_occ,otus_ITSF3_abun_time)

#Add color column for Otus with Occurrance=1
OccAbun_16sF1_time$color<-ifelse(OccAbun_16sF1_time$MeanOcc==1, "yes", "no")
OccAbun_16sF2_time$color<-ifelse(OccAbun_16sF2_time$MeanOcc==1, "yes", "no")
OccAbun_16sF3_time$color<-ifelse(OccAbun_16sF3_time$MeanOcc==1, "yes", "no")

OccAbun_ITSF1_time$color<-ifelse(OccAbun_ITSF1_time$MeanOcc==1, "yes", "no")
OccAbun_ITSF2_time$color<-ifelse(OccAbun_ITSF2_time$MeanOcc==1, "yes", "no")
OccAbun_ITSF3_time$color<-ifelse(OccAbun_ITSF3_time$MeanOcc==1, "yes", "no")

#Plot abundance-occupancy 
#Fig 4 A,B,C,E,F,G
Plot_OccAbun_16sF1_time<-ggplot(OccAbun_16sF1_time, aes(x=MeanAbun, y=MeanOcc, color=color))+geom_point()+scale_x_log10()+
  scale_size(limits=c(0,50))+
  ylab('Occupancy')+xlab("log10 Mean Relative Abundance")+
  theme(axis.title = element_text(color='black', size=17), axis.text=element_text(color='black', size=16), 
        axis.ticks=element_line(color='black'), panel.grid.major.y = element_line(color='grey50'))+
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  theme(panel.background = element_rect(fill="white", color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(panel.border = element_rect(color="black", fill=NA))+
  theme(legend.position = 'none') +
  ggtitle("Bacteria F1")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=18))+
  scale_color_manual(values=c("#00000080","#420A6880"))+
  theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))

Plot_OccAbun_16sF2_time<-ggplot(OccAbun_16sF2_time, aes(x=MeanAbun, y=MeanOcc, color=color))+geom_point()+scale_x_log10()+
  scale_size(limits=c(0,50))+
  ylab('Occupancy')+xlab("log10 Mean Relative Abundance")+
  theme(axis.title = element_text(color='black', size=17), axis.text=element_text(color='black', size=16), 
        axis.ticks=element_line(color='black'), panel.grid.major.y = element_line(color='grey50'))+
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  theme(panel.background = element_rect(fill="white", color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(panel.border = element_rect(color="black", fill=NA))+
  theme(legend.position = 'none') +
  ggtitle("Bacteria F2")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=18))+
  scale_color_manual(values=c("#00000080","#BB375480"))+
  theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))

Plot_OccAbun_16sF3_time<-ggplot(OccAbun_16sF3_time, aes(x=MeanAbun, y=MeanOcc, color=color))+geom_point()+scale_x_log10()+
  scale_size(limits=c(0,50))+
  ylab('Occupancy')+xlab("log10 Mean Relative Abundance")+
  theme(axis.title = element_text(color='black', size=17), axis.text=element_text(color='black', size=16), 
        axis.ticks=element_line(color='black'), panel.grid.major.y = element_line(color='grey50'))+
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  theme(panel.background = element_rect(fill="white", color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(panel.border = element_rect(color="black", fill=NA))+
  theme(legend.position = 'none') +
  ggtitle("Bacteria F3")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=18))+
  scale_color_manual(values=c("#00000080","#FCA50A80"))+
  theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))

Plot_OccAbun_ITSF1_time<-ggplot(OccAbun_ITSF1_time, aes(x=MeanAbun, y=MeanOcc, color=color))+geom_point()+scale_x_log10()+
  scale_size(limits=c(0,50))+
  ylab('Occupancy')+xlab("log10 Mean Relative Abundance")+
  theme(axis.title = element_text(color='black', size=17), axis.text=element_text(color='black', size=16), 
        axis.ticks=element_line(color='black'), panel.grid.major.y = element_line(color='grey50'))+
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  theme(panel.background = element_rect(fill="white", color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(panel.border = element_rect(color="black", fill=NA))+
  theme(legend.position = 'none') +
  ggtitle("Fungi F1")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=18))+
  scale_color_manual(values=c("#00000080","#420A6880"))+
  theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))

Plot_OccAbun_ITSF2_time<-ggplot(OccAbun_ITSF2_time, aes(x=MeanAbun, y=MeanOcc, color=color))+geom_point()+scale_x_log10()+
  scale_size(limits=c(0,50))+
  ylab('Occupancy')+xlab("log10 Mean Relative Abundance")+
  theme(axis.title = element_text(color='black', size=17), axis.text=element_text(color='black', size=16), 
        axis.ticks=element_line(color='black'), panel.grid.major.y = element_line(color='grey50'))+
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  theme(panel.background = element_rect(fill="white", color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(panel.border = element_rect(color="black", fill=NA))+
  theme(legend.position = 'none') +
  ggtitle("Fungi F2")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=18))+
  scale_color_manual(values=c("#00000080","#BB375480"))+
  theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))

Plot_OccAbun_ITSF3_time<-ggplot(OccAbun_ITSF3_time, aes(x=MeanAbun, y=MeanOcc, color=color))+geom_point()+scale_x_log10()+
  scale_size(limits=c(0,50))+
  ylab('Occupancy')+xlab("log10 Mean Relative Abundance")+
  theme(axis.title = element_text(color='black', size=17), axis.text=element_text(color='black', size=16), 
        axis.ticks=element_line(color='black'), panel.grid.major.y = element_line(color='grey50'))+
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  theme(panel.background = element_rect(fill="white", color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(panel.border = element_rect(color="black", fill=NA))+
  theme(legend.position = 'none') +
  ggtitle("Fungi F3")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=18))+
  scale_color_manual(values=c("#00000080","#FCA50A80"))+
  theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))



#combine plots
OccAbun_plots_time = plot_grid(Plot_OccAbun_16sF1_time,Plot_OccAbun_16sF2_time, Plot_OccAbun_16sF3_time,
                          Plot_OccAbun_ITSF1_time, Plot_OccAbun_ITSF2_time, Plot_OccAbun_ITSF3_time,
                          ncol=3, nrow=2, labels = c("A","B", "C", "D", "E", "F"), label_size = 20)
OccAbun_plots_time

ggsave("Abundance Occurance by time.png", plot=OccAbun_plots_time, device="png", width=10, height=8, units="in",dpi=600)
ggsave("Abundance Occurance by time.svg", plot=OccAbun_plots_time, device="svg", width=10, height=8, units="in",dpi=600)


#Subset OTUs with Occupancy=1
otus_16sF1_occ1_time<-subset(OccAbun_16sF1_time, MeanOcc==1)
otus_16sF2_occ1_time<-subset(OccAbun_16sF2_time, MeanOcc==1)
otus_16sF3_occ1_time<-subset(OccAbun_16sF3_time, MeanOcc==1)

otus_ITSF1_occ1_time<-subset(OccAbun_ITSF1_time, MeanOcc==1)
otus_ITSF2_occ1_time<-subset(OccAbun_ITSF2_time, MeanOcc==1)
otus_ITSF3_occ1_time<-subset(OccAbun_ITSF3_time, MeanOcc==1)

#Reorder by OTU name
otus_16sF1_occ1_time<-otus_16sF1_occ1_time[order(rownames(otus_16sF1_occ1_time)),]
otus_16sF2_occ1_time<-otus_16sF2_occ1_time[order(rownames(otus_16sF2_occ1_time)),]
otus_16sF3_occ1_time<-otus_16sF3_occ1_time[order(rownames(otus_16sF3_occ1_time)),]

otus_ITSF1_occ1_time<-otus_ITSF1_occ1_time[order(rownames(otus_ITSF1_occ1_time)),]
otus_ITSF2_occ1_time<-otus_ITSF2_occ1_time[order(rownames(otus_ITSF2_occ1_time)),]
otus_ITSF3_occ1_time<-otus_ITSF3_occ1_time[order(rownames(otus_ITSF3_occ1_time)),]

#Get taxonomic information for core taxa
otus_16sF1_occ1_time_tax<-subset(taxon_16s, rownames(taxon_16s) %in% rownames(otus_16sF1_occ1_time))
otus_16sF2_occ1_time_tax<-subset(taxon_16s, rownames(taxon_16s) %in% rownames(otus_16sF2_occ1_time))
otus_16sF3_occ1_time_tax<-subset(taxon_16s, rownames(taxon_16s) %in% rownames(otus_16sF3_occ1_time))

otus_ITSF1_occ1_time_tax<-subset(taxon_ITS, rownames(taxon_ITS) %in% rownames(otus_ITSF1_occ1_time))
otus_ITSF2_occ1_time_tax<-subset(taxon_ITS, rownames(taxon_ITS) %in% rownames(otus_ITSF2_occ1_time))
otus_ITSF3_occ1_time_tax<-subset(taxon_ITS, rownames(taxon_ITS) %in% rownames(otus_ITSF3_occ1_time))

#Bind columns
otus_16sF1_occ1_time<-bind_cols(otus_16sF1_occ1_time,otus_16sF1_occ1_time_tax)
otus_16sF2_occ1_time<-bind_cols(otus_16sF2_occ1_time,otus_16sF2_occ1_time_tax)
otus_16sF3_occ1_time<-bind_cols(otus_16sF3_occ1_time,otus_16sF3_occ1_time_tax)

otus_ITSF1_occ1_time<-bind_cols(otus_ITSF1_occ1_time,otus_ITSF1_occ1_time_tax)
otus_ITSF2_occ1_time<-bind_cols(otus_ITSF2_occ1_time,otus_ITSF2_occ1_time_tax)
otus_ITSF3_occ1_time<-bind_cols(otus_ITSF3_occ1_time,otus_ITSF3_occ1_time_tax)


#Make list of core OTUs per Facility
F1_16s_core<-rownames(otus_16sF1_occ1_time)
F2_16s_core<-rownames(otus_16sF2_occ1_time)
F3_16s_core<-rownames(otus_16sF3_occ1_time)

F1_ITS_core<-rownames(otus_ITSF1_occ1_time)
F2_ITS_core<-rownames(otus_ITSF2_occ1_time)
F3_ITS_core<-rownames(otus_ITSF3_occ1_time)


#Combine list of core OTUs for all facilities
coretaxa_16s_list<-list(F1=F1_16s_core, F2=F2_16s_core, F3=F3_16s_core)
coretaxa_ITS_list<-list(F1=F1_ITS_core, F2=F2_ITS_core, F3=F3_ITS_core)

#Make dataframe for all core taxa
core_16s<-data.frame(OTU=c(F1_16s_core, F2_16s_core, F3_16s_core), Facility=c(rep("F1",89),rep("F2",90),rep("F3",126)))
core_ITS<-data.frame(OTU=c(F1_ITS_core, F2_ITS_core, F3_ITS_core), Facility=c(rep("F1",78),rep("F2",41),rep("F3",67)))

#Get taxa information for core taxa 
core_16s_taxa<-subset(taxon_16s, rownames(taxon_16s) %in% core_16s$OTU)
core_ITS_taxa<-subset(taxon_ITS, rownames(taxon_ITS) %in% core_ITS$OTU)

#Write csv file
write.csv(core_16s_taxa, file="Core 16s.csv")
write.csv(core_ITS_taxa, file="Core ITS.csv")

#Make Venn diagrams
#Fig 4D, H
library(ggvenn)

venn_16s<-ggvenn(coretaxa_16s_list, show_percentage = FALSE,
                   fill_color = c('#420A68','#BB3754','#FCA50A'), set_name_size = 4, stroke_size = 0.5)+
  ggtitle("Bacteria")+theme(plot.title = element_text(hjust = 0.5))

venn_ITS<-ggvenn(coretaxa_ITS_list, show_percentage = FALSE,
                 fill_color = c('#420A68','#BB3754','#FCA50A'), set_name_size = 4, stroke_size = 0.5)+
  ggtitle("Fungi")+theme(plot.title = element_text(hjust = 0.5))


#Combine Venn diagrams
Venn_core = plot_grid(venn_16s, venn_ITS, 
                         ncol=2, nrow=1, labels = c("G","H"), label_size = 20, vjust = 2, hjust = -1.5)
Venn_core

ggsave("Venn_core.png", plot=Venn_core, device="png", width=8, height=4, units="in", dpi=600)
ggsave("Venn_core.svg", plot=Venn_core, device="svg", width=8, height=4, units="in", dpi=600)


#Print taxa that are shared across facilities
shared_core_16s<-intersect(intersect(F1_16s_core, F2_16s_core),F3_16s_core)
shared_core_ITS<-intersect(intersect(F1_ITS_core, F2_ITS_core),F3_ITS_core)

#Get taxa information for core taxa shared by the three facilities
shared_core_16s_taxa<-subset(taxon_16s, rownames(taxon_16s) %in% shared_core_16s)
shared_core_ITS_taxa<-subset(taxon_ITS, rownames(taxon_ITS) %in% shared_core_ITS)

#Write csv file
write.csv(shared_core_16s_taxa, file="Shared core 16s.csv")
write.csv(shared_core_ITS_taxa, file="Shared core ITS.csv")


#Transpose abundance table
otu.n0.acomp_16sF1_time<-t(otu.n0.acomp_16sF1_time)
otu.n0.acomp_16sF2_time<-t(otu.n0.acomp_16sF2_time)
otu.n0.acomp_16sF3_time<-t(otu.n0.acomp_16sF3_time)

otu.n0.acomp_ITSF1_time<-t(otu.n0.acomp_ITSF1_time)
otu.n0.acomp_ITSF2_time<-t(otu.n0.acomp_ITSF2_time)
otu.n0.acomp_ITSF3_time<-t(otu.n0.acomp_ITSF3_time)

#Reorder by rowname
otu.n0.acomp_16sF1_time <-otu.n0.acomp_16sF1_time[order(rownames(otu.n0.acomp_16sF1_time)),]
otu.n0.acomp_16sF2_time <-otu.n0.acomp_16sF2_time[order(rownames(otu.n0.acomp_16sF2_time)),]
otu.n0.acomp_16sF3_time <-otu.n0.acomp_16sF3_time[order(rownames(otu.n0.acomp_16sF3_time)),]

otu.n0.acomp_ITSF1_time <-otu.n0.acomp_ITSF1_time[order(rownames(otu.n0.acomp_ITSF1_time)),]
otu.n0.acomp_ITSF2_time <-otu.n0.acomp_ITSF2_time[order(rownames(otu.n0.acomp_ITSF2_time)),]
otu.n0.acomp_ITSF3_time <-otu.n0.acomp_ITSF3_time[order(rownames(otu.n0.acomp_ITSF3_time)),]

#subset abundance table to include only taxa with occupancy 1 in all facilities
abund_16sF1_occ1_time<-as.data.frame(subset(otu.n0.acomp_16sF1_time, rownames(otu.n0.acomp_16sF1_time) %in% shared_core_16s))
abund_16sF2_occ1_time<-as.data.frame(subset(otu.n0.acomp_16sF2_time, rownames(otu.n0.acomp_16sF2_time) %in% shared_core_16s))
abund_16sF3_occ1_time<-as.data.frame(subset(otu.n0.acomp_16sF3_time, rownames(otu.n0.acomp_16sF3_time) %in% shared_core_16s))

abund_ITSF1_occ1_time<-as.data.frame(subset(otu.n0.acomp_ITSF1_time, rownames(otu.n0.acomp_ITSF1_time) %in% shared_core_ITS))
abund_ITSF2_occ1_time<-as.data.frame(subset(otu.n0.acomp_ITSF2_time, rownames(otu.n0.acomp_ITSF2_time) %in% shared_core_ITS))
abund_ITSF3_occ1_time<-as.data.frame(subset(otu.n0.acomp_ITSF3_time, rownames(otu.n0.acomp_ITSF3_time) %in% shared_core_ITS))

#Add OTU name column
abund_16sF1_occ1_time$OTU<-rownames(abund_16sF1_occ1_time)
abund_16sF2_occ1_time$OTU<-rownames(abund_16sF2_occ1_time)
abund_16sF3_occ1_time$OTU<-rownames(abund_16sF3_occ1_time)

abund_ITSF1_occ1_time$OTU<-rownames(abund_ITSF1_occ1_time)
abund_ITSF2_occ1_time$OTU<-rownames(abund_ITSF2_occ1_time)
abund_ITSF3_occ1_time$OTU<-rownames(abund_ITSF3_occ1_time)

#Add taxon family column
abund_16sF1_occ1_time$Family<-subset(taxon_16s$Family, rownames(taxon_16s) %in% shared_core_16s)
abund_16sF2_occ1_time$Family<-subset(taxon_16s$Family, rownames(taxon_16s) %in% shared_core_16s)
abund_16sF3_occ1_time$Family<-subset(taxon_16s$Family, rownames(taxon_16s) %in% shared_core_16s)

abund_ITSF1_occ1_time$Family<-subset(taxon_ITS$Family, rownames(taxon_ITS) %in% shared_core_ITS)
abund_ITSF2_occ1_time$Family<-subset(taxon_ITS$Family, rownames(taxon_ITS) %in% shared_core_ITS)
abund_ITSF3_occ1_time$Family<-subset(taxon_ITS$Family, rownames(taxon_ITS) %in% shared_core_ITS)

#Convert to long format
library(tidyr)
abund_16sF1_occ1_time_long<-gather(abund_16sF1_occ1_time, Date, Abundance, `1/11/2019`:`4/4/2019`, factor_key=TRUE)
abund_16sF2_occ1_time_long<-gather(abund_16sF2_occ1_time, Date, Abundance, `1/11/2019`:`4/4/2019`, factor_key=TRUE)
abund_16sF3_occ1_time_long<-gather(abund_16sF3_occ1_time, Date, Abundance, `1/11/2019`:`4/4/2019`, factor_key=TRUE)

abund_ITSF1_occ1_time_long<-gather(abund_ITSF1_occ1_time, Date, Abundance, `1/11/2019`:`4/4/2019`, factor_key=TRUE)
abund_ITSF2_occ1_time_long<-gather(abund_ITSF2_occ1_time, Date, Abundance, `1/11/2019`:`4/4/2019`, factor_key=TRUE)
abund_ITSF3_occ1_time_long<-gather(abund_ITSF3_occ1_time, Date, Abundance, `1/11/2019`:`4/4/2019`, factor_key=TRUE)

#Add Facility id
abund_16sF1_occ1_time_long$Facility<-rep("F1", 950)
abund_16sF2_occ1_time_long$Facility<-rep("F2", 950)
abund_16sF3_occ1_time_long$Facility<-rep("F3", 950)

abund_ITSF1_occ1_time_long$Facility<-rep("F1", 700)
abund_ITSF2_occ1_time_long$Facility<-rep("F2", 700)
abund_ITSF3_occ1_time_long$Facility<-rep("F3", 700)

#Recalculate abundance as percentage
abund_16sF1_occ1_time_long$RA<-abund_16sF1_occ1_time_long$Abundance*100
abund_16sF2_occ1_time_long$RA<-abund_16sF2_occ1_time_long$Abundance*100
abund_16sF3_occ1_time_long$RA<-abund_16sF3_occ1_time_long$Abundance*100

abund_ITSF1_occ1_time_long$RA<-abund_ITSF1_occ1_time_long$Abundance*100
abund_ITSF2_occ1_time_long$RA<-abund_ITSF2_occ1_time_long$Abundance*100
abund_ITSF3_occ1_time_long$RA<-abund_ITSF3_occ1_time_long$Abundance*100

#Bind rows for core bacteria and fungi
abund_16s_occ1 <-bind_rows(abund_16sF1_occ1_time_long, abund_16sF2_occ1_time_long, abund_16sF3_occ1_time_long)
abund_ITS_occ1 <-bind_rows(abund_ITSF1_occ1_time_long, abund_ITSF2_occ1_time_long, abund_ITSF3_occ1_time_long)


#Remove tax indices from ITS
abund_ITS_occ1$Family_clean<-gsub("f__","",abund_ITS_occ1_mean$Family)
abund_ITS_occ1$Family_clean<-gsub("o__","",abund_ITS_occ1_mean$Family_clean)
abund_ITS_occ1$Family_clean<-gsub("p__","",abund_ITS_occ1_mean$Family_clean)
abund_ITS_occ1$Family_clean<-gsub("c__","",abund_ITS_occ1_mean$Family_clean)

#Plots of shared core by time
shared_core_16s<-ggplot(abund_16s_occ1, aes(x=Date, y=RA, group=OTU, color=OTU))+geom_point()+geom_line()+
  facet_grid(Family~Facility, scales="free_y")+
  theme(axis.text = element_text(color='black'), axis.ticks = element_line(color='black'))+
  theme(legend.position = "rigth")+ 
  theme(panel.background = element_rect(fill='grey99', color=NA), plot.background = element_rect(fill='transparent', color=NA))+
  theme(strip.background = element_rect(fill="white", color='black'), strip.text.y = element_text(color='black', angle=0, face='italic'),
        panel.border=element_rect(color='black', fill=NA))+
  theme(axis.text.x = element_text(angle=90))+
  ylab("Relative abundance (%)")+xlab("Sample collection date")

ggsave("Abundance_core_bacteria.png", plot=shared_core_16s, device="png", width=18, height=20, units="in", dpi=600)
ggsave("Abundance_core_bacteria.svg", plot=shared_core_16s, device="svg", width=18, height=20, units="in", dpi=600)

shared_core_ITS<-ggplot(abund_ITS_occ1, aes(x=Date, y=RA, group=OTU, color=OTU))+geom_point()+geom_line()+
  facet_grid(Family_clean~Facility, scales="free_y")+
  theme(axis.text = element_text(color='black'), axis.ticks = element_line(color='black'))+
  theme(legend.position = "rigth")+ 
  theme(panel.background = element_rect(fill='grey99', color=NA), plot.background = element_rect(fill='transparent', color=NA))+
  theme(strip.background = element_rect(fill="white", color='black'), strip.text.y = element_text(color='black', angle=0, face='italic'),
        panel.border=element_rect(color='black', fill=NA))+
  theme(axis.text.x = element_text(angle=90))+
  ylab("Relative abundance (%)")+xlab("Sample collection date")

ggsave("Abundance_core_fungi.png", plot=shared_core_ITS, device="png", width=18, height=20, units="in", dpi=600)
ggsave("Abundance_core_fungi.svg", plot=shared_core_ITS, device="svg", width=18, height=20, units="in", dpi=600)




#### Network analyses ####
#Subset OTUs with Occupancy>0.5
otus_16sF1_occ0.5_time<-subset(OccAbun_16sF1_time, MeanOcc>=0.5)
otus_16sF2_occ0.5_time<-subset(OccAbun_16sF2_time, MeanOcc>=0.5)
otus_16sF3_occ0.5_time<-subset(OccAbun_16sF3_time, MeanOcc>=0.5)

otus_ITSF1_occ0.5_time<-subset(OccAbun_ITSF1_time, MeanOcc>=0.5)
otus_ITSF2_occ0.5_time<-subset(OccAbun_ITSF2_time, MeanOcc>=0.5)
otus_ITSF3_occ0.5_time<-subset(OccAbun_ITSF3_time, MeanOcc>=0.5)

#Reorder by OTU name
otus_16sF1_occ0.5_time<-otus_16sF1_occ0.5_time[order(rownames(otus_16sF1_occ0.5_time)),]
otus_16sF2_occ0.5_time<-otus_16sF2_occ0.5_time[order(rownames(otus_16sF2_occ0.5_time)),]
otus_16sF3_occ0.5_time<-otus_16sF3_occ0.5_time[order(rownames(otus_16sF3_occ0.5_time)),]

otus_ITSF1_occ0.5_time<-otus_ITSF1_occ0.5_time[order(rownames(otus_ITSF1_occ0.5_time)),]
otus_ITSF2_occ0.5_time<-otus_ITSF2_occ0.5_time[order(rownames(otus_ITSF2_occ0.5_time)),]
otus_ITSF3_occ0.5_time<-otus_ITSF3_occ0.5_time[order(rownames(otus_ITSF3_occ0.5_time)),]

#Get taxonomic information for core taxa
otus_16sF1_occ0.5_time_tax<-subset(taxon_16s, rownames(taxon_16s) %in% rownames(otus_16sF1_occ0.5_time))
otus_16sF2_occ0.5_time_tax<-subset(taxon_16s, rownames(taxon_16s) %in% rownames(otus_16sF2_occ0.5_time))
otus_16sF3_occ0.5_time_tax<-subset(taxon_16s, rownames(taxon_16s) %in% rownames(otus_16sF3_occ0.5_time))

otus_ITSF1_occ0.5_time_tax<-subset(taxon_ITS, rownames(taxon_ITS) %in% rownames(otus_ITSF1_occ0.5_time))
otus_ITSF2_occ0.5_time_tax<-subset(taxon_ITS, rownames(taxon_ITS) %in% rownames(otus_ITSF2_occ0.5_time))
otus_ITSF3_occ0.5_time_tax<-subset(taxon_ITS, rownames(taxon_ITS) %in% rownames(otus_ITSF3_occ0.5_time))


#Make list of core OTUs per Facility
F1_16s_net<-rownames(otus_16sF1_occ0.5_time)
F2_16s_net<-rownames(otus_16sF2_occ0.5_time)
F3_16s_net<-rownames(otus_16sF3_occ0.5_time)

F1_ITS_net<-rownames(otus_ITSF1_occ0.5_time)
F2_ITS_net<-rownames(otus_ITSF2_occ0.5_time)
F3_ITS_net<-rownames(otus_ITSF3_occ0.5_time)


#Combine list of core OTUs for all facilities
nettaxa_16s_list<-unique(c(F1_16s_net,F2_16s_net, F3_16s_net)) #1312 nodes
nettaxa_ITS_list<-unique(c(F1_ITS_net,F2_ITS_net, F3_ITS_net)) #495 nodes


#Deal with sparcity. Only include taxa that hae an occupacy above 0.5 in the three facilities from previous analyses.
#One Network per facility, including all samples from the two years



#### Networks  ####
#Install NetCoMi package
devtools::install_github("stefpeschel/NetCoMi", dependencies = TRUE,
                         repos = c("https://cloud.r-project.org/",
                                   BiocManager::repositories()))

library(NetCoMi)

#Install package dependencies
installNetCoMiPacks()

#Make network using OTUs with 0.5 occupancy in all facilities together
#16s
#Construct network with spieceasi
net_16s<-netConstruct(phyloseq_16s_net, 
                          measure='spieceasi',
                          measurePar = list(nlambda=10,lambda.min.ratio=1e-2,pulsar.params=list(rep.num=99)),
                          normMethod = "none",
                          zeroMethod = "none",
                          sparsMethod = "none",
                          verbose = 3,
                          seed=10000)

#Analyze network. Detect hubs by highest betweenness centrality
net_16s_a<-netAnalyze(net_16s,
                      centrLCC = TRUE,
                      clustMethod = "cluster_fast_greedy",
                      hubPar = "betweenness",
                      weightDeg = FALSE, normDeg = FALSE)

#Get network summary
network_summary_16s<-summary(net_16s_a, numbNodes=5L)

#Plot network - color by cluster, highlight hubs
cols16s<-c("#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E")

#Fig 4I
plot_net_16s<-plot(net_16s_a, nodeColor="cluster", 
                   nodeSize="clr", 
                   rmSingles=TRUE, 
                   labels=FALSE,
                   nodeSizeSpread=3, 
                   title1="Bacteria", 
                   highlightHubs=TRUE,
                   showTitle=TRUE, 
                   cexTitle=0.5,
                   repulsion=1.5,
                   nodeTransp=10, 
                   hubTransp=0,
                   edgeTranspLow=90,
                   edgeTranspHigh=50,
                   colorVec=cols16s)
##Save plot as png and svg with the plot tab

#ITS
#Construct network with spieceasi
net_ITS<-netConstruct(phyloseq_ITS_net, 
                      measure='spieceasi',
                      measurePar = list(nlambda=30,lambda.min.ratio=1e-2,pulsar.params=list(rep.num=99)),
                      normMethod = "none",
                      zeroMethod = "none",
                      sparsMethod = "none",
                      verbose = 3,
                      seed=10000)

#Analyze network. Detect hubs by highest betweenness centrality
net_ITS_a<-netAnalyze(net_ITS,
                      centrLCC = TRUE,
                      clustMethod = "cluster_fast_greedy",
                      hubPar = "betweenness",
                      weightDeg = FALSE, normDeg = FALSE)

#Get network summary
network_summary_ITS<-summary(net_ITS_a, numbNodes=5L)


#Plot network - color by cluster, highlight hubs
colsITS<-c("#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02")

#Fig 4J
plot_net_ITS<-plot(net_ITS_a, nodeColor="cluster", 
                   nodeSize="clr", 
                   rmSingles=TRUE, 
                   labels=FALSE,
                   nodeSizeSpread=3, 
                   title1="Fungi", 
                   highlightHubs=TRUE,
                   showTitle=TRUE, 
                   cexTitle=0.5,
                   repulsion=0.8,
                   nodeTransp=0, 
                   hubTransp=0,
                   edgeTranspLow=80,
                   edgeTranspHigh=50,
                   colorVec=colsITS)

##Save plot as png and svg with the plot tab


#Get hub OTU names
hubs_16s<-net_16s_a$hubs$hubs1
hubs_ITS<-net_ITS_a$hubs$hubs1

#Make table with hub taxonomy
hub_16s_taxa<-subset(taxon_16s, rownames(taxon_16s) %in% hubs_16s)
hub_ITS_taxa<-subset(taxon_ITS, rownames(taxon_ITS) %in% hubs_ITS)

write.csv(hub_16s_taxa,"NetworkHubs_16s.csv")
write.csv(hub_ITS_taxa,"NetworkHubs_ITS.csv")



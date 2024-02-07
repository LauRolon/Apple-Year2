#Two-year monitoring of Lm in apple packing houses
#Analysis of core microbiomes and networks at the ASV level
#Laura Rolon
#Last updated: 06/22/22

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
setwd("/storage/work/m/mlr355/Apple/Downstream")

#### IMPORT DATA ####

##16s
asv_16s<-read.csv("ASV_16s_clean.csv", header=TRUE, row.names = 1)
taxon_16s<-as.matrix(read.csv("Taxon_16s_clean.csv", header = TRUE, row.names = 1))
metadata.16s<-read.csv("metadata_16s_clean.csv", header = TRUE, row.names = 1)

##ITS
asv_ITS<-read.csv("ASV_ITS_clean.csv", header=TRUE, row.names = 1)
taxon_ITS<-as.matrix(read.csv("Taxon_ITS_clean.csv", header = TRUE, row.names = 1))
metadata.ITS<-read.csv("metadata_ITS_clean.csv", header = TRUE, row.names = 1)

#### Prepare data for analyses ####
#Make phyloseq
phyloseq16s = phyloseq(otu_table(asv_16s, taxa_are_rows = FALSE), tax_table(taxon_16s), sample_data(metadata.16s))
phyloseqITS = phyloseq(otu_table(asv_ITS, taxa_are_rows = FALSE), tax_table(taxon_ITS), sample_data(metadata.ITS))


# Subset Phyloseq for each facility
physeq_16sF1 <- subset_samples(phyloseq16s, Facility == "F1") 
physeq_16sF2 <- subset_samples(phyloseq16s, Facility == "F2") 
physeq_16sF3 <- subset_samples(phyloseq16s, Facility == "F3") 

physeq_ITSF1 <- subset_samples(phyloseqITS, Facility == "F1") 
physeq_ITSF2 <- subset_samples(phyloseqITS, Facility == "F2") 
physeq_ITSF3 <- subset_samples(phyloseqITS, Facility == "F3") 

#Export ASV table for each facility from Phyloseq 
asv_16sF1<-t(as.data.frame(as(otu_table(physeq_16sF1), "matrix")))
asv_16sF2<-t(as.data.frame(as(otu_table(physeq_16sF2), "matrix")))
asv_16sF3<-t(as.data.frame(as(otu_table(physeq_16sF3), "matrix")))

asv_ITSF1<-t(as.data.frame(as(otu_table(physeq_ITSF1), "matrix")))
asv_ITSF2<-t(as.data.frame(as(otu_table(physeq_ITSF2), "matrix")))
asv_ITSF3<-t(as.data.frame(as(otu_table(physeq_ITSF3), "matrix")))


#Subset metadata
metadata_16sF1<-subset(metadata.16s, Facility=="F1")
metadata_16sF2<-subset(metadata.16s, Facility=="F2")
metadata_16sF3<-subset(metadata.16s, Facility=="F3")

metadata_ITSF1<-subset(metadata.ITS, Facility=="F1")
metadata_ITSF2<-subset(metadata.ITS, Facility=="F2")
metadata_ITSF3<-subset(metadata.ITS, Facility=="F3")

#Remove ASVs that have count zero in all samples - Necessary step for the zero imputation function
asv_16sF1<-asv_16sF1[ which(rowSums(asv_16sF1)>0),]
asv_16sF2<-asv_16sF2[ which(rowSums(asv_16sF2)>0),]
asv_16sF3<-asv_16sF3[ which(rowSums(asv_16sF3)>0),]

asv_ITSF1<-asv_ITSF1[ which(rowSums(asv_ITSF1)>0),]
asv_ITSF2<-asv_ITSF2[ which(rowSums(asv_ITSF2)>0),]
asv_ITSF3<-asv_ITSF3[ which(rowSums(asv_ITSF3)>0),]

#### Occupancy abundance curves for whole data set ####
#Individual sampling sites are considered independent.
#Calculate presence absence
asv_16sF1_PA <- 1*((asv_16sF1>0)==1) 
asv_16sF2_PA <- 1*((asv_16sF2>0)==1) 
asv_16sF3_PA <- 1*((asv_16sF3>0)==1) 

asv_ITSF1_PA <- 1*((asv_ITSF1>0)==1) 
asv_ITSF2_PA <- 1*((asv_ITSF2>0)==1) 
asv_ITSF3_PA <- 1*((asv_ITSF3>0)==1) 

#Calculate mean occupancy for each ASV
asv_16sF1_occ <- as.data.frame(rowSums(asv_16sF1_PA)/ncol(asv_16sF1_PA))
asv_16sF2_occ <- as.data.frame(rowSums(asv_16sF2_PA)/ncol(asv_16sF2_PA))
asv_16sF3_occ <- as.data.frame(rowSums(asv_16sF3_PA)/ncol(asv_16sF3_PA))

asv_ITSF1_occ <- as.data.frame(rowSums(asv_ITSF1_PA)/ncol(asv_ITSF1_PA))
asv_ITSF2_occ <- as.data.frame(rowSums(asv_ITSF2_PA)/ncol(asv_ITSF2_PA))
asv_ITSF3_occ <- as.data.frame(rowSums(asv_ITSF3_PA)/ncol(asv_ITSF3_PA))

#Rename column
asv_16sF1_occ <-asv_16sF1_occ %>% rename(MeanOcc=`rowSums(asv_16sF1_PA)/ncol(asv_16sF1_PA)`)
asv_16sF2_occ <-asv_16sF2_occ %>% rename(MeanOcc=`rowSums(asv_16sF2_PA)/ncol(asv_16sF2_PA)`)
asv_16sF3_occ <-asv_16sF3_occ %>% rename(MeanOcc=`rowSums(asv_16sF3_PA)/ncol(asv_16sF3_PA)`)

asv_ITSF1_occ <-asv_ITSF1_occ %>% rename(MeanOcc=`rowSums(asv_ITSF1_PA)/ncol(asv_ITSF1_PA)`)
asv_ITSF2_occ <-asv_ITSF2_occ %>% rename(MeanOcc=`rowSums(asv_ITSF2_PA)/ncol(asv_ITSF2_PA)`)
asv_ITSF3_occ <-asv_ITSF3_occ %>% rename(MeanOcc=`rowSums(asv_ITSF3_PA)/ncol(asv_ITSF3_PA)`)

#Calculate relative abundance for each ASV by Facility
#First replace zero values with a small non- zero value before clr transformation. 
#Uses Count Zero Multiplicative method to replace zeros and outputs pseudo-counts
asv.n0_16sF1<-t(cmultRepl(t(asv_16sF1), label=0, method="CZM", output="p-counts")) 
asv.n0_16sF2<-t(cmultRepl(t(asv_16sF2), label=0, method="CZM", output="p-counts")) 
asv.n0_16sF3<-t(cmultRepl(t(asv_16sF3), label=0, method="CZM", output="p-counts")) 

asv.n0_ITSF1<-t(cmultRepl(t(asv_ITSF1), label=0, method="CZM", output="p-counts")) 
asv.n0_ITSF2<-t(cmultRepl(t(asv_ITSF2), label=0, method="CZM", output="p-counts")) 
asv.n0_ITSF3<-t(cmultRepl(t(asv_ITSF3), label=0, method="CZM", output="p-counts")) 

#Note: Check the output to make sure there are no negative numbers. If samples or asv are sparse, the CZM method will 
#add a negative number that interferes with the log normalization. If necessary use function below to convert negative values
#into positives
asv_n0_16sF1<-ifelse(asv.n0_16sF1 < 0, asv.n0_16sF1*(-1), asv.n0_16sF1)
asv_n0_16sF2<-ifelse(asv.n0_16sF2 < 0, asv.n0_16sF2*(-1), asv.n0_16sF2)
asv_n0_16sF3<-ifelse(asv.n0_16sF3 < 0, asv.n0_16sF3*(-1), asv.n0_16sF3)

asv_n0_ITSF1<-ifelse(asv.n0_ITSF1 < 0, asv.n0_ITSF1*(-1), asv.n0_ITSF1)
asv_n0_ITSF2<-ifelse(asv.n0_ITSF2 < 0, asv.n0_ITSF2*(-1), asv.n0_ITSF2)
asv_n0_ITSF3<-ifelse(asv.n0_ITSF3 < 0, asv.n0_ITSF3*(-1), asv.n0_ITSF3)


#Then convert to compositions with Aitchinson 
asv.n0.acomp_16sF1<-as.data.frame(acomp(t(asv_n0_16sF1)), total=1)
asv.n0.acomp_16sF2<-as.data.frame(acomp(t(asv_n0_16sF2)), total=1)
asv.n0.acomp_16sF3<-as.data.frame(acomp(t(asv_n0_16sF3)), total=1)

asv.n0.acomp_ITSF1<-as.data.frame(acomp(t(asv_n0_ITSF1)), total=1)
asv.n0.acomp_ITSF2<-as.data.frame(acomp(t(asv_n0_ITSF2)), total=1)
asv.n0.acomp_ITSF3<-as.data.frame(acomp(t(asv_n0_ITSF3)), total=1)


#Calculate the mean relative abundance for each ASV in each facility
asv_16sF1_abun <- as.data.frame(rowSums(t(asv.n0.acomp_16sF1))/ncol(t(asv.n0.acomp_16sF1)))
asv_16sF2_abun <- as.data.frame(rowSums(t(asv.n0.acomp_16sF2))/ncol(t(asv.n0.acomp_16sF2)))
asv_16sF3_abun <- as.data.frame(rowSums(t(asv.n0.acomp_16sF3))/ncol(t(asv.n0.acomp_16sF3)))

asv_ITSF1_abun <- as.data.frame(rowSums(t(asv.n0.acomp_ITSF1))/ncol(t(asv.n0.acomp_ITSF1)))
asv_ITSF2_abun <- as.data.frame(rowSums(t(asv.n0.acomp_ITSF2))/ncol(t(asv.n0.acomp_ITSF2)))
asv_ITSF3_abun <- as.data.frame(rowSums(t(asv.n0.acomp_ITSF3))/ncol(t(asv.n0.acomp_ITSF3)))

#Rename column
asv_16sF1_abun <-asv_16sF1_abun %>% rename(MeanAbun=`rowSums(t(asv.n0.acomp_16sF1))/ncol(t(asv.n0.acomp_16sF1))`)
asv_16sF2_abun <-asv_16sF2_abun %>% rename(MeanAbun=`rowSums(t(asv.n0.acomp_16sF2))/ncol(t(asv.n0.acomp_16sF2))`)
asv_16sF3_abun <-asv_16sF3_abun %>% rename(MeanAbun=`rowSums(t(asv.n0.acomp_16sF3))/ncol(t(asv.n0.acomp_16sF3))`)

asv_ITSF1_abun <-asv_ITSF1_abun %>% rename(MeanAbun=`rowSums(t(asv.n0.acomp_ITSF1))/ncol(t(asv.n0.acomp_ITSF1))`)
asv_ITSF2_abun <-asv_ITSF2_abun %>% rename(MeanAbun=`rowSums(t(asv.n0.acomp_ITSF2))/ncol(t(asv.n0.acomp_ITSF2))`)
asv_ITSF3_abun <-asv_ITSF3_abun %>% rename(MeanAbun=`rowSums(t(asv.n0.acomp_ITSF3))/ncol(t(asv.n0.acomp_ITSF3))`)


#Bind columns
OccAbun_16sF1<-bind_cols(asv_16sF1_occ,asv_16sF1_abun)
OccAbun_16sF2<-bind_cols(asv_16sF2_occ,asv_16sF2_abun)
OccAbun_16sF3<-bind_cols(asv_16sF3_occ,asv_16sF3_abun)

OccAbun_ITSF1<-bind_cols(asv_ITSF1_occ,asv_ITSF1_abun)
OccAbun_ITSF2<-bind_cols(asv_ITSF2_occ,asv_ITSF2_abun)
OccAbun_ITSF3<-bind_cols(asv_ITSF3_occ,asv_ITSF3_abun)


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


#Subset asv with Occupancy=1
asv_16sF1_occ1<-subset(OccAbun_16sF1, MeanOcc==1)%>%
  arrange(rownames())
asv_16sF2_occ1<-subset(OccAbun_16sF2, MeanOcc==1)%>%
  arrange(rownames(asv_16sF2_occ1))
asv_16sF3_occ1<-subset(OccAbun_16sF3, MeanOcc==1)%>%
  arrange(rownames(asv_16sF3_occ1))

asv_ITSF1_occ1<-subset(OccAbun_ITSF1, MeanOcc==1)%>%
  arrange(rownames(asv_ITSF1_occ1))
asv_ITSF2_occ1<-subset(OccAbun_ITSF2, MeanOcc==1)%>%
  arrange(rownames(asv_ITSF2_occ1))
asv_ITSF3_occ1<-subset(OccAbun_ITSF3, MeanOcc==1)%>%
  arrange(rownames(asv_ITSF3_occ1))

#Add taxonomic information
asv_16sF1_occ1_tax<-subset(taxon_16s, rownames(taxon_16s) %in% rownames(asv_16sF1_occ1))
asv_16sF2_occ1_tax<-subset(taxon_16s, rownames(taxon_16s) %in% rownames(asv_16sF2_occ1))
asv_16sF3_occ1_tax<-subset(taxon_16s, rownames(taxon_16s) %in% rownames(asv_16sF3_occ1))

asv_ITSF1_occ1_tax<-subset(taxon_ITS, rownames(taxon_ITS) %in% rownames(asv_ITSF1_occ1))
asv_ITSF2_occ1_tax<-subset(taxon_ITS, rownames(taxon_ITS) %in% rownames(asv_ITSF2_occ1))
asv_ITSF3_occ1_tax<-subset(taxon_ITS, rownames(taxon_ITS) %in% rownames(asv_ITSF3_occ1))

#Bind columns
asv_16sF1_occ1<-bind_cols(asv_16sF1_occ1,asv_16sF1_occ1_tax)
asv_16sF2_occ1<-bind_cols(asv_16sF2_occ1,asv_16sF2_occ1_tax)
asv_16sF3_occ1<-bind_cols(asv_16sF3_occ1,asv_16sF3_occ1_tax)

asv_ITSF1_occ1<-bind_cols(asv_ITSF1_occ1,asv_ITSF1_occ1_tax)
asv_ITSF2_occ1<-bind_cols(asv_ITSF2_occ1,asv_ITSF2_occ1_tax)
asv_ITSF3_occ1<-bind_cols(asv_ITSF3_occ1,asv_ITSF3_occ1_tax)




#### Occupancy - abundance curves by collection time point ####
#Calculate mean occupancy by collection time. 
#We consider that the the samples collected at each facility each time are technical replicates.
#We want to identify those ASVs that were present in each facility at every timepoint throughout two years

#Merge samples taken in the same collection date
phyloseq_16sF1_time<- merge_samples(physeq_16sF1, group = "Week")
phyloseq_16sF2_time<- merge_samples(physeq_16sF2, group = "Week")
phyloseq_16sF3_time<- merge_samples(physeq_16sF3, group = "Week")

phyloseq_ITSF1_time<- merge_samples(physeq_ITSF1, group = "Week")
phyloseq_ITSF2_time<- merge_samples(physeq_ITSF2, group = "Week")
phyloseq_ITSF3_time<- merge_samples(physeq_ITSF3, group = "Week")

#Export ASV table for each facility from Phyloseq 
asv_16sF1_time<-t(as.data.frame(as(otu_table(phyloseq_16sF1_time), "matrix")))
asv_16sF2_time<-t(as.data.frame(as(otu_table(phyloseq_16sF2_time), "matrix")))
asv_16sF3_time<-t(as.data.frame(as(otu_table(phyloseq_16sF3_time), "matrix")))

asv_ITSF1_time<-t(as.data.frame(as(otu_table(phyloseq_ITSF1_time), "matrix")))
asv_ITSF2_time<-t(as.data.frame(as(otu_table(phyloseq_ITSF2_time), "matrix")))
asv_ITSF3_time<-t(as.data.frame(as(otu_table(phyloseq_ITSF3_time), "matrix")))

#Remove ASVs that have count zero in all samples - Necessary step for the zero imputation function
asv_16sF1_time<-asv_16sF1_time[ which(rowSums(asv_16sF1_time)>0),]
asv_16sF2_time<-asv_16sF2_time[ which(rowSums(asv_16sF2_time)>0),]
asv_16sF3_time<-asv_16sF3_time[ which(rowSums(asv_16sF3_time)>0),]

asv_ITSF1_time<-asv_ITSF1_time[ which(rowSums(asv_ITSF1_time)>0),]
asv_ITSF2_time<-asv_ITSF2_time[ which(rowSums(asv_ITSF2_time)>0),]
asv_ITSF3_time<-asv_ITSF3_time[ which(rowSums(asv_ITSF3_time)>0),]

# #Transpose asv table to have ASVs in rows
# asv_16sF1_time<-t(asv_16sF1_time)
# asv_16sF2_time<-t(asv_16sF2_time)
# asv_16sF3_time<-t(asv_16sF3_time)
# 
# asv_ITSF1_time<-t(asv_ITSF1_time)
# asv_ITSF2_time<-t(asv_ITSF2_time)
# asv_ITSF3_time<-t(asv_ITSF3_time)

#Calculate presence absence
asv_16sF1_time_PA <- 1*((asv_16sF1_time>0)==1) 
asv_16sF2_time_PA <- 1*((asv_16sF2_time>0)==1) 
asv_16sF3_time_PA <- 1*((asv_16sF3_time>0)==1) 

asv_ITSF1_time_PA <- 1*((asv_ITSF1_time>0)==1) 
asv_ITSF2_time_PA <- 1*((asv_ITSF2_time>0)==1) 
asv_ITSF3_time_PA <- 1*((asv_ITSF3_time>0)==1) 

#Calculate mean occupancy for each ASV
asv_16sF1_time_occ <- as.data.frame(rowSums(asv_16sF1_time_PA)/ncol(asv_16sF1_time_PA))
asv_16sF2_time_occ <- as.data.frame(rowSums(asv_16sF2_time_PA)/ncol(asv_16sF2_time_PA))
asv_16sF3_time_occ <- as.data.frame(rowSums(asv_16sF3_time_PA)/ncol(asv_16sF3_time_PA))

asv_ITSF1_time_occ <- as.data.frame(rowSums(asv_ITSF1_time_PA)/ncol(asv_ITSF1_time_PA))
asv_ITSF2_time_occ <- as.data.frame(rowSums(asv_ITSF2_time_PA)/ncol(asv_ITSF2_time_PA))
asv_ITSF3_time_occ <- as.data.frame(rowSums(asv_ITSF3_time_PA)/ncol(asv_ITSF3_time_PA))

#Rename column
asv_16sF1_time_occ <-asv_16sF1_time_occ %>% rename(MeanOcc=`rowSums(asv_16sF1_time_PA)/ncol(asv_16sF1_time_PA)`)
asv_16sF2_time_occ <-asv_16sF2_time_occ %>% rename(MeanOcc=`rowSums(asv_16sF2_time_PA)/ncol(asv_16sF2_time_PA)`)
asv_16sF3_time_occ <-asv_16sF3_time_occ %>% rename(MeanOcc=`rowSums(asv_16sF3_time_PA)/ncol(asv_16sF3_time_PA)`)

asv_ITSF1_time_occ <-asv_ITSF1_time_occ %>% rename(MeanOcc=`rowSums(asv_ITSF1_time_PA)/ncol(asv_ITSF1_time_PA)`)
asv_ITSF2_time_occ <-asv_ITSF2_time_occ %>% rename(MeanOcc=`rowSums(asv_ITSF2_time_PA)/ncol(asv_ITSF2_time_PA)`)
asv_ITSF3_time_occ <-asv_ITSF3_time_occ %>% rename(MeanOcc=`rowSums(asv_ITSF3_time_PA)/ncol(asv_ITSF3_time_PA)`)


#Calculate relative abundance for each ASV by Facility
#First replace zero values with a small non- zero value before clr transformation. 
#Uses Count Zero Multiplicative method to replace zeros and outputs pseudo-counts
asv.n0_time_16sF1<-t(cmultRepl(t(asv_16sF1_time), label=0, method="CZM", output="p-counts")) 
asv.n0_time_16sF2<-t(cmultRepl(t(asv_16sF2_time), label=0, method="CZM", output="p-counts"))
asv.n0_time_16sF3<-t(cmultRepl(t(asv_16sF3_time), label=0, method="CZM", output="p-counts"))

asv.n0_time_ITSF1<-t(cmultRepl(t(asv_ITSF1_time), label=0, method="CZM", output="p-counts"))
asv.n0_time_ITSF2<-t(cmultRepl(t(asv_ITSF2_time), label=0, method="CZM", output="p-counts"))
asv.n0_time_ITSF3<-t(cmultRepl(t(asv_ITSF3_time), label=0, method="CZM", output="p-counts")) 

#Note: Check the output to make sure there are no negative numbers. If samples or asv are sparse, the CZM method will 
#add a negative number that interferes with the log normalization. If necessary use function below to convert negative values
#into positives
asv_n0_time_16sF1<-ifelse(asv.n0_time_16sF1 < 0, asv.n0_time_16sF1*(-1), asv.n0_time_16sF1)
asv_n0_time_16sF2<-ifelse(asv.n0_time_16sF2 < 0, asv.n0_time_16sF2*(-1), asv.n0_time_16sF2)
asv_n0_time_16sF3<-ifelse(asv.n0_time_16sF3 < 0, asv.n0_time_16sF3*(-1), asv.n0_time_16sF3)

asv_n0_time_ITSF1<-ifelse(asv.n0_time_ITSF1 < 0, asv.n0_time_ITSF1*(-1), asv.n0_time_ITSF1)
asv_n0_time_ITSF2<-ifelse(asv.n0_time_ITSF2 < 0, asv.n0_time_ITSF2*(-1), asv.n0_time_ITSF2)
asv_n0_time_ITSF3<-ifelse(asv.n0_time_ITSF3 < 0, asv.n0_time_ITSF3*(-1), asv.n0_time_ITSF3)


#Then convert to compositions with Aitchinson 
asv.n0.acomp_16sF1_time<-t(as.data.frame(acomp(t(asv_n0_time_16sF1)), total=1))
asv.n0.acomp_16sF2_time<-t(as.data.frame(acomp(t(asv_n0_time_16sF2)), total=1))
asv.n0.acomp_16sF3_time<-t(as.data.frame(acomp(t(asv_n0_time_16sF3)), total=1))

asv.n0.acomp_ITSF1_time<-t(as.data.frame(acomp(t(asv_n0_time_ITSF1)), total=1))
asv.n0.acomp_ITSF2_time<-t(as.data.frame(acomp(t(asv_n0_time_ITSF2)), total=1))
asv.n0.acomp_ITSF3_time<-t(as.data.frame(acomp(t(asv_n0_time_ITSF3)), total=1))

#Calculate the mean relative abundance for each ASV in each facility
asv_16sF1_abun_time <- as.data.frame(rowSums(asv.n0.acomp_16sF1_time)/ncol(asv.n0.acomp_16sF1_time))
asv_16sF2_abun_time <- as.data.frame(rowSums(asv.n0.acomp_16sF2_time)/ncol(asv.n0.acomp_16sF2_time))
asv_16sF3_abun_time <- as.data.frame(rowSums(asv.n0.acomp_16sF3_time)/ncol(asv.n0.acomp_16sF3_time))

asv_ITSF1_abun_time <- as.data.frame(rowSums(asv.n0.acomp_ITSF1_time)/ncol(asv.n0.acomp_ITSF1_time))
asv_ITSF2_abun_time <- as.data.frame(rowSums(asv.n0.acomp_ITSF2_time)/ncol(asv.n0.acomp_ITSF2_time))
asv_ITSF3_abun_time <- as.data.frame(rowSums(asv.n0.acomp_ITSF3_time)/ncol(asv.n0.acomp_ITSF3_time))

#Rename column
asv_16sF1_abun_time <-asv_16sF1_abun_time %>% rename(MeanAbun=`rowSums(asv.n0.acomp_16sF1_time)/ncol(asv.n0.acomp_16sF1_time)`)
asv_16sF2_abun_time <-asv_16sF2_abun_time %>% rename(MeanAbun=`rowSums(asv.n0.acomp_16sF2_time)/ncol(asv.n0.acomp_16sF2_time)`)
asv_16sF3_abun_time <-asv_16sF3_abun_time %>% rename(MeanAbun=`rowSums(asv.n0.acomp_16sF3_time)/ncol(asv.n0.acomp_16sF3_time)`)

asv_ITSF1_abun_time <-asv_ITSF1_abun_time %>% rename(MeanAbun=`rowSums(asv.n0.acomp_ITSF1_time)/ncol(asv.n0.acomp_ITSF1_time)`)
asv_ITSF2_abun_time <-asv_ITSF2_abun_time %>% rename(MeanAbun=`rowSums(asv.n0.acomp_ITSF2_time)/ncol(asv.n0.acomp_ITSF2_time)`)
asv_ITSF3_abun_time <-asv_ITSF3_abun_time %>% rename(MeanAbun=`rowSums(asv.n0.acomp_ITSF3_time)/ncol(asv.n0.acomp_ITSF3_time)`)

#Bind columns
OccAbun_16sF1_time<-bind_cols(asv_16sF1_time_occ,asv_16sF1_abun_time)
OccAbun_16sF2_time<-bind_cols(asv_16sF2_time_occ,asv_16sF2_abun_time)
OccAbun_16sF3_time<-bind_cols(asv_16sF3_time_occ,asv_16sF3_abun_time)

OccAbun_ITSF1_time<-bind_cols(asv_ITSF1_time_occ,asv_ITSF1_abun_time)
OccAbun_ITSF2_time<-bind_cols(asv_ITSF2_time_occ,asv_ITSF2_abun_time)
OccAbun_ITSF3_time<-bind_cols(asv_ITSF3_time_occ,asv_ITSF3_abun_time)

#Add color column for ASVs with occurrence=1
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
        plot.background = element_rect(fill="white", color =NA)) +
  theme(panel.border = element_rect(color="black", fill=NA))+
  theme(legend.position = 'none') +
  ggtitle("Bacteria F1")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=18))+
  scale_color_manual(values=c("#00000080","#420A6880"))
  #theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))

Plot_OccAbun_16sF2_time<-ggplot(OccAbun_16sF2_time, aes(x=MeanAbun, y=MeanOcc, color=color))+geom_point()+scale_x_log10()+
  scale_size(limits=c(0,50))+
  ylab('Occupancy')+xlab("log10 Mean Relative Abundance")+
  theme(axis.title = element_text(color='black', size=17), axis.text=element_text(color='black', size=16), 
        axis.ticks=element_line(color='black'), panel.grid.major.y = element_line(color='grey50'))+
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  theme(panel.background = element_rect(fill="white", color =NA),
        plot.background = element_rect(fill="white", color =NA)) +
  theme(panel.border = element_rect(color="black", fill=NA))+
  theme(legend.position = 'none') +
  ggtitle("Bacteria F2")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=18))+
  scale_color_manual(values=c("#00000080","#BB375480"))
  #theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))

Plot_OccAbun_16sF3_time<-ggplot(OccAbun_16sF3_time, aes(x=MeanAbun, y=MeanOcc, color=color))+geom_point()+scale_x_log10()+
  scale_size(limits=c(0,50))+
  ylab('Occupancy')+xlab("log10 Mean Relative Abundance")+
  theme(axis.title = element_text(color='black', size=17), axis.text=element_text(color='black', size=16), 
        axis.ticks=element_line(color='black'), panel.grid.major.y = element_line(color='grey50'))+
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  theme(panel.background = element_rect(fill="white", color =NA),
        plot.background = element_rect(fill="white", color =NA)) +
  theme(panel.border = element_rect(color="black", fill=NA))+
  theme(legend.position = 'none') +
  ggtitle("Bacteria F3")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=18))+
  scale_color_manual(values=c("#00000080","#FCA50A80"))
  #theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))

Plot_OccAbun_ITSF1_time<-ggplot(OccAbun_ITSF1_time, aes(x=MeanAbun, y=MeanOcc, color=color))+geom_point()+scale_x_log10()+
  scale_size(limits=c(0,50))+
  ylab('Occupancy')+xlab("log10 Mean Relative Abundance")+
  theme(axis.title = element_text(color='black', size=17), axis.text=element_text(color='black', size=16), 
        axis.ticks=element_line(color='black'), panel.grid.major.y = element_line(color='grey50'))+
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  theme(panel.background = element_rect(fill="white", color =NA),
        plot.background = element_rect(fill="white", color =NA)) +
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
        plot.background = element_rect(fill="white", color =NA)) +
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
        plot.background = element_rect(fill="white", color =NA)) +
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


#Subset ASVs with Occupancy=1
asv_16sF1_occ1_time<-subset(OccAbun_16sF1_time, MeanOcc==1)
asv_16sF2_occ1_time<-subset(OccAbun_16sF2_time, MeanOcc==1)
asv_16sF3_occ1_time<-subset(OccAbun_16sF3_time, MeanOcc==1)

asv_ITSF1_occ1_time<-subset(OccAbun_ITSF1_time, MeanOcc==1)
asv_ITSF2_occ1_time<-subset(OccAbun_ITSF2_time, MeanOcc==1)
asv_ITSF3_occ1_time<-subset(OccAbun_ITSF3_time, MeanOcc==1)

#Reorder by ASV name
asv_16sF1_occ1_time<-asv_16sF1_occ1_time[order(rownames(asv_16sF1_occ1_time)),]
asv_16sF2_occ1_time<-asv_16sF2_occ1_time[order(rownames(asv_16sF2_occ1_time)),]
asv_16sF3_occ1_time<-asv_16sF3_occ1_time[order(rownames(asv_16sF3_occ1_time)),]

asv_ITSF1_occ1_time<-asv_ITSF1_occ1_time[order(rownames(asv_ITSF1_occ1_time)),]
asv_ITSF2_occ1_time<-asv_ITSF2_occ1_time[order(rownames(asv_ITSF2_occ1_time)),]
asv_ITSF3_occ1_time<-asv_ITSF3_occ1_time[order(rownames(asv_ITSF3_occ1_time)),]

#Get taxonomic information for core taxa
asv_16sF1_occ1_time_tax<-subset(taxon_16s, rownames(taxon_16s) %in% rownames(asv_16sF1_occ1_time))
asv_16sF2_occ1_time_tax<-subset(taxon_16s, rownames(taxon_16s) %in% rownames(asv_16sF2_occ1_time))
asv_16sF3_occ1_time_tax<-subset(taxon_16s, rownames(taxon_16s) %in% rownames(asv_16sF3_occ1_time))

asv_ITSF1_occ1_time_tax<-subset(taxon_ITS, rownames(taxon_ITS) %in% rownames(asv_ITSF1_occ1_time))
asv_ITSF2_occ1_time_tax<-subset(taxon_ITS, rownames(taxon_ITS) %in% rownames(asv_ITSF2_occ1_time))
asv_ITSF3_occ1_time_tax<-subset(taxon_ITS, rownames(taxon_ITS) %in% rownames(asv_ITSF3_occ1_time))

#Bind columns
asv_16sF1_occ1_time<-bind_cols(asv_16sF1_occ1_time,asv_16sF1_occ1_time_tax)
asv_16sF2_occ1_time<-bind_cols(asv_16sF2_occ1_time,asv_16sF2_occ1_time_tax)
asv_16sF3_occ1_time<-bind_cols(asv_16sF3_occ1_time,asv_16sF3_occ1_time_tax)

asv_ITSF1_occ1_time<-bind_cols(asv_ITSF1_occ1_time,asv_ITSF1_occ1_time_tax)
asv_ITSF2_occ1_time<-bind_cols(asv_ITSF2_occ1_time,asv_ITSF2_occ1_time_tax)
asv_ITSF3_occ1_time<-bind_cols(asv_ITSF3_occ1_time,asv_ITSF3_occ1_time_tax)


#Make list of core ASVs per Facility
F1_16s_core<-rownames(asv_16sF1_occ1_time)
F2_16s_core<-rownames(asv_16sF2_occ1_time)
F3_16s_core<-rownames(asv_16sF3_occ1_time)

F1_ITS_core<-rownames(asv_ITSF1_occ1_time)
F2_ITS_core<-rownames(asv_ITSF2_occ1_time)
F3_ITS_core<-rownames(asv_ITSF3_occ1_time)


#Combine list of core ASVs for all facilities
coretaxa_16s_list<-list(F1=F1_16s_core, F2=F2_16s_core, F3=F3_16s_core)
coretaxa_ITS_list<-list(F1=F1_ITS_core, F2=F2_ITS_core, F3=F3_ITS_core)

#Make dataframe for all core taxa
core_16s<-data.frame(OTU=c(F1_16s_core, F2_16s_core, F3_16s_core), Facility=c(rep("F1",69),rep("F2",139),rep("F3",71)))
core_ITS<-data.frame(OTU=c(F1_ITS_core, F2_ITS_core, F3_ITS_core), Facility=c(rep("F1",38),rep("F2",31),rep("F3",38)))

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
Venn_core = plot_grid(venn_16s, 
                      #venn_ITS, 
                         ncol=2, nrow=1, labels = c("G","H"), label_size = 20, vjust = 2, hjust = -1.5)
Venn_core

ggsave("Venn_core_16s.png", plot=Venn_core, device="png", width=8, height=4, units="in", dpi=600)
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


#### Network analyses ####
#Subset asv with Occupancy>0.5
asv_16sF1_occ0.5_time<-subset(OccAbun_16sF1_time, MeanOcc>=0.5)
asv_16sF2_occ0.5_time<-subset(OccAbun_16sF2_time, MeanOcc>=0.5)
asv_16sF3_occ0.5_time<-subset(OccAbun_16sF3_time, MeanOcc>=0.5)

asv_ITSF1_occ0.5_time<-subset(OccAbun_ITSF1_time, MeanOcc>=0.5)
asv_ITSF2_occ0.5_time<-subset(OccAbun_ITSF2_time, MeanOcc>=0.5)
asv_ITSF3_occ0.5_time<-subset(OccAbun_ITSF3_time, MeanOcc>=0.5)

#Reorder by ASV name
asv_16sF1_occ0.5_time<-asv_16sF1_occ0.5_time[order(rownames(asv_16sF1_occ0.5_time)),]
asv_16sF2_occ0.5_time<-asv_16sF2_occ0.5_time[order(rownames(asv_16sF2_occ0.5_time)),]
asv_16sF3_occ0.5_time<-asv_16sF3_occ0.5_time[order(rownames(asv_16sF3_occ0.5_time)),]

asv_ITSF1_occ0.5_time<-asv_ITSF1_occ0.5_time[order(rownames(asv_ITSF1_occ0.5_time)),]
asv_ITSF2_occ0.5_time<-asv_ITSF2_occ0.5_time[order(rownames(asv_ITSF2_occ0.5_time)),]
asv_ITSF3_occ0.5_time<-asv_ITSF3_occ0.5_time[order(rownames(asv_ITSF3_occ0.5_time)),]

#Get taxonomic information for core taxa
asv_16sF1_occ0.5_time_tax<-subset(taxon_16s, rownames(taxon_16s) %in% rownames(asv_16sF1_occ0.5_time))
asv_16sF2_occ0.5_time_tax<-subset(taxon_16s, rownames(taxon_16s) %in% rownames(asv_16sF2_occ0.5_time))
asv_16sF3_occ0.5_time_tax<-subset(taxon_16s, rownames(taxon_16s) %in% rownames(asv_16sF3_occ0.5_time))

asv_ITSF1_occ0.5_time_tax<-subset(taxon_ITS, rownames(taxon_ITS) %in% rownames(asv_ITSF1_occ0.5_time))
asv_ITSF2_occ0.5_time_tax<-subset(taxon_ITS, rownames(taxon_ITS) %in% rownames(asv_ITSF2_occ0.5_time))
asv_ITSF3_occ0.5_time_tax<-subset(taxon_ITS, rownames(taxon_ITS) %in% rownames(asv_ITSF3_occ0.5_time))


#Make list of core ASVs per Facility
F1_16s_net<-rownames(asv_16sF1_occ0.5_time)
F2_16s_net<-rownames(asv_16sF2_occ0.5_time)
F3_16s_net<-rownames(asv_16sF3_occ0.5_time)

F1_ITS_net<-rownames(asv_ITSF1_occ0.5_time)
F2_ITS_net<-rownames(asv_ITSF2_occ0.5_time)
F3_ITS_net<-rownames(asv_ITSF3_occ0.5_time)


#Combine list of core ASVs for all facilities
nettaxa_16s_list<-unique(c(F1_16s_net,F2_16s_net, F3_16s_net)) #1785 nodes
nettaxa_ITS_list<-unique(c(F1_ITS_net,F2_ITS_net, F3_ITS_net)) #348 nodes

asv_16s.T<-t(asv_16s)
asv_ITS.T<-t(asv_ITS)

asvs_16s_net<-subset(asv_16s.T, rownames(asv_16s.T) %in% nettaxa_16s_list)
asvs_ITS_net<-subset(asv_ITS.T, rownames(asv_ITS.T) %in% nettaxa_ITS_list)


#Deal with sparcity. Only include taxa that have an occupacy above 0.5 in the three facilities from previous analyses.
#One Network per facility, including all samples from the two years

#Make phyloseq for Network analyses
phyloseq_16s_net<-phyloseq(otu_table(asvs_16s_net, taxa_are_rows = TRUE), tax_table(taxon_16s), sample_data(metadata.16s))
phyloseq_ITS_net<-phyloseq(otu_table(asvs_ITS_net, taxa_are_rows = TRUE), tax_table(taxon_ITS), sample_data(metadata.ITS))

#### Networks  ####
#Install NetCoMi package
devtools::install_github("stefpeschel/NetCoMi", dependencies = TRUE,
                         repos = c("https://cloud.r-project.org/",
                                   BiocManager::repositories()))

library(NetCoMi)

#Install package dependencies
installNetCoMiPacks()

#Make network using ASVs with 0.5 occupancy in all facilities together
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
cols16s<-c("#F37B59","#D89000","#AFA100","#00BF7D","#00BF4F","#BF80FF","#FF689F")

#Fig 4I

##Save plot as png and svg with the plot tab

png(file="Network_Bacteria2.png", width=500, height=400, res=600)
plot(net_16s_a, nodeColor="cluster", 
                   nodeSize="clr", 
                   rmSingles=TRUE, 
                   #labels=NULL,
                    labelScale=TRUE,
                   nodeSizeSpread=1, 
                   highlightHubs=TRUE,
                   showTitle=FALSE, 
     cexLabels=0,
                        repulsion=1.5,
                   nodeTransp=10, 
                   hubTransp=0,
                   edgeTranspLow=90,
                   edgeTranspHigh=50,
                   colorVec=cols16s)
dev.off()

svg(file="Network_Bacteria.svg", width=500, height=400)
plot(net_16s_a, nodeColor="cluster", 
     nodeSize="clr", 
     rmSingles=TRUE, 
     #labels=NULL,
     labelScale=TRUE,
     nodeSizeSpread=1, 
     highlightHubs=TRUE,
     showTitle=FALSE, 
     cexLabels=0,
     repulsion=1.5,
     nodeTransp=10, 
     hubTransp=0,
     edgeTranspLow=90,
     edgeTranspHigh=50,
     colorVec=cols16s)
dev.off()


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
png(file="Network_Fungi.png", width=500, height=400, res=600)
plot(net_ITS_a, nodeColor="cluster", 
     nodeSize="clr", 
     rmSingles=TRUE, 
     #labels=NULL,
     labelScale=TRUE,
     nodeSizeSpread=1, 
     highlightHubs=TRUE,
     showTitle=FALSE, 
     cexLabels=0,
     repulsion=1.5,
     nodeTransp=10, 
     hubTransp=0,
     edgeTranspLow=90,
     edgeTranspHigh=50,
     colorVec=colsITS)
dev.off()

svg(file="Network_Fungi.svg", width=500, height=400)
plot(net_ITS_a, nodeColor="cluster", 
     nodeSize="clr", 
     rmSingles=TRUE, 
     #labels=NULL,
     labelScale=TRUE,
     nodeSizeSpread=1, 
     highlightHubs=TRUE,
     showTitle=FALSE, 
     cexLabels=0,
     repulsion=1.5,
     nodeTransp=10, 
     hubTransp=0,
     edgeTranspLow=90,
     edgeTranspHigh=50,
     colorVec=colsITS)
dev.off()

##Save plot as png and svg with the plot tab


#Get hub ASV names
hubs_16s<-net_16s_a$hubs$hubs1
hubs_ITS<-net_ITS_a$hubs$hubs1

#Make table with hub taxonomy
hub_16s_taxa<-subset(taxon_16s, rownames(taxon_16s) %in% hubs_16s)
hub_ITS_taxa<-subset(taxon_ITS, rownames(taxon_ITS) %in% hubs_ITS)

write.csv(hub_16s_taxa,"NetworkHubs_16s_2.csv")
write.csv(hub_ITS_taxa,"NetworkHubs_ITS.csv")



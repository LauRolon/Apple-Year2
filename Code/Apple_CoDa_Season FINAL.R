#Two-year monitoring of Lm in apple packing houses
#Analysis of microbiomes and mycobiomes using CoDa approach
#Laura Rolon
#Last updated: 08/26/21

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
library(dendextend)
library(readxl)
library(gplots)
library(pairwiseAdonis)
library(SpadeR)
library(psych)
library(plotly)
library(htmlwidgets)
library(svglite)

set.seed(336)

#Set working directory to where files are located
setwd("C:/Users/lau_r/Google Drive/Penn State/Research/File for R/Apple/CoDa/By season")

#### IMPORT DATA ####

#Import OTU table - 16s (SILVA version 132)
otus_16s<-as.data.frame(import_mothur(mothur_shared_file = 'apple16s_132.shared'))

#Import taxonomy table -16s (SILVA version 132)
taxon_16s <- as.data.frame(import_mothur(mothur_constaxonomy_file = 'apple16s_132.cons.taxonomy'))
colnames(taxon_16s) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus")
TAX_16s =tax_table(as.matrix(taxon_16s))

#Import metadata -16s
metadata_16s <-read.csv("metadata_apple_16s.csv", header=TRUE, row.names=1)
META_16s = sample_data(metadata_16s)


#Import OTU table - ITS
otus_ITS<-as.data.frame(import_mothur(mothur_shared_file = 'appleITS.shared'))

#Import taxonomy table -ITS
taxon_ITS <- as.data.frame(import_mothur(mothur_constaxonomy_file = 'appleITS.cons.taxonomy'))
colnames(taxon_ITS) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus")
TAX_ITS =tax_table(as.matrix(taxon_ITS))

#Import metadata -ITS
metadata_ITS <-read.csv("metadata_apple_ITS.csv", header=TRUE, row.names=1)
META_ITS = sample_data(metadata_ITS)



#### RAREFACTION CURVES ####
#Figure S2
#Add ggrare function to R. Code available at https://rdrr.io/github/gauravsk/ranacapa/src/R/ggrare.R
phyloseq_rare_apple16sY1 = phyloseq(otu_table(otus_16sY1, taxa_are_rows = TRUE), TAX_16sY1, META_16sY1)
phyloseq_rare_appleITSY1 = phyloseq(otu_table(otus_ITSY1, taxa_are_rows = TRUE), TAX_ITSY1, META_ITSY1)
phyloseq_rare_apple16sY2 = phyloseq(otu_table(otus_16sY2, taxa_are_rows = TRUE), TAX_16sY2, META_16sY2)
phyloseq_rare_appleITSY2 = phyloseq(otu_table(otus_ITSY2, taxa_are_rows = TRUE), TAX_ITSY2, META_ITSY2)
rare_apple16sY1 <- ggrare(phyloseq_rare_apple16sY1, step = 100, se=TRUE, color="Facility")
rare_appleITSY1 <- ggrare(phyloseq_rare_appleITSY1, step = 100, se=TRUE, color="Facility")
rare_apple16sY2 <- ggrare(phyloseq_rare_apple16sY2, step = 100, se=TRUE, color="Facility")
rare_appleITSY2 <- ggrare(phyloseq_rare_appleITSY2, step = 100, se=TRUE, color="Facility")

#Plot by facility and year
rare_16sY1_byfacility_plot_full <- rare_apple16sY1 + facet_grid(Facility~.) + 
  theme(strip.background = element_rect(color='black', fill='white', size=0.8, linetype = 'solid'))+
  xlab("Number of reads (x 10,000)") + ylab("Number of unique OTUs") +
  scale_x_continuous(breaks= seq(0,500000, 10000), labels = function(x){x/10000}) + ylim(0,3000)+
  theme(axis.text.x = element_text(angle=90), axis.text = element_text(color='black', size=8)) +
  theme(panel.background = element_rect(fill='transparent', color = NA),
        plot.background = element_rect(fill = 'transparent',color = NA), 
        panel.border = element_rect(color='black',fill = NA))+
  theme(panel.grid.major = element_line(color='gray'))+
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend= -Inf) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend= -Inf) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend= -Inf) +
  scale_fill_viridis_d(begin=0.2, end=0.8, option='inferno')+
  scale_color_viridis_d(begin=0.2, end=0.8, option='inferno')+
  theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))
rare_16sY1_byfacility_plot_full
ggsave("Rarefaction_16sY1.png", plot =rare_16sY1_byfacility_plot_full, device="png", width=5, height=4, units="in",dpi=600)

rare_ITSY1_byfacility_plot_full <- rare_appleITSY1 + facet_grid(Facility~.) + 
  theme(strip.background = element_rect(color='black', fill='white', size=0.8, linetype = 'solid'))+
  xlab("Number of reads (x 10,000)") + ylab("Number of unique OTUs") +
  scale_x_continuous(breaks= seq(0,500000, 10000), labels = function(x){x/10000}) + ylim(0,1500)+
  theme(axis.text.x = element_text(angle=90), axis.text = element_text(color='black', size=8)) +
  theme(panel.background = element_rect(fill='transparent', color = NA),
        plot.background = element_rect(fill = 'transparent',color = NA), 
        panel.border = element_rect(color='black',fill = NA))+
  theme(panel.grid.major = element_line(color='gray'))+
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend= -Inf) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend= -Inf) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend= -Inf) +
  scale_fill_viridis_d(begin=0.2, end=0.8, option='inferno')+
  scale_color_viridis_d(begin=0.2, end=0.8, option='inferno')+
  theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))
rare_ITSY1_byfacility_plot_full
ggsave("Rarefaction_ITSY1.png", plot =rare_ITSY1_byfacility_plot_full, device="png", width=5, height=4, units="in",dpi=600)

rare_16sY2_byfacility_plot_full <- rare_apple16sY2 + facet_grid(Facility~.) + 
  theme(strip.background = element_rect(color='black', fill='white', size=0.8, linetype = 'solid'))+
  xlab("Number of reads (x 10,000)") + ylab("Number of unique OTUs") +
  scale_x_continuous(breaks= seq(0,500000, 20000), labels = function(x){x/10000}) + ylim(0,3000)+
  theme(axis.text.x = element_text(angle=90), axis.text = element_text(color='black', size=8)) +
  theme(panel.background = element_rect(fill='transparent', color = NA),
        plot.background = element_rect(fill = 'transparent',color = NA), 
        panel.border = element_rect(color='black',fill = NA))+
  theme(panel.grid.major = element_line(color='gray'))+
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend= -Inf) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend= -Inf) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend= -Inf) +
  scale_fill_viridis_d(begin=0.2, end=0.8, option='inferno')+
  scale_color_viridis_d(begin=0.2, end=0.8, option='inferno')+
  theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))
rare_16sY2_byfacility_plot_full
ggsave("Rarefaction_16sY2.png", plot =rare_16sY2_byfacility_plot_full, device="png", width=5, height=4, units="in",dpi=600)

rare_ITSY2_byfacility_plot_full <- rare_appleITSY2 + facet_grid(Facility~.) + 
  theme(strip.background = element_rect(color='black', fill='white', size=0.8, linetype = 'solid'))+
  xlab("Number of reads (x 10,000)") + ylab("Number of unique OTUs") +
  scale_x_continuous(breaks= seq(0,500000, 10000), labels = function(x){x/10000}) + ylim(0,1500)+
  theme(axis.text.x = element_text(angle=90), axis.text = element_text(color='black', size=8)) +
  theme(panel.background = element_rect(fill='transparent', color = NA),
        plot.background = element_rect(fill = 'transparent',color = NA), 
        panel.border = element_rect(color='black',fill = NA))+
  theme(panel.grid.major = element_line(color='gray'))+
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend= -Inf) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend= -Inf) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend= -Inf) +
  scale_fill_viridis_d(begin=0.2, end=0.8, option='inferno')+
  scale_color_viridis_d(begin=0.2, end=0.8, option='inferno')+
  theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))
rare_ITSY2_byfacility_plot_full
ggsave("Rarefaction_ITSY2.png", plot =rare_ITSY2_byfacility_plot_full, device="png", width=5, height=4, units="in",dpi=600)

#Combine rarefaction plots
Rarefaction = plot_grid(rare_16sY1_byfacility_plot_full+ theme(legend.position = "none"),
                     rare_16sY2_byfacility_plot_full+ theme(legend.position = "none"),
                     rare_ITSY1_byfacility_plot_full+ theme(legend.position = "none"),
                     rare_ITSY2_byfacility_plot_full+ theme(legend.position = "none"),
                     ncol=2, nrow=2, labels = c("A","B","C","D"), label_size = 20)
Rarefaction
ggsave("Rarefaction_combined.png", plot =Rarefaction, device="png", width=8, height=10, units="in",dpi=600)

#Zoomed in Rarefaction plots
rare_16sY1_cropped <- rare_apple16sY1 + facet_grid(Facility~.) + 
  theme(strip.background = element_rect(color='black', fill='white', size=0.8, linetype = 'solid'))+
  xlab("Number of reads (x 10,000)") + ylab("Number of unique OTUs") +
  scale_x_continuous(limits=c(0,50000), breaks= seq(0,50000, 5000) , labels = function(x){x/10000}) + ylim(0,2000)+
  theme(axis.text.x = element_text(size=8, color='black'), axis.text.y = element_text(size=8, color='black'),
        axis.text=element_text(size=8)) + 
  theme(panel.background = element_rect(fill='transparent', color = NA),
        plot.background = element_rect(fill = 'transparent',color = NA), 
        panel.border = element_rect(color='black',fill = NA))+
  theme(panel.grid.major = element_line(color='gray'),panel.grid.minor.x = element_line(color='gray'))+
  scale_fill_viridis_d(begin=0.2, end=0.8, option='inferno')+
  scale_color_viridis_d(begin=0.2, end=0.8, option='inferno')+
  theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))
rare_16sY1_cropped

rare_16sY2_cropped <- rare_apple16sY2 + facet_grid(Facility~.) + 
  theme(strip.background = element_rect(color='black', fill='white', size=0.8, linetype = 'solid'))+
  xlab("Number of reads (x 10,000)") + ylab("Number of unique OTUs") +
  scale_x_continuous(limits=c(0,50000), breaks= seq(0,50000, 5000) , labels = function(x){x/10000}) + ylim(0,2000)+
  theme(axis.text.x = element_text(size=8, color='black'), axis.text.y = element_text(size=8, color='black'),
        axis.text=element_text(size=8)) +
  theme(panel.background = element_rect(fill='transparent', color = NA),
        plot.background = element_rect(fill = 'transparent',color = NA), 
        panel.border = element_rect(color='black',fill = NA))+
  theme(panel.grid.major = element_line(color='gray'),panel.grid.minor.x = element_line(color='gray'))+
  scale_fill_viridis_d(begin=0.2, end=0.8, option='inferno')+
  scale_color_viridis_d(begin=0.2, end=0.8, option='inferno')+
  theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))
rare_16sY2_cropped

rare_ITSY1_cropped <- rare_appleITSY1 + facet_grid(Facility~.) + 
  theme(strip.background = element_rect(color='black', fill='white', size=0.8, linetype = 'solid'))+
  xlab("Number of reads (x 10,000)") + ylab("Number of unique OTUs") +
  scale_x_continuous(limits=c(0,50000), breaks= seq(0,50000, 5000) , labels = function(x){x/10000}) + ylim(0,1500)+
  theme(axis.text.x = element_text(size=8, color='black'), axis.text.y = element_text(size=8, color='black'),
        axis.text=element_text(size=8)) +
  theme(panel.background = element_rect(fill='transparent', color = NA),
        plot.background = element_rect(fill = 'transparent',color = NA), 
        panel.border = element_rect(color='black',fill = NA))+
  theme(panel.grid.major = element_line(color='gray'),panel.grid.minor.x = element_line(color='gray'))+
  scale_fill_viridis_d(begin=0.2, end=0.8, option='inferno')+
  scale_color_viridis_d(begin=0.2, end=0.8, option='inferno')+
  theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))
rare_ITSY1_cropped

rare_ITSY2_cropped <- rare_appleITSY2 + facet_grid(Facility~.) + 
  theme(strip.background = element_rect(color='black', fill='white', size=0.8, linetype = 'solid'))+
  xlab("Number of reads (x 10,000)") + ylab("Number of unique OTUs") +
  scale_x_continuous(limits=c(0,50000), breaks= seq(0,50000, 5000) , labels = function(x){x/10000}) + ylim(0,1500)+
  theme(axis.text.x = element_text(size=8, color='black'), axis.text.y = element_text(size=8, color='black'),
        axis.text=element_text(size=8)) +
  theme(panel.background = element_rect(fill='transparent', color = NA),
        plot.background = element_rect(fill = 'transparent',color = NA), 
        panel.border = element_rect(color='black',fill = NA))+
  theme(panel.grid.major = element_line(color='gray'),panel.grid.minor.x = element_line(color='gray'))+
  scale_fill_viridis_d(begin=0.2, end=0.8, option='inferno')+
  scale_color_viridis_d(begin=0.2, end=0.8, option='inferno')+
  theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))
rare_ITSY2_cropped

#Combine zoomed-in rarefaction plots
Rarefaction_cropped = plot_grid(rare_16sY1_cropped+ theme(legend.position = "none"),
                        rare_16sY2_cropped+ theme(legend.position = "none"),
                        rare_ITSY1_cropped+ theme(legend.position = "none"),
                        rare_ITSY2_cropped+ theme(legend.position = "none"),
                        ncol=2, nrow=2, labels = c("A","B","C","D"), label_size = 20)
Rarefaction_cropped
ggsave("Rarefaction_combined_cropped.png", plot =Rarefaction_cropped, device="png", width=8, height=10, units="in",dpi=600)



#### ESTIMATE % DISCOVERED DIVERSITY ####
#Code based on Taejung's work
##Calculate estimated richness for each sample.

## Original script adopted from "https://cran.r-project.org/web/packages/SpadeR/SpadeR.pdf"

## SpadeR::ChaoSpecies function is used to estimate richness.


richness_estimate = function(otu,...) {
  
  options(warn = -1)
  
  b = data.frame(matrix(nrow=as.matrix(dim(otu))[2,], ncol=3))
  
  colnames(b) <- c("Chao1 Estimates","Observed OTUs", "%Covered Species")
  
  for (i in 1:as.matrix(dim(otu))[2,]) {
    
    a =SpadeR::ChaoSpecies(otu[,i], datatype="abundance", k=10, conf=0.95)
    
    b[i,1]= as.numeric(a$Species_table[3,1])
    
    b[i,2]= apply(as.data.frame(a$Basic_data_information),2,as.numeric)[2,2]
    
    b[i,3]= (b[i,2]/b[i,1])*100
    
    rownames(b) <- colnames(otu) }
  
  print(b)
  
}


#Estimate richness using Chao1 index and % discovery
spadeR_16s_estimate <- richness_estimate(otus_16s)
spadeR_ITS_estimate <- richness_estimate(otus_ITS)


#Add year column
spadeR_16sY1_estimate$Year<-rep("Y1", 117)
spadeR_16sY2_estimate$Year<-rep("Y2", 107)
spadeR_ITSY1_estimate$Year<-rep("Y1", 117)
spadeR_ITSY2_estimate$Year<-rep("Y2", 107)

spadeR_16s_estimate$Year<-c(rep("Y1",117),rep("Y2",107))
spadeR_16s_estimate$Seq<-rep("16s",224)


#Merge with metadata
spadeR_16sY1<-bind_cols(spadeR_16sY1_estimate, metadata_16sY1)
spadeR_16sY2<-bind_cols(spadeR_16sY2_estimate, metadata_16sY2)
spadeR_ITSY1<-bind_cols(spadeR_ITSY1_estimate, metadata_ITSY1)
spadeR_ITSY2<-bind_cols(spadeR_ITSY2_estimate, metadata_ITSY2)

#Combine Y1 and Y2 data
spadeR_16sY1Y2<-bind_rows(spadeR_16sY1,spadeR_16sY2)
spadeR_ITSY1Y2<-bind_rows(spadeR_ITSY1,spadeR_ITSY2)

#Statistical analysis
#Summary statistics
describeBy(spadeR_16sY1Y2$`%Covered Species`, group=spadeR_16sY1Y2$Year, mat = TRUE) #By year
describeBy(spadeR_16sY1Y2$`%Covered Species`, list(spadeR_16sY1Y2$Year, spadeR_16sY1Y2$Facility), mat = TRUE) #By year and facility

describeBy(spadeR_ITSY1Y2$`%Covered Species`, group=spadeR_ITSY1Y2$Year, mat = TRUE) #By year
describeBy(spadeR_ITSY1Y2$`%Covered Species`, list(spadeR_ITSY1Y2$Year, spadeR_ITSY1Y2$Facility), mat = TRUE) #By year and facility

#t-test
t.test_16s<-t.test(`%Covered Species`~Year, data=spadeR_16sY1Y2) #t = -11.873, df = 136.94, p-value < 2.2e-16
t.test_ITS<-t.test(`%Covered Species`~Year, data=spadeR_ITSY1Y2) #t = 8.4661, df = 216.01, p-value = 3.888e-15



#### COMPOSITIONAL ANALYSIS OF MICROBIOME AT OTU LEVEL ####
#Based on Microbiome Analysis in R. Chap 10.

#Step 1: Convert OTU table to appropriate format
#Following step requires samples on rows and OTUs in columns
head(t(otus_16s)) #check that data in the correct format: samples on rows and OTUs in columns
head(t(otus_ITS))


#Step 2: Replace zero values before clr transformation. 
#Uses Count Zero Multiplicative method to replace zeros and outputs pseudo-counts
otu.n0_16s<-t(cmultRepl(t(otus_16s), label=0, method="CZM", output="p-counts")) #1740394 corrected values
otu.n0_ITS<-t(cmultRepl(t(otus_ITS), label=0, method="CZM", output="p-counts")) #989240  corrected values

head(otu.n0_16s) #output table needs to have samples in columns and OTUs in rows
head(otu.n0_ITS)

#Step 3: Convert data to proportions
otu.n0_16s_prop<-apply(otu.n0_16s, 2, function(x) {x/sum(x)})
otu.n0_ITS_prop<-apply(otu.n0_ITS, 2, function(x) {x/sum(x)})

#Step 4: Perform abundance and sample filtering and deal sparsity
#Filter the data to remove all taxa with less than 0.00001% abundance in any sample
otu.n0_16s_prop_f<-otu.n0_16s[apply(otu.n0_16s_prop, 1, min) > 0.0000001, ]
otu.n0_ITS_prop_f<-otu.n0_ITS[apply(otu.n0_ITS_prop, 1, min) > 0.0000001, ]

head(otu.n0_16s_prop_f) #Check that samples are on columns and OTUs in rows
head(otu.n0_ITS_prop_f)

#Step 5: perform CLR transformation
otu.n0.clr_16s<-t(apply(otu.n0_16s_prop_f, 2, function(x){log(x)-mean(log(x))}))
otu.n0.clr_ITS<-t(apply(otu.n0_ITS_prop_f, 2, function(x){log(x)-mean(log(x))}))

head(otu.n0.clr_16s) #Check output table. Samples should be in rows and OTUs in columns
head(otu.n0.clr_ITS)

#Step 6: Perform Singular Value Decomposition (PCA)
pc.clr_16s<-prcomp(otu.n0.clr_16s)
pc.clr_ITS<-prcomp(otu.n0.clr_ITS)

png("Screeplot - PCA.png", width = 400, height = 300, units = 'px')
par(mar=c(2,2,2,2))
par(mfrow=c(1,2))
screeplot(pc.clr_16s, type='lines', main="Bacteria")
screeplot(pc.clr_ITS, type='lines', main="Fungi ")
dev.off()

#Calculate total variance of the data
mvar.clr_16s<-mvar(otu.n0.clr_16s)
mvar.clr_ITS<-mvar(otu.n0.clr_ITS)

#Display results - 16s 
row_16s<-rownames(otu.n0.clr_16s) #Make vector with sample names
pc_out_16s<-as.data.frame(pc.clr_16s$x[,1:2]) #Get PC1 and PC2 
pc_out_meta_16s<-as.data.frame(bind_cols(pc_out_16s,metadata_16s)) #Add metadata information
row.names(pc_out_meta_16s)<-row_16s #Add rownames to dataframe
pc_out_meta_16s$Facility<-as.factor(pc_out_meta_16s$Facility)
pc_out_meta_16s$L..monocytogenes<-as.factor(pc_out_meta_16s$L..monocytogenes)
pc_out_meta_16s$FacY<-as.factor(pc_out_meta_16s$FacY)
pc_out_meta_16s$Year<-as.factor(pc_out_meta_16s$Year)
pc_out_meta_16s$Week<-as.factor(pc_out_meta_16s$Week)
pc_out_meta_16s$Month<-as.factor(pc_out_meta_16s$Month)

# Make PCA plot - First 2 axis- color by facility/Y and shape by L..mono
#Fig 2A
PCA_16s <- ggplot(pc_out_meta_16s, aes(x=PC1,y=PC2, color=Facility, shape=Year))+
  geom_point(size=3)+
  scale_shape_manual(values = c(15,12))+
  theme(legend.text=element_text(size=15,color='black'), legend.title= element_text(size=15, face='italic')) +
  theme(legend.position = 'right')+
  theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  theme(axis.title = element_text(size=15,color='black')) +
  theme(panel.background = element_rect(fill='transparent', color = NA),
        plot.background = element_rect(fill = 'transparent',color = NA), 
        panel.border = element_rect(color='black',fill = NA,size=1))+
  scale_x_continuous(name = paste("PC1: ", round(pc.clr_16s$sdev[1]^2/mvar.clr_16s*100, digits=1), "%", sep="")) +
  scale_y_continuous(name = paste("PC2: ", round(pc.clr_16s$sdev[2]^2/mvar.clr_16s*100, digits=1), "%", sep="")) +
  ggtitle("Bacteria", subtitle = "PCA by Facility and Season")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=16))+
  scale_color_viridis_d(begin = 0.2, end = 0.8, option='inferno')
PCA_16s
ggsave("PCA_Bacteria.svg", plot =PCA_16s, device="svg", width=6, height=5, units="in",dpi=600)

#Display results -ITS
row_ITS<-rownames(otu.n0.clr_ITS) #Make vector with sample names
pc_out_ITS<-as.data.frame(pc.clr_ITS$x[,1:2]) #Get PC1 and PC2
pc_out_meta_ITS<-as.data.frame(bind_cols(pc_out_ITS,metadata_ITS)) #Add metadata information
row.names(pc_out_meta_ITS)<-row_ITS #Add rownames to dataframe
pc_out_meta_ITS$Facility<-as.factor(pc_out_meta_ITS$Facility)
pc_out_meta_ITS$L..monocytogenes<-as.factor(pc_out_meta_ITS$L..monocytogenes)
pc_out_meta_ITS$FacY<-as.factor(pc_out_meta_ITS$FacY)
pc_out_meta_ITS$Year<-as.factor(pc_out_meta_ITS$Year)
pc_out_meta_ITS$Week<-as.factor(pc_out_meta_ITS$Week)
pc_out_meta_ITS$Month<-as.factor(pc_out_meta_ITS$Month)

# Make PCA plot - First 2 axis- color by facility/Y and shape by L..mono
#Fig 2D
PCA_ITS <- ggplot(pc_out_meta_ITS, aes(x=PC1,y=PC2, color=Facility, shape=Year))+
  geom_point(size=3)+
  scale_shape_manual(values = c(15,12))+
  theme(legend.text=element_text(size=15,color='black'), legend.title= element_text(size=13, face='italic')) +
  theme(legend.position = 'right')+
  theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  theme(axis.title = element_text(size=15,color='black')) +
  theme(panel.background = element_rect(fill='transparent', color = NA),
        plot.background = element_rect(fill = 'transparent',color = NA),
        panel.border = element_rect(color='black',fill = NA,size=1))+
  scale_x_continuous(name = paste("PC1: ", round(pc.clr_ITS$sdev[1]^2/mvar.clr_ITS*100, digits=1), "%", sep="")) +
  scale_y_continuous(name = paste("PC2: ", round(pc.clr_ITS$sdev[2]^2/mvar.clr_ITS*100, digits=1), "%", sep="")) +
    ggtitle("Fungi", subtitle = "PCA by Facility and Season")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=16))+
  scale_color_viridis_d(begin = 0.2, end = 0.8, option='inferno')
PCA_ITS
ggsave("PCA_Fungi.svg", plot =PCA_ITS, device="svg", width=6, height=5, units="in",dpi=600)


# PERMANOVA #
#Calculate Aitchinson distance
dist_16s<-dist(otu.n0.clr_16s, method='euclidean')
dist_ITS<-dist(otu.n0.clr_ITS, method='euclidean')

#16s-
permanova_16s<-pairwise.adonis2(dist_16s~Year:Facility+Facility+Year, data=metadata_16s, perm = 999, p.adjust.m = 'bonferroni')
permanova_16s
p.adjust(permanova_16s$`Y2_vs_Y1`$`Pr(>F)`, method = 'bonferroni')

#ITS
permanova_ITS<-pairwise.adonis2(dist_ITS~Year:Facility+Facility+Year, data=metadata_ITS, perm = 999, p.adjust.m = 'bonferroni')
permanova_ITS
p.adjust(permanova_ITS$`Y2_vs_Y1`$`Pr(>F)`, method = 'bonferroni')


##NOTE: Permanova showed significant interaction effect for Year:Facility. 
##Continue by (1) splitting OTU table by Year to see Facility effect within each year and
##(2) splitting OTU table by Facility to see Seasonal effect within each Facility

#### Composition of microbiota by year - Bubble plots by facility and year #### 
#OTU level plots

#Transform sample counts into compositions
otu.n0.acomp_16s<-as.data.frame(acomp(t(otu.n0_16s)), total=1)
otu.n0.acomp_ITS<-as.data.frame(acomp(t(otu.n0_ITS)), total=1)

#Make Phyloseq object for each year
OTU_16s <- otu_table(otu.n0.acomp_16s, taxa_are_rows = FALSE)
phyloseq_16s = phyloseq(OTU_16s, TAX_16s, META_16s)

OTU_ITS <- otu_table(otu.n0.acomp_ITS, taxa_are_rows = FALSE)
phyloseq_ITS = phyloseq(OTU_ITS, TAX_ITS, META_ITS)


#Make long format table from Phyloseq object
otu_16s_long <- phyloseq_16s %>%  
  transform_sample_counts(function(x) {x * 100} ) %>% #Recalculates the relative abundance 
  psmelt() %>%  #Melts to long format
  arrange(desc(Abundance))

otu_ITS_long <- phyloseq_ITS %>%  
  transform_sample_counts(function(x) {x * 100} ) %>% #Recalculates the relative abundance 
  psmelt() %>%  #Melts to long format
  arrange(desc(Abundance))


#Save OTU composition as .csv
write.csv(otu_16s_long, "otu_16s.csv")
write.csv(otu_ITS_long, "otu_ITS.csv")

#Make vector with OTU names above 5% relative abundance in at least one sample
otu_16s_over5<-unique(c(otu_16s_long$OTU[which(otu_16s_long$Abundance >=5)]))
otu_ITS_over5<-unique(c(otu_ITS_long$OTU[which(otu_ITS_long$Abundance >=5)])) 

#Filter table to obtain only OTUs with over 5% in at least one sample
otu_16s_over5abund <- filter(otu_16s_long, OTU %in% otu_16s_over5)
otu_ITS_over5abund <- filter(otu_ITS_long, OTU %in% otu_ITS_over5)

### Calculate mean relative abundance by Facility each OTU
otu_16s_over5abund_mean<-otu_16s_over5abund%>%
  group_by(OTU, Facility, Family, Year)%>%
  summarize(Mean=mean(Abundance))%>%
  arrange(desc(Mean))

otu_ITS_over5abund_mean<-otu_ITS_over5abund%>%
  group_by(OTU, Facility, Family, Year)%>%
  summarize(Mean=mean(Abundance))%>%
  arrange(desc(Mean))


#Remove tax rank labels from ITS Family labels
otu_ITS_over5abund_mean$Family_clean<-gsub("f__","",otu_ITS_over5abund_mean$Family)
otu_ITS_over5abund_mean$Family_clean<-gsub("o__","",otu_ITS_over5abund_mean$Family_clean)
otu_ITS_over5abund_mean$Family_clean<-gsub("p__","",otu_ITS_over5abund_mean$Family_clean)
otu_ITS_over5abund_mean$Family_clean<-gsub("c__","",otu_ITS_over5abund_mean$Family_clean)


#Bubbleplots by facility at the OTU level
#Fig S3-A
#16s
bubbleplot_otu16s<-ggplot(otu_16s_over5abund_mean, aes(x=Facility, y=reorder(OTU,desc(OTU)), size=Mean, color=Facility))+
  geom_point()+facet_grid(Family~Year, scales = "free", space = 'free')+
  scale_size(limits=c(0,50),name = "Relative abundance (%)")+
  theme(axis.title = element_blank(), axis.text.x=element_text(color='black', size=10), 
        axis.ticks=element_line(color='black'),
        axis.text.y = element_text(color='black', size=10)) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, nrow = 1)) +
  theme(legend.position = "bottom")+
  theme(panel.background = element_rect(fill="grey99", color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_rect(fill="white", color = "black"), 
        strip.text.y = element_text(size=10, angle=0, face = 'italic'), 
        panel.border = element_rect(color="black", fill=NA))+
  theme(legend.text=element_text(size=10,color='black'), legend.title= element_text(size=10)) +
  ggtitle("Bacteria")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
  scale_color_viridis_d(begin = 0.2, end = 0.8, option='inferno')+
  theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))
bubbleplot_otu16s
ggsave("Bubbleplot_otu16s.svg", plot=bubbleplot_otu16s, device="svg", width=8, height=10, units="in", dpi=600)

#ITS
#Fig S3-B
bubbleplot_otuITS<-ggplot(otu_ITS_over5abund_mean, aes(x=Facility, y=reorder(OTU,desc(OTU)), size=Mean, color=Facility))+
  geom_point()+facet_grid(Family_clean~Year, scales = "free", space = 'free')+
  scale_size(limits=c(0,50), name = "Relative abundance (%)")+
  theme(axis.title = element_blank(), axis.text.x=element_text(color='black', size=16), 
        axis.ticks=element_line(color='black'),
        axis.text.y = element_text(color='black', size=16)) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  theme(panel.background = element_rect(fill="grey99", color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_rect(fill="white", color = "black"), 
        strip.text.y = element_text(size=16, angle=0, face = 'italic'), 
        panel.border = element_rect(color="black", fill=NA))+
  theme(legend.text=element_text(size=15,color='black'), legend.title= element_text(size=16)) +
  ggtitle("Bacteria")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=18))+
  scale_color_viridis_d(begin = 0.2, end = 0.8, option='inferno')+
  theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))
bubbleplot_otuITS
ggsave("Bubbleplot_otuITS.svg", plot=bubbleplot_otuITS, device="svg", width=8, height=10, units="in", dpi=600)


#### Heatmap of microbiota by presence of Lm at OTU level by year ####


### Calculate mean relative abundance by Facility each OTU
otu_16s_over5abund_mean_Lm<-otu_16s_over5abund%>%
  group_by(OTU, Facility, Family, Year, L..monocytogenes)%>%
  summarize(Mean=mean(Abundance))%>%
  arrange(desc(Mean))

otu_ITS_over5abund_mean_Lm<-otu_ITS_over5abund%>%
  group_by(OTU, Facility, Family, Year, L..monocytogenes)%>%
  summarize(Mean=mean(Abundance))%>%
  arrange(desc(Mean))


#Remove tax rank labels from ITS Family labels
otu_ITS_over5abund_mean_Lm$Family_clean<-gsub("f__","",otu_ITS_over5abund_mean_Lm$Family)
otu_ITS_over5abund_mean_Lm$Family_clean<-gsub("o__","",otu_ITS_over5abund_mean_Lm$Family_clean)
otu_ITS_over5abund_mean_Lm$Family_clean<-gsub("p__","",otu_ITS_over5abund_mean_Lm$Family_clean)
otu_ITS_over5abund_mean_Lm$Family_clean<-gsub("c__","",otu_ITS_over5abund_mean_Lm$Family_clean)


#Heatmaps by facility at the OTU level
#Fig 5A
#16s
heatmap_otu16s<-ggplot(otu_16s_over5abund_mean_Lm, aes(x=L..monocytogenes, y=reorder(OTU,desc(OTU)), fill=Mean))+
  geom_tile(color='black')+facet_grid(Family~Facility+Year, scales = "free", space = 'free')+
  scale_size(limits=c(0,50),name = "Relative abundance (%)")+
  theme(axis.title = element_blank(), axis.text.x=element_text(color='black', size=10), 
        axis.ticks=element_line(color='black'),
        axis.text.y = element_text(color='black', size=10)) + 
  theme(legend.position = "bottom")+
  theme(panel.background = element_rect(fill="grey99", color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_rect(fill="white", color = "black"), 
        strip.text.y = element_text(size=10, angle=0, face = 'italic'), 
        panel.border = element_rect(color="black", fill=NA))+
  theme(legend.text=element_text(size=10,color='black'), legend.title= element_text(size=10)) +
  ggtitle("Bacteria")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
  scale_fill_viridis_c(begin = 1, end = 0, option='inferno')+
  theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))
heatmap_otu16s
ggsave("Heatmap_otu16s.svg", plot=heatmap_otu16s, device="svg", width=8, height=10, units="in", dpi=600)
ggsave("Heatmap_otu16s.png", plot=heatmap_otu16s, device="png", width=8, height=10, units="in", dpi=600)

#ITS
#Fig 6A
heatmap_otuITS<-ggplot(otu_ITS_over5abund_mean_Lm, aes(x=L..monocytogenes, y=reorder(OTU,desc(OTU)), fill=Mean))+
  geom_tile(color='black')+facet_grid(Family~Facility+Year, scales = "free", space = 'free')+
  scale_size(limits=c(0,50),name = "Relative abundance (%)")+
  theme(axis.title = element_blank(), axis.text.x=element_text(color='black', size=10), 
        axis.ticks=element_line(color='black'),
        axis.text.y = element_text(color='black', size=10)) + 
  theme(legend.position = "bottom")+
  theme(panel.background = element_rect(fill="grey99", color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_rect(fill="white", color = "black"), 
        strip.text.y = element_text(size=10, angle=0, face = 'italic'), 
        panel.border = element_rect(color="black", fill=NA))+
  theme(legend.text=element_text(size=10,color='black'), legend.title= element_text(size=10)) +
  ggtitle("Fungi")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=18))+
  scale_fill_viridis_c(begin = 1, end = 0, option='inferno')+
  theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))
heatmap_otuITS
ggsave("heatmapt_otuITS.svg", plot=heatmap_otuITS, device="svg", width=8, height=10, units="in", dpi=600)
ggsave("heatmapt_otuITS.png", plot=heatmap_otuITS, device="png", width=8, height=10, units="in", dpi=600)



#### Compositional analysis by year at the OTU level ####

#Make phyloseq
OTU_16s <- otu_table(otus_16s, taxa_are_rows = TRUE)
phyloseq_16s = phyloseq(OTU_16s, TAX_16s, META_16s)
TREE_16s = rtree(ntaxa(phyloseq_16s), rooted=TRUE, tip.label = taxa_names(phyloseq_16s))  
phyloseq16s <- phyloseq(OTU_16s, TAX_16s, TREE_16s, META_16s)

OTU_ITS <- otu_table(otus_ITS, taxa_are_rows = TRUE)
phyloseq_ITS = phyloseq(OTU_ITS, TAX_ITS, META_ITS)
TREE_ITS = rtree(ntaxa(phyloseq_ITS), rooted=TRUE, tip.label = taxa_names(phyloseq_ITS))
phyloseqITS <- phyloseq(OTU_ITS, TAX_ITS, TREE_ITS, META_ITS)
  
  
# Subset Phyloseq for each year
physeq_16sY1 <- subset_samples(phyloseq16s, Year == "Y1") 
physeq_ITSY1 <- subset_samples(phyloseqITS, Year == "Y1")

physeq_16sY2 <- subset_samples(phyloseq16s, Year == "Y2") 
physeq_ITSY2 <- subset_samples(phyloseqITS, Year == "Y2")

#Export OTU table for each facility from Phyloseq 
otus_16sY1<-as.data.frame(as(otu_table(physeq_16sY1), "matrix"))
otus_ITSY1<-as.data.frame(as(otu_table(physeq_ITSY1), "matrix"))

otus_16sY2<-as.data.frame(as(otu_table(physeq_16sY2), "matrix"))
otus_ITSY2<-as.data.frame(as(otu_table(physeq_ITSY2), "matrix"))

#Subset metadata
metadata_16sY1<-subset(metadata_16s, Year=="Y1")
metadata_16sY2<-subset(metadata_16s, Year=="Y2")

metadata_ITSY1<-subset(metadata_ITS, Year=="Y1")
metadata_ITSY2<-subset(metadata_ITS, Year=="Y2")

#Remove OTUs that have count zero in all samples - Necessary step for the zero imputation function
otus_16sY1<-otus_16sY1[ which(rowSums(otus_16sY1)>0),]
otus_16sY2<-otus_16sY2[ which(rowSums(otus_16sY2)>0),]

otus_ITSY1<-otus_ITSY1[ which(rowSums(otus_ITSY1)>0),]
otus_ITSY2<-otus_ITSY2[ which(rowSums(otus_ITSY2)>0),]

#Step 1: Convert OTU table to appropriate format
#Following step requires samples on rows and OTUs in columns
head(t(otus_16sY1)) 
head(t(otus_ITSY1))

head(t(otus_16sY2)) 
head(t(otus_ITSY2))

#Step 2: Replace zero values with a small non- zero value before clr transformation. 
#Uses Count Zero Multiplicative method to replace zeros and outputs pseudo-counts
otu.n0_16sY1<-t(cmultRepl(t(otus_16sY1), label=0, method="CZM", output="p-counts")) #349844  corrected values
otu.n0_ITSY1<-t(cmultRepl(t(otus_ITSY1), label=0, method="CZM", output="p-counts")) #184651  corrected values

otu.n0_16sY2<-t(cmultRepl(t(otus_16sY2), label=0, method="CZM", output="p-counts")) #633448   corrected values
otu.n0_ITSY2<-t(cmultRepl(t(otus_ITSY2), label=0, method="CZM", output="p-counts")) #159654   corrected values

#output table needs to have samples in columns and OTUs in rows
head(otu.n0_16sY1) 
head(otu.n0_ITSY1)

head(otu.n0_16sY2) 
head(otu.n0_ITSY2)

#Step 3: Convert data to proportions
otu.n0_16sY1_prop<-apply(otu.n0_16sY1, 2, function(x) {x/sum(x)})
otu.n0_ITSY1_prop<-apply(otu.n0_ITSY1, 2, function(x) {x/sum(x)})

otu.n0_16sY2_prop<-apply(otu.n0_16sY2, 2, function(x) {x/sum(x)})
otu.n0_ITSY2_prop<-apply(otu.n0_ITSY2, 2, function(x) {x/sum(x)})

head(otu.n0_16sY1_prop)
head(otu.n0_ITSY1_prop)

head(otu.n0_16sY2_prop)
head(otu.n0_ITSY2_prop)

#Step 4: Perform abundance and sample filtering and deal sparsity
#Filter the data to remove all taxa with less than 0.00001% abundance in any sample
otu.n0_16sY1_prop_f<-otu.n0_16sY1[apply(otu.n0_16sY1_prop, 1, min) > 0.0000001, ]
otu.n0_ITSY1_prop_f<-otu.n0_ITSY1[apply(otu.n0_ITSY1_prop, 1, min) > 0.0000001, ]

otu.n0_16sY2_prop_f<-otu.n0_16sY2[apply(otu.n0_16sY2_prop, 1, min) > 0.0000001, ]
otu.n0_ITSY2_prop_f<-otu.n0_ITSY2[apply(otu.n0_ITSY2_prop, 1, min) > 0.0000001, ]

#Check that samples are on columns and OTUs in rows
head(otu.n0_16sY1_prop_f) 
head(otu.n0_ITSY1_prop_f)

head(otu.n0_16sY2_prop_f)
head(otu.n0_ITSY2_prop_f)

#Step 5: perform CLR transformation
otu.n0.clr_16sY1<-t(apply(otu.n0_16sY1_prop_f, 2, function(x){log(x)-mean(log(x))}))
otu.n0.clr_ITSY1<-t(apply(otu.n0_ITSY1_prop_f, 2, function(x){log(x)-mean(log(x))}))

otu.n0.clr_16sY2<-t(apply(otu.n0_16sY2_prop_f, 2, function(x){log(x)-mean(log(x))}))
otu.n0.clr_ITSY2<-t(apply(otu.n0_ITSY2_prop_f, 2, function(x){log(x)-mean(log(x))}))

#Check output table. Samples should be in rows and OTUs in columns
head(otu.n0.clr_16sY1) 
head(otu.n0.clr_ITSY1)

head(otu.n0.clr_16sY2) 
head(otu.n0.clr_ITSY2)

#Step 6: Perform Singular Value Decomposition (PCA)
pc.clr_16sY1<-prcomp(otu.n0.clr_16sY1)
pc.clr_ITSY1<-prcomp(otu.n0.clr_ITSY1)

pc.clr_16sY2<-prcomp(otu.n0.clr_16sY2)
pc.clr_ITSY2<-prcomp(otu.n0.clr_ITSY2)

png("Screeplot - PCA by year .png", width = 1000, height = 500, units = 'px')
par(mar=c(2,2,2,2))
par(mfrow=c(2,2))
screeplot(pc.clr_16sY1, type='lines', main="Bacteria Y1")
screeplot(pc.clr_16sY2, type='lines', main="Bacteria Y2")
screeplot(pc.clr_ITSY1, type='lines', main="Fungi Y1")
screeplot(pc.clr_ITSY2, type='lines', main="Fungi Y2")
dev.off()

#Calculate total variance of the data
mvar.clr_16sY1<-mvar(otu.n0.clr_16sY1)
mvar.clr_ITSY1<-mvar(otu.n0.clr_ITSY1)

mvar.clr_16sY2<-mvar(otu.n0.clr_16sY2)
mvar.clr_ITSY2<-mvar(otu.n0.clr_ITSY2)

#Display results - 16sY1
row_16sY1<-rownames(otu.n0.clr_16sY1) #Make vector with sample names
pc_out_16sY1<-as.data.frame(pc.clr_16sY1$x[,1:2]) #Get PC1 and PC2 
pc_out_meta_16sY1<-as.data.frame(bind_cols(pc_out_16sY1,metadata_16sY1)) #Add metadata information
row.names(pc_out_meta_16sY1)<-row_16sY1 #Add rownames to dataframe
pc_out_meta_16sY1$Facility<-as.factor(pc_out_meta_16sY1$Facility)
pc_out_meta_16sY1$L..monocytogenes<-as.factor(pc_out_meta_16sY1$L..monocytogenes)

# Make PCA plot - First 2 axis- color by facility and shape by L..mono
#Fig 2B
PCA_16sY1<- ggplot(pc_out_meta_16sY1, aes(x=PC1,y=PC2, color=Facility, shape=L..monocytogenes))+
  geom_point(size=3)+
  theme(legend.text=element_text(size=15,color='black'), legend.title= element_text(size=15, face='italic')) +
  theme(legend.position = 'right')+
  theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  theme(axis.title = element_text(size=15,color='black')) +
  theme(panel.background = element_rect(fill='transparent', color = NA),
        plot.background = element_rect(fill = 'transparent',color = NA), 
        panel.border = element_rect(color='black',fill = NA,size=1))+
  scale_x_continuous(name = paste("PC1: ", round(pc.clr_16sY1$sdev[1]^2/mvar.clr_16sY1*100, digits=1), "%", sep="")) +
  scale_y_continuous(name = paste("PC2: ", round(pc.clr_16sY1$sdev[2]^2/mvar.clr_16sY1*100, digits=1), "%", sep="")) +
  ggtitle("Bacteria - Y1", subtitle = "PCA by Facility and Lm presence")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=16))+
  scale_color_viridis_d(begin = 0.2, end = 0.8, option='inferno')
PCA_16sY1
ggsave("PCA_16sY1.svg", plot =PCA_16sY1, device="svg", width=6, height=5, units="in",dpi=600)

#Display results - 16sY2 
row_16sY2<-rownames(otu.n0.clr_16sY2) #Make vector with sample names
pc_out_16sY2<-as.data.frame(pc.clr_16sY2$x[,1:2]) #Get PC1 and PC2 
pc_out_meta_16sY2<-as.data.frame(bind_cols(pc_out_16sY2,metadata_16sY2)) #Add metadata information
row.names(pc_out_meta_16sY2)<-row_16sY2 #Add rownames to dataframe
pc_out_meta_16sY2$Facility<-as.factor(pc_out_meta_16sY2$Facility)
pc_out_meta_16sY2$L..monocytogenes<-as.factor(pc_out_meta_16sY2$L..monocytogenes)

# Make PCA plot - First 2 axis- color by facility and shape by L..mono
#Fig 2C
PCA_16sY2 <- ggplot(pc_out_meta_16sY2, aes(x=PC1,y=PC2, color=Facility, shape=L..monocytogenes))+
  geom_point(size=3)+
  theme(legend.text=element_text(size=15,color='black'), legend.title= element_text(size=15, face='italic')) +
  theme(legend.position = 'right')+
  theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  theme(axis.title = element_text(size=15,color='black')) +
  theme(panel.background = element_rect(fill='transparent', color = NA),
        plot.background = element_rect(fill = 'transparent',color = NA), 
        panel.border = element_rect(color='black',fill = NA,size=1))+
  scale_x_continuous(name = paste("PC1: ", round(pc.clr_16sY2$sdev[1]^2/mvar.clr_16sY2*100, digits=1), "%", sep="")) +
  scale_y_continuous(name = paste("PC2: ", round(pc.clr_16sY2$sdev[2]^2/mvar.clr_16sY2*100, digits=1), "%", sep="")) +
  ggtitle("Bacteria - Y2", subtitle = "PCA by Facility and Lm presence")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=16))+
  scale_color_viridis_d(begin = 0.2, end = 0.8, option='inferno')
PCA_16sY2
ggsave("PCA_16sY2.svg", plot =PCA_16sY2, device="svg", width=6, height=5, units="in",dpi=600)

#Display results - ITSY1
row_ITSY1<-rownames(otu.n0.clr_ITSY1) #Make vector with sample names
pc_out_ITSY1<-as.data.frame(pc.clr_ITSY1$x[,1:2]) #Get PC1 and PC2 
pc_out_meta_ITSY1<-as.data.frame(bind_cols(pc_out_ITSY1,metadata_ITSY1)) #Add metadata information
row.names(pc_out_meta_ITSY1)<-row_ITSY1 #Add rownames to dataframe
pc_out_meta_ITSY1$Facility<-as.factor(pc_out_meta_ITSY1$Facility)
pc_out_meta_ITSY1$L..monocytogenes<-as.factor(pc_out_meta_ITSY1$L..monocytogenes)

# Make PCA plot - First 2 axis- color by facility and shape by L..mono
#Fig 2E
PCA_ITSY1 <- ggplot(pc_out_meta_ITSY1, aes(x=PC1,y=PC2, color=Facility, shape=L..monocytogenes))+
  geom_point(size=3)+
  theme(legend.text=element_text(size=15,color='black'), legend.title= element_text(size=15, face='italic')) +
  theme(legend.position = 'right')+
  theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  theme(axis.title = element_text(size=15,color='black')) +
  theme(panel.background = element_rect(fill='transparent', color = NA),
        plot.background = element_rect(fill = 'transparent',color = NA), 
        panel.border = element_rect(color='black',fill = NA,size=1))+
  scale_x_continuous(name = paste("PC1: ", round(pc.clr_ITSY1$sdev[1]^2/mvar.clr_ITSY1*100, digits=1), "%", sep="")) +
  scale_y_continuous(name = paste("PC2: ", round(pc.clr_ITSY1$sdev[2]^2/mvar.clr_ITSY1*100, digits=1), "%", sep="")) +
  ggtitle("Fungi - Y1", subtitle = "PCA by Facility and Lm presence")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=16))+
  scale_color_viridis_d(begin = 0.2, end = 0.8, option='inferno')
PCA_ITSY1
ggsave("PCA_ITSY1.svg", plot =PCA_ITSY1, device="svg", width=6, height=5, units="in",dpi=600)

#Display results - ITSY2
row_ITSY2<-rownames(otu.n0.clr_ITSY2) #Make vector with sample names
pc_out_ITSY2<-as.data.frame(pc.clr_ITSY2$x[,1:2]) #Get PC1 and PC2 
pc_out_meta_ITSY2<-as.data.frame(bind_cols(pc_out_ITSY2,metadata_ITSY2)) #Add metadata information
row.names(pc_out_meta_ITSY2)<-row_ITSY2 #Add rownames to dataframe
pc_out_meta_ITSY2$Facility<-as.factor(pc_out_meta_ITSY2$Facility)
pc_out_meta_ITSY2$L..monocytogenes<-as.factor(pc_out_meta_ITSY2$L..monocytogenes)

# Make PCA plot - First 2 axis- color by facility and shape by L..mono
#Fig 2E
PCA_ITSY2 <- ggplot(pc_out_meta_ITSY2, aes(x=PC1,y=PC2, color=Facility, shape=L..monocytogenes))+
  geom_point(size=3)+
  theme(legend.text=element_text(size=15,color='black'), legend.title= element_text(size=15, face='italic')) +
  theme(legend.position = 'right')+
  theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  theme(axis.title = element_text(size=15,color='black')) +
  theme(panel.background = element_rect(fill='transparent', color = NA),
        plot.background = element_rect(fill = 'transparent',color = NA), 
        panel.border = element_rect(color='black',fill = NA,size=1))+
  scale_x_continuous(name = paste("PC1: ", round(pc.clr_ITSY2$sdev[1]^2/mvar.clr_ITSY2*100, digits=1), "%", sep="")) +
  scale_y_continuous(name = paste("PC2: ", round(pc.clr_ITSY2$sdev[2]^2/mvar.clr_ITSY2*100, digits=1), "%", sep="")) +
  ggtitle("Fungi - Y2", subtitle = "PCA by Facility and Lm presence")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=16))+
  scale_color_viridis_d(begin = 0.2, end = 0.8, option='inferno')
PCA_ITSY2
ggsave("PCA_ITSY2.svg", plot =PCA_ITSY2, device="svg", width=6, height=5, units="in",dpi=600)


# PERMANOVA #
#Using Aitchinson distances calculated previously for the dendrogram

#Two-way anova by facility and Lmono presence
#16s

permanova_16sY1<-pairwise.adonis2(dist_16sY1~L..monocytogenes:Facility+Facility+L..monocytogenes, data=metadata_16sY1, perm = 999, p.adjust.m = 'bonferroni')
permanova_16sY1
p.adjust(permanova_16sY1$`+_vs_-`$`Pr(>F)`, method = 'bonferroni')

permanova_16sY2<-pairwise.adonis2(dist_16sY2~L..monocytogenes:Facility+Facility+L..monocytogenes, data=metadata_16sY2, perm = 999, p.adjust.m = 'bonferroni')
permanova_16sY2
p.adjust(permanova_16sY2$`+_vs_-`$`Pr(>F)`, method = 'bonferroni')


#ITS
permanova_ITSY1<-pairwise.adonis2(dist_ITSY1~L..monocytogenes:Facility+Facility+L..monocytogenes, data=metadata_ITSY1, perm = 999, p.adjust.m = 'bonferroni')
permanova_ITSY1
p.adjust(permanova_ITSY1$`+_vs_-`$`Pr(>F)`, method = 'bonferroni')

permanova_ITSY2<-pairwise.adonis2(dist_ITSY2~L..monocytogenes:Facility+Facility+L..monocytogenes, data=metadata_ITSY2, perm = 999, p.adjust.m = 'bonferroni')
permanova_ITSY2
p.adjust(permanova_ITSY2$`+_vs_-`$`Pr(>F)`, method = 'bonferroni')

#NOTE: For the 4 datasets, there was a significant difference by facilities. Run one-way anova to see where the differences are.

#Pairwize one-way anova by facility
pairwise.adonis(dist_16sY1, factors=metadata_16sY1$Facility, perm = 999, p.adjust.m = 'bonferroni')
pairwise.adonis(dist_16sY2, factors=metadata_16sY2$Facility, perm = 999, p.adjust.m = 'bonferroni')

pairwise.adonis(dist_ITSY1, factors=metadata_ITSY1$Facility, perm = 999, p.adjust.m = 'bonferroni')
pairwise.adonis(dist_ITSY2, factors=metadata_ITSY2$Facility, perm = 999, p.adjust.m = 'bonferroni')



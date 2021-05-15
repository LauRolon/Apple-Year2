#Two-year monitoring of Lm in apple packing houses
#Analysis of microbiomes and mycobiomes using CoDa approach
#Laura Rolon
#Last updated: 05/04/21

#Load packages
library(ape)
library(vegan)
library(ggplot2)
library(phyloseq)
library(cowplot)
library(tidyr)
library(dplyr)
library(ggpubr)
library(compositions)
library(zCompositions)
library(ALDEx2)
library(viridis)
library(dendextend)
library(readxl)
library(gplots)
library(BiocParallel)
library(pairwiseAdonis)
library(SpadeR)
library(psych)
library(plotly)
library(htmlwidgets)


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

#Import taxonomy table -ITSY1
taxon_ITS <- as.data.frame(import_mothur(mothur_constaxonomy_file = 'appleITS.cons.taxonomy'))
colnames(taxon_ITS) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus")
TAX_ITS =tax_table(as.matrix(taxon_ITS))

#Import metadata -ITS
metadata_ITS <-read.csv("metadata_apple_ITS.csv", header=TRUE, row.names=1)
META_ITS = sample_data(metadata_ITS)



#### RAREFACTION CURVES ####
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
spadeR_16sv138_estimate <- richness_estimate(otus_16sv138)
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
# 
# #Violin plots for alpha diversity (Chao1 index)
# alpha_shannon_ITSY1 <- ggviolin(spadeR_16sY1Y2, x = "Year", y = "Chao1 Estimates", add = "boxplot", 
#                                 fill= "Year" ) +
#   theme(axis.text.x = element_text(size=25), axis.text.y = element_text(size=25)) +
#   theme(axis.title = element_blank())+
#   theme(plot.margin=margin(t=1, b=0.5, l=0.5, r=0.5, unit = 'cm'))+
#   ylim(0,10)+
#   scale_fill_manual(values=c("#EF8A62","#D0747F","#67A9CF"))+
#   theme(panel.background = element_rect(fill='transparent', color = NA),
#         plot.background = element_rect(fill = 'transparent',color = NA), 
#         panel.border = element_rect(color='black',fill = NA,size=1))
# alpha_shannon_ITSY1


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
PCA_16s <- ggplot(pc_out_meta_16s, aes(x=PC1,y=PC2, color=Facility, shape=Year))+
  geom_point(size=3)+
  theme(legend.text=element_text(size=13,color='black'), legend.title= element_text(size=13, face='italic')) +
  theme(legend.position = 'right')+
  theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  theme(axis.title = element_text(size=13,color='black')) +
  theme(panel.background = element_rect(fill='transparent', color = NA),
        plot.background = element_rect(fill = 'transparent',color = NA), 
        panel.border = element_rect(color='black',fill = NA,size=1))+
   scale_x_continuous(name = paste("PC1: ", round(pc.clr_16s$sdev[1]^2/mvar.clr_16s*100, digits=1), "%", sep="")) +
  scale_y_continuous(name = paste("PC2: ", round(pc.clr_16s$sdev[2]^2/mvar.clr_16s*100, digits=1), "%", sep="")) +
  ggtitle("Bacteria - ", subtitle = "PCA by Facility and Season")+theme(plot.title = element_text(hjust = 0.5, face = 'bold'))+
  scale_color_viridis_d(begin = 0.2, end = 0.8, option='inferno')
PCA_16s


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
PCA_ITS <- ggplot(pc_out_meta_ITS, aes(x=PC1,y=PC2, color=Facility, shape=Year))+
  geom_point(size=3)+
  theme(legend.text=element_text(size=13,color='black'), legend.title= element_text(size=13, face='italic')) +
  theme(legend.position = 'right')+
  theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  theme(axis.title = element_text(size=13,color='black')) +
  theme(panel.background = element_rect(fill='transparent', color = NA),
        plot.background = element_rect(fill = 'transparent',color = NA),
        panel.border = element_rect(color='black',fill = NA,size=1))+
  ggtitle("Fungi", subtitle = "PCA by Facility and Season")+theme(plot.title = element_text(hjust = 0.5, face = 'bold'))+
  scale_color_viridis_d(begin = 0.2, end = 0.8, option='inferno')
PCA_ITS

#Combine bubble plots
PCA_all = plot_grid(PCA_16s,PCA_ITS, 
                        ncol=2, nrow=1)
PCA_all


ggsave("PCA_all.png", plot =PCA_all, device="png", width=10, height=5, units="in",dpi=600)


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
PCA_16sY1<- ggplot(pc_out_meta_16sY1, aes(x=PC1,y=PC2, color=Facility, shape=L..monocytogenes))+
  geom_point(size=3)+
  theme(legend.text=element_text(size=13,color='black'), legend.title= element_text(size=13, face='italic')) +
  theme(legend.position = 'right')+
  theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  theme(axis.title = element_text(size=13,color='black')) +
  theme(panel.background = element_rect(fill='transparent', color = NA),
        plot.background = element_rect(fill = 'transparent',color = NA), 
        panel.border = element_rect(color='black',fill = NA,size=1))+
  scale_x_continuous(name = paste("PC1: ", round(pc.clr_16sY1$sdev[1]^2/mvar.clr_16sY1*100, digits=1), "%", sep="")) +
  scale_y_continuous(name = paste("PC2: ", round(pc.clr_16sY1$sdev[2]^2/mvar.clr_16sY1*100, digits=1), "%", sep="")) +
  ggtitle("Bacteria - Y1", subtitle = "PCA by Facility and Lm presence")+theme(plot.title = element_text(hjust = 0.5, face = 'bold'))+
  scale_color_viridis_d(begin = 0.2, end = 0.8, option='inferno')
PCA_16sY1

#Display results - 16sY2 
row_16sY2<-rownames(otu.n0.clr_16sY2) #Make vector with sample names
pc_out_16sY2<-as.data.frame(pc.clr_16sY2$x[,1:2]) #Get PC1 and PC2 
pc_out_meta_16sY2<-as.data.frame(bind_cols(pc_out_16sY2,metadata_16sY2)) #Add metadata information
row.names(pc_out_meta_16sY2)<-row_16sY2 #Add rownames to dataframe
pc_out_meta_16sY2$Facility<-as.factor(pc_out_meta_16sY2$Facility)
pc_out_meta_16sY2$L..monocytogenes<-as.factor(pc_out_meta_16sY2$L..monocytogenes)

# Make PCA plot - First 2 axis- color by facility and shape by L..mono
PCA_16sY2 <- ggplot(pc_out_meta_16sY2, aes(x=PC1,y=PC2, color=Facility, shape=L..monocytogenes))+
  geom_point(size=3)+
  theme(legend.text=element_text(size=13,color='black'), legend.title= element_text(size=13, face='italic')) +
  theme(legend.position = 'right')+
  theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  theme(axis.title = element_text(size=13,color='black')) +
  theme(panel.background = element_rect(fill='transparent', color = NA),
        plot.background = element_rect(fill = 'transparent',color = NA), 
        panel.border = element_rect(color='black',fill = NA,size=1))+
  scale_x_continuous(name = paste("PC1: ", round(pc.clr_16sY2$sdev[1]^2/mvar.clr_16sY2*100, digits=1), "%", sep="")) +
  scale_y_continuous(name = paste("PC2: ", round(pc.clr_16sY2$sdev[2]^2/mvar.clr_16sY2*100, digits=1), "%", sep="")) +
  ggtitle("Bacteria - Y2", subtitle = "PCA by Facility and Lm presence")+theme(plot.title = element_text(hjust = 0.5, face = 'bold'))+
  scale_color_viridis_d(begin = 0.2, end = 0.8, option='inferno')
PCA_16sY2

#Display results - ITSY1
row_ITSY1<-rownames(otu.n0.clr_ITSY1) #Make vector with sample names
pc_out_ITSY1<-as.data.frame(pc.clr_ITSY1$x[,1:2]) #Get PC1 and PC2 
pc_out_meta_ITSY1<-as.data.frame(bind_cols(pc_out_ITSY1,metadata_ITSY1)) #Add metadata information
row.names(pc_out_meta_ITSY1)<-row_ITSY1 #Add rownames to dataframe
pc_out_meta_ITSY1$Facility<-as.factor(pc_out_meta_ITSY1$Facility)
pc_out_meta_ITSY1$L..monocytogenes<-as.factor(pc_out_meta_ITSY1$L..monocytogenes)

# Make PCA plot - First 2 axis- color by facility and shape by L..mono
PCA_ITSY1 <- ggplot(pc_out_meta_ITSY1, aes(x=PC1,y=PC2, color=Facility, shape=L..monocytogenes))+
  geom_point(size=3)+
  theme(legend.text=element_text(size=13,color='black'), legend.title= element_text(size=13, face='italic')) +
  theme(legend.position = 'right')+
  theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  theme(axis.title = element_text(size=13,color='black')) +
  theme(panel.background = element_rect(fill='transparent', color = NA),
        plot.background = element_rect(fill = 'transparent',color = NA), 
        panel.border = element_rect(color='black',fill = NA,size=1))+
  scale_x_continuous(name = paste("PC1: ", round(pc.clr_ITSY1$sdev[1]^2/mvar.clr_ITSY1*100, digits=1), "%", sep="")) +
  scale_y_continuous(name = paste("PC2: ", round(pc.clr_ITSY1$sdev[2]^2/mvar.clr_ITSY1*100, digits=1), "%", sep="")) +
  ggtitle("Fungi - Y1", subtitle = "PCA by Facility and Lm presence")+theme(plot.title = element_text(hjust = 0.5, face = 'bold'))+
  scale_color_viridis_d(begin = 0.2, end = 0.8, option='inferno')
PCA_ITSY1

#Display results - ITSY2
row_ITSY2<-rownames(otu.n0.clr_ITSY2) #Make vector with sample names
pc_out_ITSY2<-as.data.frame(pc.clr_ITSY2$x[,1:2]) #Get PC1 and PC2 
pc_out_meta_ITSY2<-as.data.frame(bind_cols(pc_out_ITSY2,metadata_ITSY2)) #Add metadata information
row.names(pc_out_meta_ITSY2)<-row_ITSY2 #Add rownames to dataframe
pc_out_meta_ITSY2$Facility<-as.factor(pc_out_meta_ITSY2$Facility)
pc_out_meta_ITSY2$L..monocytogenes<-as.factor(pc_out_meta_ITSY2$L..monocytogenes)

# Make PCA plot - First 2 axis- color by facility and shape by L..mono
PCA_ITSY2 <- ggplot(pc_out_meta_ITSY2, aes(x=PC1,y=PC2, color=Facility, shape=L..monocytogenes))+
  geom_point(size=3)+
  theme(legend.text=element_text(size=13,color='black'), legend.title= element_text(size=13, face='italic')) +
  theme(legend.position = 'right')+
  theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  theme(axis.title = element_text(size=13,color='black')) +
  theme(panel.background = element_rect(fill='transparent', color = NA),
        plot.background = element_rect(fill = 'transparent',color = NA), 
        panel.border = element_rect(color='black',fill = NA,size=1))+
  scale_x_continuous(name = paste("PC1: ", round(pc.clr_ITSY2$sdev[1]^2/mvar.clr_ITSY2*100, digits=1), "%", sep="")) +
  scale_y_continuous(name = paste("PC2: ", round(pc.clr_ITSY2$sdev[2]^2/mvar.clr_ITSY2*100, digits=1), "%", sep="")) +
  ggtitle("Fungi - Y2", subtitle = "PCA by Facility and Lm presence")+theme(plot.title = element_text(hjust = 0.5, face = 'bold'))+
  scale_color_viridis_d(begin = 0.2, end = 0.8, option='inferno')
PCA_ITSY2

PCA_byyear = plot_grid(PCA_16sY1, PCA_16sY2, PCA_ITSY1,  PCA_ITSY2, 
                               ncol=2, nrow=2)
PCA_byyear

ggsave("PC_byyear.png", plot =PCA_byyear, device="png", width=12, height=10, units="in",dpi=600)



## Hierarchical Clustering #
#Generate the distance matrix based on Aitchinson simplex
dist_16sY1<-dist(otu.n0.clr_16sY1, method='euclidean')
dist_ITSY1<-dist(otu.n0.clr_ITSY1, method='euclidean')

dist_16sY2<-dist(otu.n0.clr_16sY2, method='euclidean')
dist_ITSY2<-dist(otu.n0.clr_ITSY2, method='euclidean')

#Cluster the data
hc_16sY1<-hclust(dist_16sY1, method = 'ward.D2')
hc_ITSY1<-hclust(dist_ITSY1, method = 'ward.D2')

hc_16sY2<-hclust(dist_16sY2, method = 'ward.D2')
hc_ITSY2<-hclust(dist_ITSY2, method = 'ward.D2')


#Plot the dendrogram
metadata_16sY1$Facility<-as.factor(metadata_16sY1$Facility) #Make Facility a factor
metadata_16sY1$L..monocytogenes<-as.factor(metadata_16sY1$L..monocytogenes) #make L.mono a factor

metadata_16sY2$Facility<-as.factor(metadata_16sY2$Facility) #Make Facility a factor
metadata_16sY2$L..monocytogenes<-as.factor(metadata_16sY2$L..monocytogenes) #make L.mono a factor

metadata_ITSY1$Facility<-as.factor(metadata_ITSY1$Facility) #Make Facility a factor
metadata_ITSY1$L..monocytogenes<-as.factor(metadata_ITSY1$L..monocytogenes) #make L.mono a factor

metadata_ITSY2$Facility<-as.factor(metadata_ITSY2$Facility) #Make Facility a factor
metadata_ITSY2$L..monocytogenes<-as.factor(metadata_ITSY2$L..monocytogenes) #make L.mono a factor

#Make vector with colors y facility and order them as in dendrogram
cols_fac <- c('#420A68','#BB3754','#FCA50A') #colors for facilities
col_16sY1_fac <- (cols_fac[metadata_16sY1$Facility])[order.dendrogram(as.dendrogram(hc_16sY1))]
col_ITSY1_fac <- (cols_fac[metadata_ITSY1$Facility])[order.dendrogram(as.dendrogram(hc_ITSY1))]

col_16sY2_fac <- (cols_fac[metadata_16sY2$Facility])[order.dendrogram(as.dendrogram(hc_16sY2))]
col_ITSY2_fac <- (cols_fac[metadata_ITSY2$Facility])[order.dendrogram(as.dendrogram(hc_ITSY2))]

#Make vector with colors L. monocytogenes and order them as in dendrogram
cols_Lm <- c('#F8CD9C','#EA7580') #colors for year
col_16sY1_Lm <- (cols_Lm[metadata_16sY1$L..monocytogenes])[order.dendrogram(as.dendrogram(hc_16sY1))]
col_ITSY1_Lm <- (cols_Lm[metadata_ITSY1$L..monocytogenes])[order.dendrogram(as.dendrogram(hc_ITSY1))]

col_16sY2_Lm <- (cols_Lm[metadata_16sY2$L..monocytogenes])[order.dendrogram(as.dendrogram(hc_16sY2))]
col_ITSY2_Lm <- (cols_Lm[metadata_ITSY2$L..monocytogenes])[order.dendrogram(as.dendrogram(hc_ITSY2))]


#Make dendrogram object
dendro_16sY1<- as.dendrogram(hc_16sY1) %>%
  set("labels", "")
dendro_16sY2<- as.dendrogram(hc_16sY2) %>%
  set("labels", "")
dendro_ITSY1<- as.dendrogram(hc_ITSY1) %>%
  set("labels", "")
dendro_ITSY2<- as.dendrogram(hc_ITSY2) %>%
  set("labels", "")


#plot and save
png("Dendrogram - 16sY1.png", width = 8, height = 4, units = 'in', res=600)
par(mar=c(3,6,3,4), xpd=T)
plot(dendro_16sY1, horiz=FALSE, main= "Bacteria Y1 - ", axes=FALSE)
colored_bars(colors = cbind(col_16sY1_Lm,col_16sY1_fac), dend=dendro_16sY1, rowLabels = c('L. monocytogenes','Facility'), 
             sort_by_labels_order = FALSE)
legend("topright", inset = c(-0.1,-0.1), legend = c('-','+','F1','F2','F3'), 
       fill = c('#F8CD9C','#EA7580','#420A68','#BB3754','#FCA50A'), bty='n')
dev.off() 

png("Dendrogram - 16sY2.png", width = 8, height = 4, units = 'in', res=600)
par(mar=c(3,6,3,4), xpd=T)
plot(dendro_16sY2, horiz=FALSE, main= "Bacteria Y2 - ", axes=FALSE)
colored_bars(colors = cbind(col_16sY2_Lm,col_16sY2_fac), dend=dendro_16sY2, rowLabels = c('L. monocytogenes','Facility'), 
             sort_by_labels_order = FALSE)
legend("topright", inset = c(-0.1,-0.1), legend = c('-','+','F1','F2','F3'), 
       fill = c('#F8CD9C','#EA7580','#420A68','#BB3754','#FCA50A'), bty='n')
dev.off()  

png("Dendrogram - ITSY1.png", width = 8, height = 4, units = 'in', res=600)
par(mar=c(3,6,3,4), xpd=T)
plot(dendro_ITSY1, horiz=FALSE, main= "Fungi Y1", axes=FALSE)
colored_bars(colors = cbind(col_ITSY1_Lm,col_ITSY1_fac), dend=dendro_ITSY1, rowLabels = c('L. monocytogenes','Facility'), 
             sort_by_labels_order = FALSE)
legend("topright", inset = c(-0.1,-0.1), legend = c('-','+','F1','F2','F3'), 
       fill = c('#F8CD9C','#EA7580','#420A68','#BB3754','#FCA50A'), bty='n')
dev.off() 

png("Dendrogram - ITSY2.png", width = 8, height = 4, units = 'in', res=600)
par(mar=c(3,6,3,4), xpd=T)
plot(dendro_ITSY2, horiz=FALSE, main= "Fungi Y2 ", axes=FALSE)
colored_bars(colors = cbind(col_ITSY2_Lm,col_ITSY2_fac), dend=dendro_ITSY2, rowLabels = c('L. monocytogenes','Facility'), 
             sort_by_labels_order = FALSE)
legend("topright", inset = c(-0.1,-0.1), legend = c('-','+','F1','F2','F3'), 
       fill = c('#F8CD9C','#EA7580','#420A68','#BB3754','#FCA50A'), bty='n')
dev.off()  

 

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


#### Microbiota and mycobiota composition plots ####
META_16sY1 = sample_data(metadata_16sY1)
META_16sY2 = sample_data(metadata_16sY2)

META_ITSY1 = sample_data(metadata_ITSY1)
META_ITSY2 = sample_data(metadata_ITSY2)

#Note: used compositional approach to transform the sample counts to compositions. 

#Transform sample counts into compositions
otu.n0.acomp_16sY1<-as.data.frame(acomp(t(otu.n0_16sY1)), total=1)
otu.n0.acomp_16sY2<-as.data.frame(acomp(t(otu.n0_16sY2)), total=1)

otu.n0.acomp_ITSY1<-as.data.frame(acomp(t(otu.n0_ITSY1)), total=1)
otu.n0.acomp_ITSY2<-as.data.frame(acomp(t(otu.n0_ITSY2)), total=1)


#OTU level plots
#Make Phyloseq object for eash year
OTU_16sY1 <- otu_table(otu.n0.acomp_16sY1, taxa_are_rows = FALSE)
phyloseq_16sY1 = phyloseq(OTU_16sY1, TAX_16s, META_16sY1)
TREE_16sY1 = rtree(ntaxa(phyloseq_16sY1), rooted=TRUE, tip.label = taxa_names(phyloseq_16sY1))  
phyloseq16sY1 <- phyloseq(OTU_16sY1, TAX_16s, TREE_16sY1, META_16sY1)

OTU_16sY2 <- otu_table(otu.n0.acomp_16sY2, taxa_are_rows = FALSE)
phyloseq_16sY2 = phyloseq(OTU_16sY2, TAX_16s, META_16sY2)
TREE_16sY2 = rtree(ntaxa(phyloseq_16sY2), rooted=TRUE, tip.label = taxa_names(phyloseq_16sY2))  
phyloseq16sY2 <- phyloseq(OTU_16sY2, TAX_16s, TREE_16sY2, META_16sY2)

OTU_ITSY1 <- otu_table(otu.n0.acomp_ITSY1, taxa_are_rows = FALSE)
phyloseq_ITSY1 = phyloseq(OTU_ITSY1, TAX_ITS, META_ITSY1)
TREE_ITSY1 = rtree(ntaxa(phyloseq_ITSY1), rooted=TRUE, tip.label = taxa_names(phyloseq_ITSY1))  
phyloseqITSY1 <- phyloseq(OTU_ITSY1, TAX_ITS, TREE_ITSY1, META_ITSY1)

OTU_ITSY2 <- otu_table(otu.n0.acomp_ITSY2, taxa_are_rows = FALSE)
phyloseq_ITSY2 = phyloseq(OTU_ITSY2, TAX_ITS, META_ITSY2)
TREE_ITSY2 = rtree(ntaxa(phyloseq_ITSY2), rooted=TRUE, tip.label = taxa_names(phyloseq_ITSY2))  
phyloseqITSY2 <- phyloseq(OTU_ITSY2, TAX_ITS, TREE_ITSY2, META_ITSY2)


#Make long format table from Phyloseq object
otu_16sY1_long <- phyloseq16sY1 %>%  
  transform_sample_counts(function(x) {x * 100} ) %>% #Recalculates the relative abundance 
  psmelt() %>%  #Melts to long format
  arrange(desc(Abundance))

otu_16sY2_long <- phyloseq16sY2 %>%  
  transform_sample_counts(function(x) {x * 100} ) %>% #Recalculates the relative abundance 
  psmelt() %>%  #Melts to long format
  arrange(desc(Abundance))

otu_ITSY1_long <- phyloseqITSY1 %>%  
  transform_sample_counts(function(x) {x * 100} ) %>% #Recalculates the relative abundance 
  psmelt() %>%  #Melts to long format
  arrange(desc(Abundance))

otu_ITSY2_long <- phyloseqITSY2 %>%  
  transform_sample_counts(function(x) {x * 100} ) %>% #Recalculates the relative abundance 
  psmelt() %>%  #Melts to long format
  arrange(desc(Abundance))

#Save OTU composition as .csv
write.csv(otu_16sY1_long, "otu_16sY1.csv")
write.csv(otu_16sY2_long, "otu_16sY2.csv")
write.csv(otu_ITSY1_long, "otu_ITSY1.csv")
write.csv(otu_ITSY2_long, "otu_ITSY2.csv")

#Make vector with OTU names above 5% relative abundance in at least one sample
otu_16sY1_over5<-unique(c(otu_16sY1_long$OTU[which(otu_16sY1_long$Abundance >=5)]))
otu_16sY2_over5<-unique(c(otu_16sY2_long$OTU[which(otu_16sY2_long$Abundance >=5)])) 
otu_ITSY1_over5<-unique(c(otu_ITSY1_long$OTU[which(otu_ITSY1_long$Abundance >=5)])) 
otu_ITSY2_over5<-unique(c(otu_ITSY2_long$OTU[which(otu_ITSY2_long$Abundance >=5)]))

#Filter table to obtain only OTUs with over 5% in at least one sample
otu_16sY1_over5abund <- filter(otu_16sY1_long, OTU %in% otu_16sY1_over5)
otu_16sY2_over5abund <- filter(otu_16sY2_long, OTU %in% otu_16sY2_over5)
otu_ITSY1_over5abund <- filter(otu_ITSY1_long, OTU %in% otu_ITSY1_over5)
otu_ITSY2_over5abund <- filter(otu_ITSY2_long, OTU %in% otu_ITSY2_over5)

### Bubbleplots by facility at the OTU level
#16s
bubbleplot_otu16sY1<-ggplot(otu_16sY1_over5abund, aes(x=Facility, y=reorder(OTU,desc(OTU)), size=Abundance, color=Facility))+
  geom_point()+facet_grid(Family~., scales = "free", space = 'free')+
  scale_size(limits=c(0,100), name = "Relative abundance (%)")+
  theme(axis.title = element_blank(), axis.text.x=element_text(color='black', size=8), axis.ticks=element_line(color='black'),
        axis.text.y = element_text(color='black', size=6)) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  theme(panel.background = element_rect(fill="grey99", color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_rect(fill="white", color = "black"), strip.text.y = element_text(size=6, angle=0), 
        panel.border = element_rect(color="black", fill=NA))+
  ggtitle("Bacteria Y1 - ")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
  scale_color_viridis_d(begin = 0.2, end = 0.8, option='inferno')+
  theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))
bubbleplot_otu16sY1

bubbleplot_otu16sY2<-ggplot(otu_16sY2_over5abund, aes(x=Facility, y=reorder(OTU,desc(OTU)), size=Abundance, color=Facility))+
  geom_point()+facet_grid(Family~., scales = "free", space = 'free')+
  scale_size(limits=c(0,100), name = "Relative abundance (%)")+
  theme(axis.title = element_blank(), axis.text.x=element_text(color='black', size=8), axis.ticks=element_line(color='black'),
        axis.text.y = element_text(color='black', size=6)) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  theme(panel.background = element_rect(fill="grey99", color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_rect(fill="white", color = "black"), strip.text.y = element_text(size=6, angle=0), 
        panel.border = element_rect(color="black", fill=NA))+
  ggtitle("Bacteria Y2 - ")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
  scale_color_viridis_d(begin = 0.2, end = 0.8, option='inferno')+
  theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))
bubbleplot_otu16sY2

Bubbleplots_otu16s = plot_grid(bubbleplot_otu16sY1+theme(legend.position = 'right'),
                            bubbleplot_otu16sY2+theme(legend.position = 'right'), 
                            ncol=2, nrow=1, labels = c("A","B"), label_size = 20, vjust = 2, hjust = -1.5)
Bubbleplots_otu16s

ggsave("Bubbleplot_otu16s.png", plot=Bubbleplots_otu16s, device="png", width=12, height=10, units="in", dpi=600)


#ITS
bubbleplot_otuITSY1<-ggplot(otu_ITSY1_over5abund, aes(x=Facility, y=reorder(OTU,desc(OTU)), size=Abundance, color=Facility))+
  geom_point()+facet_grid(Family~., scales = "free", space = 'free')+
  scale_size(limits=c(0,100), name = "Relative abundance (%)")+
  theme(axis.title = element_blank(), axis.text.x=element_text(color='black', size=8), axis.ticks=element_line(color='black'),
        axis.text.y = element_text(color='black', size=6)) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  theme(panel.background = element_rect(fill="grey99", color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_rect(fill="white", color = "black"), strip.text.y = element_text(size=6, angle=0), 
        panel.border = element_rect(color="black", fill=NA))+
  ggtitle("Fungi Y1")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
  scale_color_viridis_d(begin = 0.2, end = 0.8, option='inferno')+
  theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))
bubbleplot_otuITSY1

bubbleplot_otuITSY2<-ggplot(otu_ITSY2_over5abund, aes(x=Facility, y=reorder(OTU,desc(OTU)), size=Abundance, color=Facility))+
  geom_point()+facet_grid(Family~., scales = "free", space = 'free')+
  scale_size(limits=c(0,100), name = "Relative abundance (%)")+
  theme(axis.title = element_blank(), axis.text.x=element_text(color='black', size=8), axis.ticks=element_line(color='black'),
        axis.text.y = element_text(color='black', size=6)) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  theme(panel.background = element_rect(fill="grey99", color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_rect(fill="white", color = "black"), strip.text.y = element_text(size=6, angle=0), 
        panel.border = element_rect(color="black", fill=NA))+
  ggtitle("Fungi Y2")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
  scale_color_viridis_d(begin = 0.2, end = 0.8, option='inferno')+
  theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))
bubbleplot_otuITSY2

Bubbleplots_otuITS = plot_grid(bubbleplot_otuITSY1+theme(legend.position = 'right'),
                               bubbleplot_otuITSY2+theme(legend.position = 'right'), 
                               ncol=2, nrow=1, labels = c("C","D"), label_size = 20, vjust = 2, hjust = -1.5)
Bubbleplots_otuITS

ggsave("Bubbleplot_otuITS.png", plot=Bubbleplots_otuITS, device="png", width=12, height=10, units="in", dpi=600)


#### Core microbiome at the OTU level by facility ####
#With the microbiome package
library(microbiome)
library(knitr)

# Subset data for each facility
F1_16sY1 <- subset_samples(phyloseq16sY1, Facility == "F1") 
F2_16sY1 <- subset_samples(phyloseq16sY1, Facility == "F2") 
F3_16sY1 <- subset_samples(phyloseq16sY1, Facility == "F3") 

F1_16sY2 <- subset_samples(phyloseq16sY2, Facility == "F1") 
F2_16sY2 <- subset_samples(phyloseq16sY2, Facility == "F2") 
F3_16sY2 <- subset_samples(phyloseq16sY2, Facility == "F3") 

F1_ITSY1 <- subset_samples(phyloseqITSY1, Facility == "F1") 
F2_ITSY1 <- subset_samples(phyloseqITSY1, Facility == "F2") 
F3_ITSY1 <- subset_samples(phyloseqITSY1, Facility == "F3") 

F1_ITSY2 <- subset_samples(phyloseqITSY2, Facility == "F1") 
F2_ITSY2 <- subset_samples(phyloseqITSY2, Facility == "F2") 
F3_ITSY2 <- subset_samples(phyloseqITSY2, Facility == "F3")




#Core microbiota by facility
#Families detected in at least 75% of the samples with a relative abundance threshold value above 0.01%
F1_16sY1_core <- core(F1_16sY1, detection = 0.001, prevalence = 0.75) 
F2_16sY1_core <- core(F2_16sY1, detection = 0.001, prevalence = 0.75)
F3_16sY1_core <- core(F3_16sY1, detection = 0.001, prevalence = 0.75)

F1_16sY2_core <- core(F1_16sY2, detection = 0.001, prevalence = 0.75) 
F2_16sY2_core <- core(F2_16sY2, detection = 0.001, prevalence = 0.75)
F3_16sY2_core <- core(F3_16sY2, detection = 0.001, prevalence = 0.75)

F1_ITSY1_core <- core(F1_ITSY1, detection = 0.001, prevalence = 0.75) 
F2_ITSY1_core <- core(F2_ITSY1, detection = 0.001, prevalence = 0.75)
F3_ITSY1_core <- core(F3_ITSY1, detection = 0.001, prevalence = 0.75)

F1_ITSY2_core <- core(F1_ITSY2, detection = 0.001, prevalence = 0.75) 
F2_ITSY2_core <- core(F2_ITSY2, detection = 0.001, prevalence = 0.75)
F3_ITSY2_core <- core(F3_ITSY2, detection = 0.001, prevalence = 0.75)


# get the taxonomy data
F1_16sY1_tax <- as.data.frame(tax_table(F1_16sY1_core)) #17 taxa
F2_16sY1_tax <- as.data.frame(tax_table(F2_16sY1_core)) #12 taxa
F3_16sY1_tax <- as.data.frame(tax_table(F3_16sY1_core)) #12 taxa

F1_16sY2_tax <- as.data.frame(tax_table(F1_16sY2_core)) #16 taxa
F2_16sY2_tax <- as.data.frame(tax_table(F2_16sY2_core)) #31 taxa
F3_16sY2_tax <- as.data.frame(tax_table(F3_16sY2_core)) #26 taxa

F1_ITSY1_tax <- as.data.frame(tax_table(F1_ITSY1_core)) #18 taxa
F2_ITSY1_tax <- as.data.frame(tax_table(F2_ITSY1_core)) #16 taxa
F3_ITSY1_tax <- as.data.frame(tax_table(F3_ITSY1_core)) #12 taxa

F1_ITSY2_tax <- as.data.frame(tax_table(F1_ITSY2_core)) #25 taxa
F2_ITSY2_tax <- as.data.frame(tax_table(F2_ITSY2_core)) #20 taxa
F3_ITSY2_tax <- as.data.frame(tax_table(F3_ITSY2_core)) #12 taxa

#Extract family
F1_16sY1_family<-F1_16sY1_tax$Family
F2_16sY1_family<-F2_16sY1_tax$Family
F3_16sY1_family<-F3_16sY1_tax$Family

F1_16sY2_family<-F1_16sY2_tax$Family
F2_16sY2_family<-F2_16sY2_tax$Family
F3_16sY2_family<-F3_16sY2_tax$Family

F1_ITSY1_family<-F1_ITSY1_tax$Family
F2_ITSY1_family<-F2_ITSY1_tax$Family
F3_ITSY1_family<-F3_ITSY1_tax$Family

F1_ITSY2_family<-F1_ITSY2_tax$Family
F2_ITSY2_family<-F2_ITSY2_tax$Family
F3_ITSY2_family<-F3_ITSY2_tax$Family

#Make list with all core taxa for each category
coretaxa_16sY1_list<-list(F1=F1_16sY1_family, F2=F2_16sY1_family, F3=F3_16sY1_family)
coretaxa_16sY2_list<-list(F1=F1_16sY2_family, F2=F2_16sY2_family, F3=F3_16sY2_family)
coretaxa_ITSY1_list<-list(F1=F1_ITSY1_family, F2=F2_ITSY1_family, F3=F3_ITSY1_family)
coretaxa_ITSY2_list<-list(F1=F1_ITSY2_family, F2=F2_ITSY2_family, F3=F3_ITSY2_family)

#Make Venn diagrams
library(ggvenn)

venn_16sY1<-ggvenn(coretaxa_16sY1_list, show_percentage = FALSE,
                   fill_color = c('#420A68','#BB3754','#FCA50A'), set_name_size = 4, stroke_size = 0.5)+
  ggtitle("Bacteria Y1")

venn_16sY2<-ggvenn(coretaxa_16sY2_list, show_percentage = FALSE,
                   fill_color = c('#420A68','#BB3754','#FCA50A'), set_name_size = 4, stroke_size = 0.5)+
  ggtitle("Bacteria Y2")

venn_ITSY1<-ggvenn(coretaxa_ITSY1_list, show_percentage = FALSE,
                   fill_color = c('#420A68','#BB3754','#FCA50A'), set_name_size = 4, stroke_size = 0.5)+
  ggtitle("Fungi Y1")

venn_ITSY2<-ggvenn(coretaxa_ITSY2_list, show_percentage = FALSE,
                   fill_color = c('#420A68','#BB3754','#FCA50A'), set_name_size = 4, stroke_size = 0.5)+
  ggtitle("Fungi Y2")


#Combine Venn diagrams
Venn_16s_Fac = plot_grid(venn_16sY1, venn_16sY2, 
                         ncol=2, nrow=1, labels = c("A","B"), label_size = 20, vjust = 2, hjust = -1.5)
Venn_16s_Fac

ggsave("Venn_core16s_facility.png", plot=Venn_16s_Fac, device="png", width=8, height=4, units="in", dpi=600)


Venn_ITS_Fac = plot_grid(venn_ITSY1, venn_ITSY2, 
                         ncol=2, nrow=1, labels = c("A","B"), label_size = 20, vjust = 2, hjust = -1.5)
Venn_ITS_Fac

ggsave("Venn_coreITS_facility.png", plot=Venn_ITS_Fac, device="png", width=8, height=4, units="in", dpi=600)

#Print taxa that are shared across facilities
intersect(intersect(F1_16sY1_family, F2_16sY1_family),F3_16sY1_family)
intersect(intersect(F1_16sY2_family, F2_16sY2_family),F3_16sY2_family)
intersect(intersect(F1_ITSY1_family, F2_ITSY1_family),F3_ITSY1_family)
intersect(intersect(F1_ITSY2_family, F2_ITSY2_family),F3_ITSY2_family)

#Compare core taxa by year for each facility
coretaxa_16sF1_list<-list(Y1=F1_16sY1_family, Y2=F1_16sY2_family)
coretaxa_16sF2_list<-list(Y1=F2_16sY1_family, Y2=F2_16sY2_family)
coretaxa_16sF3_list<-list(Y1=F3_16sY1_family, Y2=F3_16sY2_family)

coretaxa_ITSF1_list<-list(Y1=F1_ITSY1_family, Y2=F1_ITSY2_family)
coretaxa_ITSF2_list<-list(Y1=F2_ITSY1_family, Y2=F2_ITSY2_family)
coretaxa_ITSF3_list<-list(Y1=F3_ITSY1_family, Y2=F3_ITSY2_family)

#Make Venn diagrams

venn_16sF1<-ggvenn(coretaxa_16sF1_list, show_percentage = FALSE,
                   fill_color = c('#420A68','#8D6CA4'), set_name_size = 4, stroke_size = 0.5)

venn_16sF2<-ggvenn(coretaxa_16sF2_list, show_percentage = FALSE,
                   fill_color = c('#BB3754','#D68798'), set_name_size = 4, stroke_size = 0.5)

venn_16sF3<-ggvenn(coretaxa_16sF3_list, show_percentage = FALSE,
                   fill_color = c('#FCA50A','#FDC96C'), set_name_size = 4, stroke_size = 0.5)


venn_ITSF1<-ggvenn(coretaxa_ITSF1_list, show_percentage = FALSE,
                   fill_color = c('#420A68','#8D6CA4'), set_name_size = 4, stroke_size = 0.5)

venn_ITSF2<-ggvenn(coretaxa_ITSF2_list, show_percentage = FALSE,
                   fill_color = c('#BB3754','#D68798'), set_name_size = 4, stroke_size = 0.5)

venn_ITSF3<-ggvenn(coretaxa_ITSF3_list, show_percentage = FALSE,
                   fill_color = c('#FCA50A','#FDC96C'), set_name_size = 4, stroke_size = 0.5)


#Combine Venn diagrams by year
Venn_16s_Year = plot_grid(venn_16sF1, venn_16sF2, venn_16sF3,
                          ncol=3, nrow=2, labels = c("F1","F2","F3"), label_size = 20, vjust = 2, hjust = -1.5)
Venn_16s_Year

ggsave("Venn_core16s_year.png", plot=Venn_16s_Year, device="png", width=8, height=4, units="in", dpi=600)


Venn_ITS_Year = plot_grid(venn_ITSF1, venn_ITSF2, venn_ITSF3,
                          ncol=3, nrow=2, labels = c("F1","F2","F3"), label_size = 20, vjust = 2, hjust = -1.5)
Venn_ITS_Year

ggsave("Venn_coreITS_year.png", plot=Venn_ITS_Year, device="png", width=8, height=4, units="in", dpi=600)

#Print taxa that are shared across seasons for each facility
intersect(F1_16sY1_family,F1_16sY2_family)
intersect(F2_16sY1_family,F2_16sY2_family)
intersect(F3_16sY1_family,F3_16sY2_family)

intersect(F1_ITSY1_family,F1_ITSY2_family)
intersect(F2_ITSY1_family,F2_ITSY2_family)
intersect(F3_ITSY1_family,F3_ITSY2_family)



#### Core microbiome at the familiy by facility ####
#With the microbiome package
library(microbiome)
library(knitr)

#Make phyloseq with relative abundances at the Family level
phyloseq_family_Y1 <- phyloseq16sY1 %>%  
  tax_glom(taxrank = "Family")

phyloseq_family_16sY2 <- phyloseq16sY2 %>%  
  tax_glom(taxrank = "Family")

phyloseq_family_ITSY1 <- phyloseqITSY1 %>%  
  tax_glom(taxrank = "Family")

phyloseq_family_ITSY2 <- phyloseqITSY2 %>%  
  tax_glom(taxrank = "Family")

# Subset data for each facility
F1_16sY1 <- subset_samples(phyloseq_family_16sY1, Facility == "F1") 
F2_16sY1 <- subset_samples(phyloseq_family_16sY1, Facility == "F2") 
F3_16sY1 <- subset_samples(phyloseq_family_16sY1, Facility == "F3") 

F1_16sY2 <- subset_samples(phyloseq_family_16sY2, Facility == "F1") 
F2_16sY2 <- subset_samples(phyloseq_family_16sY2, Facility == "F2") 
F3_16sY2 <- subset_samples(phyloseq_family_16sY2, Facility == "F3") 

F1_ITSY1 <- subset_samples(phyloseq_family_ITSY1, Facility == "F1") 
F2_ITSY1 <- subset_samples(phyloseq_family_ITSY1, Facility == "F2") 
F3_ITSY1 <- subset_samples(phyloseq_family_ITSY1, Facility == "F3") 

F1_ITSY2 <- subset_samples(phyloseq_family_ITSY2, Facility == "F1") 
F2_ITSY2 <- subset_samples(phyloseq_family_ITSY2, Facility == "F2") 
F3_ITSY2 <- subset_samples(phyloseq_family_ITSY2, Facility == "F3") 


#Core microbiota by facility
#Families detected in at least 75% of the samples with a relative abundance threshold value above 0.01%
F1_16sY1_core <- core(F1_16sY1, detection = 0.01, prevalence = 0.75) 
F2_16sY1_core <- core(F2_16sY1, detection = 0.01, prevalence = 0.75)
F3_16sY1_core <- core(F3_16sY1, detection = 0.01, prevalence = 0.75)

F1_16sY2_core <- core(F1_16sY2, detection = 0.01, prevalence = 0.75) 
F2_16sY2_core <- core(F2_16sY2, detection = 0.01, prevalence = 0.75)
F3_16sY2_core <- core(F3_16sY2, detection = 0.01, prevalence = 0.75)

F1_ITSY1_core <- core(F1_ITSY1, detection = 0.01, prevalence = 0.75) 
F2_ITSY1_core <- core(F2_ITSY1, detection = 0.01, prevalence = 0.75)
F3_ITSY1_core <- core(F3_ITSY1, detection = 0.01, prevalence = 0.75)

F1_ITSY2_core <- core(F1_ITSY2, detection = 0.01, prevalence = 0.75) 
F2_ITSY2_core <- core(F2_ITSY2, detection = 0.01, prevalence = 0.75)
F3_ITSY2_core <- core(F3_ITSY2, detection = 0.01, prevalence = 0.75)


# get the taxonomy data
F1_16sY1_tax <- as.data.frame(tax_table(F1_16sY1_core)) #9 taxa
F2_16sY1_tax <- as.data.frame(tax_table(F2_16sY1_core)) #8 taxa
F3_16sY1_tax <- as.data.frame(tax_table(F3_16sY1_core)) #9 taxa

F1_16sY2_tax <- as.data.frame(tax_table(F1_16sY2_core)) #13 taxa
F2_16sY2_tax <- as.data.frame(tax_table(F2_16sY2_core)) #16 taxa
F3_16sY2_tax <- as.data.frame(tax_table(F3_16sY2_core)) #13 taxa

F1_ITSY1_tax <- as.data.frame(tax_table(F1_ITSY1_core)) #8 taxa
F2_ITSY1_tax <- as.data.frame(tax_table(F2_ITSY1_core)) #3 taxa
F3_ITSY1_tax <- as.data.frame(tax_table(F3_ITSY1_core)) #8 taxa

F1_ITSY2_tax <- as.data.frame(tax_table(F1_ITSY2_core)) #7 taxa
F2_ITSY2_tax <- as.data.frame(tax_table(F2_ITSY2_core)) #6 taxa
F3_ITSY2_tax <- as.data.frame(tax_table(F3_ITSY2_core)) #8 taxa

#Extract family
F1_16sY1_family<-F1_16sY1_tax$Family
F2_16sY1_family<-F2_16sY1_tax$Family
F3_16sY1_family<-F3_16sY1_tax$Family

F1_16sY2_family<-F1_16sY2_tax$Family
F2_16sY2_family<-F2_16sY2_tax$Family
F3_16sY2_family<-F3_16sY2_tax$Family

F1_ITSY1_family<-F1_ITSY1_tax$Family
F2_ITSY1_family<-F2_ITSY1_tax$Family
F3_ITSY1_family<-F3_ITSY1_tax$Family

F1_ITSY2_family<-F1_ITSY2_tax$Family
F2_ITSY2_family<-F2_ITSY2_tax$Family
F3_ITSY2_family<-F3_ITSY2_tax$Family

#Make list with all core taxa for each category
coretaxa_16sY1_list<-list(F1=F1_16sY1_family, F2=F2_16sY1_family, F3=F3_16sY1_family)
coretaxa_16sY2_list<-list(F1=F1_16sY2_family, F2=F2_16sY2_family, F3=F3_16sY2_family)
coretaxa_ITSY1_list<-list(F1=F1_ITSY1_family, F2=F2_ITSY1_family, F3=F3_ITSY1_family)
coretaxa_ITSY2_list<-list(F1=F1_ITSY2_family, F2=F2_ITSY2_family, F3=F3_ITSY2_family)

#Make Venn diagrams
library(ggvenn)

venn_16sY1<-ggvenn(coretaxa_16sY1_list, show_percentage = FALSE,
                   fill_color = c('#420A68','#BB3754','#FCA50A'), set_name_size = 4, stroke_size = 0.5)+
  ggtitle("Bacteria Y1")

venn_16sY2<-ggvenn(coretaxa_16sY2_list, show_percentage = FALSE,
                   fill_color = c('#420A68','#BB3754','#FCA50A'), set_name_size = 4, stroke_size = 0.5)+
  ggtitle("Bacteria Y2")

venn_ITSY1<-ggvenn(coretaxa_ITSY1_list, show_percentage = FALSE,
                   fill_color = c('#420A68','#BB3754','#FCA50A'), set_name_size = 4, stroke_size = 0.5)+
  ggtitle("Fungi Y1")

venn_ITSY2<-ggvenn(coretaxa_ITSY2_list, show_percentage = FALSE,
                   fill_color = c('#420A68','#BB3754','#FCA50A'), set_name_size = 4, stroke_size = 0.5)+
  ggtitle("Fungi Y2")


#Combine Venn diagrams
Venn_16s_Fac = plot_grid(venn_16sY1, venn_16sY2, 
                         ncol=2, nrow=1, labels = c("A","B"), label_size = 20, vjust = 2, hjust = -1.5)
Venn_16s_Fac

ggsave("Venn_core16s_facility.png", plot=Venn_16s_Fac, device="png", width=8, height=4, units="in", dpi=600)


Venn_ITS_Fac = plot_grid(venn_ITSY1, venn_ITSY2, 
                         ncol=2, nrow=1, labels = c("A","B"), label_size = 20, vjust = 2, hjust = -1.5)
Venn_ITS_Fac

ggsave("Venn_coreITS_facility.png", plot=Venn_ITS_Fac, device="png", width=8, height=4, units="in", dpi=600)

#Print taxa that are shared across facilities
intersect(intersect(F1_16sY1_family, F2_16sY1_family),F3_16sY1_family)
intersect(intersect(F1_16sY2_family, F2_16sY2_family),F3_16sY2_family)
intersect(intersect(F1_ITSY1_family, F2_ITSY1_family),F3_ITSY1_family)
intersect(intersect(F1_ITSY2_family, F2_ITSY2_family),F3_ITSY2_family)

#Compare core taxa by year for each facility
coretaxa_16sF1_list<-list(Y1=F1_16sY1_family, Y2=F1_16sY2_family)
coretaxa_16sF2_list<-list(Y1=F2_16sY1_family, Y2=F2_16sY2_family)
coretaxa_16sF3_list<-list(Y1=F3_16sY1_family, Y2=F3_16sY2_family)

coretaxa_ITSF1_list<-list(Y1=F1_ITSY1_family, Y2=F1_ITSY2_family)
coretaxa_ITSF2_list<-list(Y1=F2_ITSY1_family, Y2=F2_ITSY2_family)
coretaxa_ITSF3_list<-list(Y1=F3_ITSY1_family, Y2=F3_ITSY2_family)

#Make Venn diagrams

venn_16sF1<-ggvenn(coretaxa_16sF1_list, show_percentage = FALSE,
                   fill_color = c('#420A68','#8D6CA4'), set_name_size = 4, stroke_size = 0.5)

venn_16sF2<-ggvenn(coretaxa_16sF2_list, show_percentage = FALSE,
                   fill_color = c('#BB3754','#D68798'), set_name_size = 4, stroke_size = 0.5)

venn_16sF3<-ggvenn(coretaxa_16sF3_list, show_percentage = FALSE,
                   fill_color = c('#FCA50A','#FDC96C'), set_name_size = 4, stroke_size = 0.5)


venn_ITSF1<-ggvenn(coretaxa_ITSF1_list, show_percentage = FALSE,
                   fill_color = c('#420A68','#8D6CA4'), set_name_size = 4, stroke_size = 0.5)

venn_ITSF2<-ggvenn(coretaxa_ITSF2_list, show_percentage = FALSE,
                   fill_color = c('#BB3754','#D68798'), set_name_size = 4, stroke_size = 0.5)

venn_ITSF3<-ggvenn(coretaxa_ITSF3_list, show_percentage = FALSE,
                   fill_color = c('#FCA50A','#FDC96C'), set_name_size = 4, stroke_size = 0.5)


#Combine Venn diagrams by year
Venn_16s_Year = plot_grid(venn_16sF1, venn_16sF2, venn_16sF3,
                          ncol=3, nrow=2, labels = c("F1","F2","F3"), label_size = 20, vjust = 2, hjust = -1.5)
Venn_16s_Year

ggsave("Venn_core16s_year.png", plot=Venn_16s_Year, device="png", width=8, height=4, units="in", dpi=600)


Venn_ITS_Year = plot_grid(venn_ITSF1, venn_ITSF2, venn_ITSF3,
                          ncol=3, nrow=2, labels = c("F1","F2","F3"), label_size = 20, vjust = 2, hjust = -1.5)
Venn_ITS_Year

ggsave("Venn_coreITS_year.png", plot=Venn_ITS_Year, device="png", width=8, height=4, units="in", dpi=600)

#Print taxa that are shared across seasons for each facility
intersect(F1_16sY1_family,F1_16sY2_family)
intersect(F2_16sY1_family,F2_16sY2_family)
intersect(F3_16sY1_family,F3_16sY2_family)

intersect(F1_ITSY1_family,F1_ITSY2_family)
intersect(F2_ITSY1_family,F2_ITSY2_family)
intersect(F3_ITSY1_family,F3_ITSY2_family)

## Core microbiomes using Poisson distribution
#OTU level



# #Make barplots with abundance of core microbiome by facility
# #Open in Excel the files family_16sY1.csv,  family_16sY2.csv, family_ITSY1.csv, family_ITSY2.csv for the core bactera shared by all facilities in 2 years.
# #Copy and paste in a new sheet and save as .csv
# 
# family_16sY1_core<-read.csv('core_16sY1.csv', header = TRUE)
# family_16sY2_core<-read.csv('core_16sY2.csv', header = TRUE)
# family_ITSY1_core<-read.csv('core_ITSY1.csv', header = TRUE)
# family_ITSY2_core<-read.csv('core_ITSY2.csv', header = TRUE)
# 
# #Add year for all observations
# family_16sY1_core$Year<-rep("Y1", 1170)
# family_16sY2_core$Year<-rep("Y2", 1070)
# family_ITSY1_core$Year<-rep("Y1", 936)
# family_ITSY2_core$Year<-rep("Y2", 856)
# 
# #Combine the two years data
# family_16s_core<-bind_rows(family_16sY1_core, family_16sY2_core)
# family_ITS_core<-bind_rows(family_ITSY1_core, family_ITSY2_core)
# 
# #Plot
# #16s
# boxplot_16s_core<-ggplot(family_16s_core, aes(x=Facility, y=Abundance))+
#   geom_boxplot(color='grey')+ facet_wrap(~Family+Year, scales = 'free_y', labeller = label_wrap_gen(), nrow=5, ncol=4)+
#   geom_jitter(aes(color=Facility))+
#   scale_size(range=c(.1,10), name = "Relative abundance (%)")+
#   theme(axis.title.x = element_blank(), axis.text.x=element_text(color='black', size=10), axis.ticks=element_line(color='black'),
#         axis.title.y=element_text(size=15),axis.text.y = element_text(color='black')) + 
#   guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
#   ylab("Relative Abundance") + 
#   theme(panel.background = element_rect(fill="grey99", color =NA),
#         plot.background = element_rect(fill="transparent", color =NA)) +
#   theme(strip.background= element_rect(fill="white", color = "black"), strip.text = element_text(size=7, face='bold'), 
#         panel.border = element_rect(color="black", fill=NA))+
#   ggtitle("Core Bacteria")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
#   scale_color_viridis_d(begin = 0.2, end = 0.8, option='inferno')+
#   theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))
# boxplot_16s_core
# 
# ggsave("Boxplot_core_16s.png", plot=boxplot_16s_core, device="png", width=10, height=10, units="in", dpi=600)
# 
# #ITS
# boxplot_ITS_core<-ggplot(family_ITS_core, aes(x=Facility, y=Abundance))+
#   geom_boxplot(color='grey')+ facet_wrap(~Family+Year, scales = 'free_y', labeller = label_wrap_gen(), nrow=4, ncol=4)+
#   geom_jitter(aes(color=Facility))+
#   scale_size(range=c(.1,10), name = "Relative abundance (%)")+
#   theme(axis.title.x = element_blank(), axis.text.x=element_text(color='black', size=10), axis.ticks=element_line(color='black'),
#         axis.title.y=element_text(size=15),axis.text.y = element_text(color='black')) + 
#   guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
#   ylab("Relative Abundance") + 
#   theme(panel.background = element_rect(fill="grey99", color =NA),
#         plot.background = element_rect(fill="transparent", color =NA)) +
#   theme(strip.background= element_rect(fill="white", color = "black"), strip.text = element_text(size=7, face='bold'), 
#         panel.border = element_rect(color="black", fill=NA))+
#   ggtitle("Core Fungi")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
#   scale_color_viridis_d(begin = 0.2, end = 0.8, option='inferno')+
#   theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))
# boxplot_ITS_core
# 
# ggsave("Boxplot_core_ITS.png", plot=boxplot_ITS_core, device="png", width=10, height=10, units="in", dpi=600)


#### Network analysis ####

# Using SpiecEasi method - adapted from turorial: http://psbweb05.psb.ugent.be/conet/microbialnetworks/spieceasi.php

#Install SpiecEasi package
#library(devtools)
#install_github("zdk123/SpiecEasi")
library(SpiecEasi)

#Install packag seqtime
#install_github("hallucigenia-sparsa/seqtime") 
library(seqtime)

F1_16sY1_count <- subset_samples(physeq_16sY1, Facility == "F1") 
F2_16sY1_count <- subset_samples(physeq_16sY1, Facility == "F2") 
F3_16sY1_count <- subset_samples(physeq_16sY1, Facility == "F3") 

#Extract otu table and taxon table from each phyloseq by facility
otus_F1_16sY1=otu_table(F1_16sY1_count)
taxa_F1_16sY1=tax_table(F1_16sY1_count)

otus_F2_16sY1=otu_table(F2_16sY1_count)
taxa_F2_16sY1=tax_table(F2_16sY1_count)

otus_F3_16sY1=otu_table(F3_16sY1_count)
taxa_F3_16sY1=tax_table(F3_16sY1_count)

#Filter OTU table to reduce sparcity of the dataset
F1_16sY1_filterobj=filterTaxonMatrix(otus_F1_16sY1, minocc=20, keepSum = TRUE, return.filtered.indices = TRUE)
F1_16sY1_otus.f=F1_16sY1_filterobj$mat
F1_16sY1_taxa.f=taxa_F1_16sY1[setdiff(1:nrow(taxa_F1_16sY1),F1_16sY1_filterobj$filtered.indices),]
F1_16sY1_taxa.f<-F1_16sY1_taxa.f[,-6]
F1_16sY1_dummyTaxonomy=c(rep("Dummy", 5))
F1_16sY1_taxa.f=rbind(F1_16sY1_taxa.f,F1_16sY1_dummyTaxonomy)
rownames(F1_16sY1_taxa.f)[nrow(F1_16sY1_taxa.f)]="OTU_Filt"
rownames(F1_16sY1_otus.f)[nrow(F1_16sY1_otus.f)]="0TU_Filt"
tail(F1_16sY1_otus.f)

F2_16sY1_filterobj=filterTaxonMatrix(otus_F2_16sY1, minocc=20, keepSum = TRUE, return.filtered.indices = TRUE)
F2_16sY1_otus.f=F2_16sY1_filterobj$mat
F2_16sY1_taxa.f=taxa_F2_16sY1[setdiff(1:nrow(taxa_F2_16sY1),F2_16sY1_filterobj$filtered.indices),]
F2_16sY1_taxa.f<-F2_16sY1_taxa.f[,-6]
F2_16sY1_dummyTaxonomy=c(rep("Dummy", 5))
F2_16sY1_taxa.f=rbind(F2_16sY1_taxa.f,F2_16sY1_dummyTaxonomy)
rownames(F2_16sY1_taxa.f)[nrow(F2_16sY1_taxa.f)]="OTU_Filt"
rownames(F2_16sY1_otus.f)[nrow(F2_16sY1_otus.f)]="0TU_Filt"
tail(F2_16sY1_otus.f)

F3_16sY1_filterobj=filterTaxonMatrix(otus_F3_16sY1, minocc=20, keepSum = TRUE, return.filtered.indices = TRUE)
F3_16sY1_otus.f=F3_16sY1_filterobj$mat
F3_16sY1_taxa.f=taxa_F3_16sY1[setdiff(1:nrow(taxa_F3_16sY1),F3_16sY1_filterobj$filtered.indices),]
F3_16sY1_taxa.f<-F3_16sY1_taxa.f[,-6]
F3_16sY1_dummyTaxonomy=c(rep("Dummy", 5))
F3_16sY1_taxa.f=rbind(F3_16sY1_taxa.f,F3_16sY1_dummyTaxonomy)
rownames(F3_16sY1_taxa.f)[nrow(F3_16sY1_taxa.f)]="OTU_Filt"
rownames(F3_16sY1_otus.f)[nrow(F3_16sY1_otus.f)]="0TU_Filt"
tail(F3_16sY1_otus.f)


#Create Phyloseq object
phyloseq_spiec.easi_F1_16sY1 = phyloseq(otu_table(F1_16sY1_otus.f, taxa_are_rows = TRUE), tax_table(F1_16sY1_taxa.f))
phyloseq_spiec.easi_F2_16sY1 = phyloseq(otu_table(F2_16sY1_otus.f, taxa_are_rows = TRUE), tax_table(F2_16sY1_taxa.f))
phyloseq_spiec.easi_F3_16sY1 = phyloseq(otu_table(F3_16sY1_otus.f, taxa_are_rows = TRUE), tax_table(F3_16sY1_taxa.f))


#Use SPIEC-EASI to create the covariance matrix -Using pulsar with stars. 
spiec.out_F1_16sY1=spiec.easi(phyloseq_spiec.easi_F1_16sY1, method="mb",icov.select.params=list(rep.num=20)) 
spiec.out_F2_16sY1=spiec.easi(phyloseq_spiec.easi_F2_16sY1, method="mb",icov.select.params=list(rep.num=20))  
spiec.out_F3_16sY1=spiec.easi(phyloseq_spiec.easi_F3_16sY1, method="mb",icov.select.params=list(rep.num=20))  



#Network Analysis
betaMat_F1_16sY1=as.matrix(symBeta(getOptBeta(spiec.out_F1_16sY1)))
betaMat_F2_16sY1=as.matrix(symBeta(getOptBeta(spiec.out_F2_16sY1)))
betaMat_F3_16sY1=as.matrix(symBeta(getOptBeta(spiec.out_F3_16sY1)))

#Calculate Edge numbers in SPIEC-EASI
positive_F1_16sY1=length(betaMat_F1_16sY1[betaMat_F1_16sY1>0])/2
negative_F1_16sY1=length(betaMat_F1_16sY1[betaMat_F1_16sY1<0])/2
total_F1_16sY1=length(betaMat_F1_16sY1[betaMat_F1_16sY1!=0])/2

positive_F2_16sY1=length(betaMat_F2_16sY1[betaMat_F2_16sY1>0])/2
negative_F2_16sY1=length(betaMat_F2_16sY1[betaMat_F2_16sY1<0])/2
total_F2_16sY1=length(betaMat_F2_16sY1[betaMat_F2_16sY1!=0])/2

positive_F3_16sY1=length(betaMat_F3_16sY1[betaMat_F3_16sY1>0])/2
negative_F3_16sY1=length(betaMat_F3_16sY1[betaMat_F3_16sY1<0])/2
total_F3_16sY1=length(betaMat_F3_16sY1[betaMat_F3_16sY1!=0])/2

#Colors for the edges
otu.ids_F1_16sY1<-colnames(spiec.out_F1_16sY1[[1]]$data)
otu.ids_F2_16sY1<-colnames(spiec.out_F2_16sY1[[1]]$data)
otu.ids_F3_16sY1<-colnames(spiec.out_F3_16sY1[[1]]$data)

edge_cols_F1_16sY1 <-ifelse(betaMat_F1_16sY1>0, 'forestgreen', 'red')[upper.tri(betaMat_F1_16sY1) & betaMat_F1_16sY1!=0]
edge_cols_F2_16sY1 <-ifelse(betaMat_F2_16sY1>0, 'forestgreen', 'red')[upper.tri(betaMat_F2_16sY1) & betaMat_F2_16sY1!=0]
edge_cols_F3_16sY1 <-ifelse(betaMat_F3_16sY1>0, 'forestgreen', 'red')[upper.tri(betaMat_F3_16sY1) & betaMat_F3_16sY1!=0]

ig2.mb_F1_16sY1 <- adj2igraph(getRefit(spiec.out_F1_16sY1),  rmEmptyNodes=TRUE,
                     vertex.attr=list(name=taxa_names(phyloseq_spiec.easi_F1_16sY1)),
                     edge.attr=list(color= edge_cols_F1_16sY1))

ig2.mb_F2_16sY1 <- adj2igraph(getRefit(spiec.out_F2_16sY1),  rmEmptyNodes=TRUE,
                              vertex.attr=list(name=taxa_names(phyloseq_spiec.easi_F2_16sY1)),
                              edge.attr=list(color= edge_cols_F2_16sY1))

ig2.mb_F3_16sY1 <- adj2igraph(getRefit(spiec.out_F3_16sY1),  rmEmptyNodes=TRUE,
                              vertex.attr=list(name=taxa_names(phyloseq_spiec.easi_F3_16sY1)),
                              edge.attr=list(color= edge_cols_F3_16sY1))


#How many nodes connected at specific rank
library(intergraph)
library(network)

getrank_F1_16sY1 ="Family"
nb_nodes_F1_16sY1 <- vcount(ig2.mb_F1_16sY1)
tax_table(phyloseq_spiec.easi_F1_16sY1) <- tax_table(phyloseq_spiec.easi_F1_16sY1)[,getrank]
otu_ids_F1_16sY1 <- V(ig2.mb_F1_16sY1)$name
idx_F1_16sY1 <- which(row.names(tax_table(phyloseq_spiec.easi_F1_16sY1)) %in% otu_ids_F1_16sY1)
taxanames_F1_16sY1 <- as.character(rownames(tax_table(phyloseq_spiec.easi_F1_16sY1)))[idx_F1_16sY1]
ig2_F1_16sY1 <- asNetwork(ig2.mb_F1_16sY1)
network.vertex.names(ig2_F1_16sY1) <- taxanames_F1_16sY1
net_F1_16sY1 <- ig2_F1_16sY1
net_F1_16sY1 %v% getrank_F1_16sY1  = as.character(taxanames_F1_16sY1)

getrank_F2_16sY1 ="Family"
nb_nodes_F2_16sY1 <- vcount(ig2.mb_F2_16sY1)
tax_table(phyloseq_spiec.easi_F2_16sY1) <- tax_table(phyloseq_spiec.easi_F2_16sY1)[,getrank]
otu_ids_F2_16sY1 <- V(ig2.mb_F2_16sY1)$name
idx_F2_16sY1 <- which(row.names(tax_table(phyloseq_spiec.easi_F2_16sY1)) %in% otu_ids_F2_16sY1)
taxanames_F2_16sY1 <- as.character(rownames(tax_table(phyloseq_spiec.easi_F2_16sY1)))[idx_F2_16sY1]
ig2_F2_16sY1 <- asNetwork(ig2.mb_F2_16sY1)
network.vertex.names(ig2_F2_16sY1) <- taxanames_F2_16sY1
net_F2_16sY1 <- ig2_F2_16sY1
net_F2_16sY1 %v% getrank_F2_16sY1  = as.character(taxanames_F2_16sY1)

getrank_F3_16sY1 ="Family"
nb_nodes_F3_16sY1 <- vcount(ig2.mb_F3_16sY1)
tax_table(phyloseq_spiec.easi_F3_16sY1) <- tax_table(phyloseq_spiec.easi_F3_16sY1)[,getrank]
otu_ids_F3_16sY1 <- V(ig2.mb_F3_16sY1)$name
idx_F3_16sY1 <- which(row.names(tax_table(phyloseq_spiec.easi_F3_16sY1)) %in% otu_ids_F3_16sY1)
taxanames_F3_16sY1 <- as.character(rownames(tax_table(phyloseq_spiec.easi_F3_16sY1)))[idx_F3_16sY1]
ig2_F3_16sY1 <- asNetwork(ig2.mb_F3_16sY1)
network.vertex.names(ig2_F3_16sY1) <- taxanames_F3_16sY1
net_F3_16sY1 <- ig2_F3_16sY1
net_F3_16sY1 %v% getrank_F3_16sY1  = as.character(taxanames_F3_16sY1)


#Plot the network
library(GGally)

network_plot_F1_16sY1 <-  ggnet2(net_F1_16sY1,
             color = getrank_F1_16sY1 ,
             alpha = 0.75,
             size = 6, 
             edge.size=1,
             edge.color="color",
             edge.alpha = 0.5,
             label = TRUE, 
             label.size = 4)+theme(legend.position = 'none')+ggtitle("Bacteria Y1 - F1")
network_plot_F1_16sY1

network_plot_F2_16sY1 <-  ggnet2(net_F2_16sY1,
                                 color = getrank_F2_16sY1 ,
                                 alpha = 0.75,
                                 size = 6, 
                                 edge.size=1,
                                 edge.color="color",
                                 edge.alpha = 0.5,
                                 label = TRUE, 
                                 label.size = 4)+theme(legend.position = 'none')+ggtitle("Bacteria Y1 - F2")
network_plot_F2_16sY1

network_plot_F3_16sY1 <-  ggnet2(net_F3_16sY1,
                                 color = getrank_F3_16sY1 ,
                                 alpha = 0.75,
                                 size = 6, 
                                 edge.size=1,
                                 edge.color="color",
                                 edge.alpha = 0.5,
                                 label = TRUE, 
                                 label.size = 4)+theme(legend.position = 'none')+ggtitle("Bacteria Y1 - F3")
network_plot_F3_16sY1

#Cluster the network and list the Taxa present in each cluster
spiec.graph_F1_16sY1=adj2igraph(getRefit(spiec.out_F1_16sY1), vertex.attr=list(name=taxa_names(phyloseq_spiec.easi_F1_16sY1)))
plot_network(spiec.graph_F1_16sY1, phyloseq_spiec.easi_F1_16sY1, type='taxa', color="Family", label=NULL)


clusters_F1_16sY1=cluster_fast_greedy(spiec.graph_F1_16sY1)
clusterOneIndices_F1_16sY1=which(clusters_F1_16sY1$membership==1)
clusterOneOtus_F1_16sY1=clusters_F1_16sY1$names[clusterOneIndices_F1_16sY1]

sort(table(getTaxonomy(clusterOneOtus_F1_16sY1,F1_16sY1_taxa.f,useRownames = TRUElevel="Family",)),decreasing = TRUE)

clusterTwoIndices=which(clusters$membership==2)
clusterTwoOtus=clusters$names[clusterTwoIndices]
















#### Differential abundance by facility at Family level ####

#Make Phyloseq with count data
OTU_16sY1__count <- otu_table(otus_16sY1_, taxa_are_rows = TRUE)
phyloseq_16sY1__count = phyloseq(OTU_16sY1__count, TAX_16s, META_16sY1)
TREE_16sY1__count = rtree(ntaxa(phyloseq_16sY1_), rooted=TRUE, tip.label = taxa_names(phyloseq_16sY1_))  
phyloseq16sY1__count <- phyloseq(OTU_16sY1__count, TAX_16s, TREE_16sY1__count, META_16sY1)

OTU_16sY2__count <- otu_table(otus_16sY2_, taxa_are_rows = TRUE)
phyloseq_16sY2__count = phyloseq(OTU_16sY2__count, TAX_16s, META_16sY2)
TREE_16sY2__count = rtree(ntaxa(phyloseq_16sY2_), rooted=TRUE, tip.label = taxa_names(phyloseq_16sY2_))  
phyloseq16sY2__count <- phyloseq(OTU_16sY2__count, TAX_16s, TREE_16sY2__count, META_16sY2)


#Make phyloseq collapse Phyloseq to family level
phyloseq_family_16sY1__count <- phyloseq16sY1__count %>%  
  tax_glom(taxrank = "Family")

phyloseq_family_16sY2__count <- phyloseq16sY2__count %>%  
  tax_glom(taxrank = "Family")


# Subset Phyloseq with read counts at family level for each facility
F1_16sY1_count <- subset_samples(phyloseq_family_16sY1__count, Facility == "F1") 
F2_16sY1_count <- subset_samples(phyloseq_family_16sY1__count, Facility == "F2") 
F3_16sY1_count <- subset_samples(phyloseq_family_16sY1__count, Facility == "F3") 

F1_16sY2_count <- subset_samples(phyloseq_family_16sY2__count, Facility == "F1") 
F2_16sY2_count <- subset_samples(phyloseq_family_16sY2__count, Facility == "F2") 
F3_16sY2_count <- subset_samples(phyloseq_family_16sY2__count, Facility == "F3") 

# F1_ITSY1_count <- subset_samples(phyloseqITSY1_family, Facility == "F1") 
# F2_ITSY1_count <- subset_samples(phyloseqITSY1_family, Facility == "F2") 
# F3_ITSY1_count <- subset_samples(phyloseqITSY1_family, Facility == "F3") 
# 
# F1_ITSY2_count <- subset_samples(phyloseqITSY2_family, Facility == "F1") 
# F2_ITSY2_count <- subset_samples(phyloseqITSY2_family, Facility == "F2") 
# F3_ITSY2_count <- subset_samples(phyloseqITSY2_family, Facility == "F3") 

#Take otu table for each facility from Phyloseq
family_F116sY1_otu<-as.data.frame(as(otu_table(F1_16sY1_count), "matrix"))
family_F216sY1_otu<-as.data.frame(as(otu_table(F2_16sY1_count), "matrix"))
family_F316sY1_otu<-as.data.frame(as(otu_table(F3_16sY1_count), "matrix"))

family_F116sY2_otu<-as.data.frame(as(otu_table(F1_16sY2_count), "matrix"))
family_F216sY2_otu<-as.data.frame(as(otu_table(F2_16sY2_count), "matrix"))
family_F316sY2_otu<-as.data.frame(as(otu_table(F3_16sY2_count), "matrix"))

# family_F1ITSY1_otu<-as.data.frame(as(otu_table(F1_ITSY1_count), "matrix"))
# family_F2ITSY1_otu<-as.data.frame(as(otu_table(F2_ITSY1_count), "matrix"))
# family_F3ITSY1_otu<-as.data.frame(as(otu_table(F3_ITSY1_count), "matrix"))
# 
# family_F1ITSY2_otu<-as.data.frame(as(otu_table(F1_ITSY2_count), "matrix"))
# family_F2ITSY2_otu<-as.data.frame(as(otu_table(F2_ITSY2_count), "matrix"))
# family_F3ITSY2_otu<-as.data.frame(as(otu_table(F3_ITSY2_count), "matrix"))

#Merge OTU tables for each pair of facilities
family_otu_16sY1_F1F2<-as.data.frame(bind_cols(family_F116sY1_otu, family_F216sY1_otu)) #Otu table for aldex2 needs to have otus in rows
family_otu_16sY1_F1F3<-as.data.frame(bind_cols(family_F116sY1_otu, family_F316sY1_otu)) #Otu table for aldex2 needs to have otus in rows
family_otu_16sY1_F2F3<-as.data.frame(bind_cols(family_F216sY1_otu, family_F316sY1_otu)) #Otu table for aldex2 needs to have otus in rows

family_otu_16sY2_F1F2<-as.data.frame(bind_cols(family_F116sY2_otu, family_F216sY2_otu)) #Otu table for aldex2 needs to have otus in rows
family_otu_16sY2_F1F3<-as.data.frame(bind_cols(family_F116sY2_otu, family_F316sY2_otu)) #Otu table for aldex2 needs to have otus in rows
family_otu_16sY2_F2F3<-as.data.frame(bind_cols(family_F216sY2_otu, family_F316sY2_otu)) #Otu table for aldex2 needs to have otus in rows

# family_otu_ITSY1_F1F2<-as.data.frame(bind_cols(family_F1ITSY1_otu, family_F2ITSY1_otu)) #Otu table for aldex2 needs to have otus in rows
# family_otu_ITSY1_F1F3<-as.data.frame(bind_cols(family_F1ITSY1_otu, family_F3ITSY1_otu)) #Otu table for aldex2 needs to have otus in rows
# family_otu_ITSY1_F2F3<-as.data.frame(bind_cols(family_F2ITSY1_otu, family_F3ITSY1_otu)) #Otu table for aldex2 needs to have otus in rows
# 
# family_otu_ITSY2_F1F2<-as.data.frame(bind_cols(family_F1ITSY2_otu, family_F2ITSY2_otu)) #Otu table for aldex2 needs to have otus in rows
# family_otu_ITSY2_F1F3<-as.data.frame(bind_cols(family_F1ITSY2_otu, family_F3ITSY2_otu)) #Otu table for aldex2 needs to have otus in rows
# family_otu_ITSY2_F2F3<-as.data.frame(bind_cols(family_F2ITSY2_otu, family_F3ITSY2_otu)) #Otu table for aldex2 needs to have otus in rows

#Subset metadata by pairs of facilities. Aldex needs a category with the same amount of samples as input
metadata_F116sY1<-metadata_16sY1[ which(metadata_16sY1$Facility =="F1"),]
metadata_F216sY1<-metadata_16sY1[ which(metadata_16sY1$Facility =="F2"),]
metadata_F316sY1<-metadata_16sY1[ which(metadata_16sY1$Facility =="F3"),]

metadata_F116sY2<-metadata_16sY2[ which(metadata_16sY2$Facility =="F1"),]
metadata_F216sY2<-metadata_16sY2[ which(metadata_16sY2$Facility =="F2"),]
metadata_F316sY2<-metadata_16sY2[ which(metadata_16sY2$Facility =="F3"),]

# metadata_F1ITSY1<-metadata_ITSY1[ which(metadata_ITSY1$Facility =="F1"),]
# metadata_F2ITSY1<-metadata_ITSY1[ which(metadata_ITSY1$Facility =="F2"),]
# metadata_F3ITSY1<-metadata_ITSY1[ which(metadata_ITSY1$Facility =="F3"),]
# 
# metadata_F1ITSY2<-metadata_ITSY2[ which(metadata_ITSY2$Facility =="F1"),]
# metadata_F2ITSY2<-metadata_ITSY2[ which(metadata_ITSY2$Facility =="F2"),]
# metadata_F3ITSY2<-metadata_ITSY2[ which(metadata_ITSY2$Facility =="F3"),]

#Change Facility to character vector - as.factor affects aldex.effect function
metadata_F116sY1$Facility <- as.character(metadata_F116sY1$Facility)
metadata_F216sY1$Facility <- as.character(metadata_F216sY1$Facility)
metadata_F316sY1$Facility <- as.character(metadata_F316sY1$Facility)

metadata_F116sY2$Facility <- as.character(metadata_F116sY2$Facility)
metadata_F216sY2$Facility <- as.character(metadata_F216sY2$Facility)
metadata_F316sY2$Facility <- as.character(metadata_F316sY2$Facility)

# metadata_F1ITSY1$Facility <- as.character(metadata_F1ITSY1$Facility)
# metadata_F2ITSY1$Facility <- as.character(metadata_F2ITSY1$Facility)
# metadata_F3ITSY1$Facility <- as.character(metadata_F3ITSY1$Facility)
# 
# metadata_F1ITSY2$Facility <- as.character(metadata_F1ITSY2$Facility)
# metadata_F2ITSY2$Facility <- as.character(metadata_F2ITSY2$Facility)
# metadata_F3ITSY2$Facility <- as.character(metadata_F3ITSY2$Facility)

#Bind metadata by pair of facilities
family_metadata_16sY1_F1F2<-bind_rows(metadata_F116sY1, metadata_F216sY1) 
family_metadata_16sY1_F1F3<-bind_rows(metadata_F116sY1, metadata_F316sY1)
family_metadata_16sY1_F2F3<-bind_rows(metadata_F216sY1, metadata_F316sY1)

family_metadata_16sY2_F1F2<-bind_rows(metadata_F116sY2, metadata_F216sY2) 
family_metadata_16sY2_F1F3<-bind_rows(metadata_F116sY2, metadata_F316sY2)
family_metadata_16sY2_F2F3<-bind_rows(metadata_F216sY2, metadata_F316sY2)

# family_metadata_ITSY1_F1F2<-bind_rows(metadata_F1ITSY1, metadata_F2ITSY1) 
# family_metadata_ITSY1_F1F3<-bind_rows(metadata_F1ITSY1, metadata_F3ITSY1)
# family_metadata_ITSY1_F2F3<-bind_rows(metadata_F2ITSY1, metadata_F3ITSY1)
# 
# family_metadata_ITSY2_F1F2<-bind_rows(metadata_F1ITSY2, metadata_F2ITSY2) 
# family_metadata_ITSY2_F1F3<-bind_rows(metadata_F1ITSY2, metadata_F3ITSY2)
# family_metadata_ITSY2_F2F3<-bind_rows(metadata_F2ITSY2, metadata_F3ITSY2)

#Differential abundance analysis 
#Based on Gloor et al 2016. Annals of Epidemiology 26:322-329.
#Supplemental material page 19.

#Generate 128 Dirichlet distributed Monte Carlo instances and center-log ratio transform them
#OTU table needs to be as OTU in rows. Only two conditions can be compared at a time.

#Compare between samples that were positive and negative
#At family level
Aldex_family_16sY1_F1F2.clr<-aldex.clr(family_otu_16sY1_F1F2, mc.samples = 128, conds = family_metadata_16sY1_F1F2$Facility)
Aldex_family_16sY1_F1F3.clr<-aldex.clr(family_otu_16sY1_F1F3, mc.samples = 128, conds = family_metadata_16sY1_F1F3$Facility)
Aldex_family_16sY1_F2F3.clr<-aldex.clr(family_otu_16sY1_F2F3, mc.samples = 128, conds = family_metadata_16sY1_F2F3$Facility)

Aldex_family_16sY2_F1F2.clr<-aldex.clr(family_otu_16sY2_F1F2, mc.samples = 128, conds = family_metadata_16sY2_F1F2$Facility)
Aldex_family_16sY2_F1F3.clr<-aldex.clr(family_otu_16sY2_F1F3, mc.samples = 128, conds = family_metadata_16sY2_F1F3$Facility)
Aldex_family_16sY2_F2F3.clr<-aldex.clr(family_otu_16sY2_F2F3, mc.samples = 128, conds = family_metadata_16sY2_F2F3$Facility)

# Aldex_family_ITSY1_F1F2.clr<-aldex.clr(family_otu_ITSY1_F1F2, mc.samples = 128, conds = family_metadata_ITSY1_F1F2$Facility)
# Aldex_family_ITSY1_F1F3.clr<-aldex.clr(family_otu_ITSY1_F1F3, mc.samples = 128, conds = family_metadata_ITSY1_F1F3$Facility)
# Aldex_family_ITSY1_F2F3.clr<-aldex.clr(family_otu_ITSY1_F2F3, mc.samples = 128, conds = family_metadata_ITSY1_F2F3$Facility)
# 
# Aldex_family_ITSY2_F1F2.clr<-aldex.clr(family_otu_ITSY2_F1F2, mc.samples = 128, conds = family_metadata_ITSY2_F1F2$Facility)
# Aldex_family_ITSY2_F1F3.clr<-aldex.clr(family_otu_ITSY2_F1F3, mc.samples = 128, conds = family_metadata_ITSY2_F1F3$Facility)
# Aldex_family_ITSY2_F2F3.clr<-aldex.clr(family_otu_ITSY2_F2F3, mc.samples = 128, conds = family_metadata_ITSY2_F2F3$Facility)


#Calculate the expected effect size
#At family level
Aldex_family_16sY1_F1F2.e<-aldex.effect(Aldex_family_16sY1_F1F2.clr)
Aldex_family_16sY1_F1F3.e<-aldex.effect(Aldex_family_16sY1_F1F3.clr)
Aldex_family_16sY1_F2F3.e<-aldex.effect(Aldex_family_16sY1_F2F3.clr)

Aldex_family_16sY2_F1F2.e<-aldex.effect(Aldex_family_16sY2_F1F2.clr)
Aldex_family_16sY2_F1F3.e<-aldex.effect(Aldex_family_16sY2_F1F3.clr)
Aldex_family_16sY2_F2F3.e<-aldex.effect(Aldex_family_16sY2_F2F3.clr)

# Aldex_family_ITSY1_F1F2.e<-aldex.effect(Aldex_family_ITSY1_F1F2.clr)
# Aldex_family_ITSY1_F1F3.e<-aldex.effect(Aldex_family_ITSY1_F1F3.clr)
# Aldex_family_ITSY1_F2F3.e<-aldex.effect(Aldex_family_ITSY1_F2F3.clr)
# 
# Aldex_family_ITSY2_F1F2.e<-aldex.effect(Aldex_family_ITSY2_F1F2.clr)
# Aldex_family_ITSY2_F1F3.e<-aldex.effect(Aldex_family_ITSY2_F1F3.clr)
# Aldex_family_ITSY2_F2F3.e<-aldex.effect(Aldex_family_ITSY2_F2F3.clr)

#Generate the P and Benjamini-Hochberg corrected P values.
#At family level
Aldex_family_16sY1_F1F2.t<-aldex.ttest(Aldex_family_16sY1_F1F2.clr)
Aldex_family_16sY1_F1F3.t<-aldex.ttest(Aldex_family_16sY1_F1F3.clr)
Aldex_family_16sY1_F2F3.t<-aldex.ttest(Aldex_family_16sY1_F2F3.clr)

Aldex_family_16sY2_F1F2.t<-aldex.ttest(Aldex_family_16sY2_F1F2.clr)
Aldex_family_16sY2_F1F3.t<-aldex.ttest(Aldex_family_16sY2_F1F3.clr)
Aldex_family_16sY2_F2F3.t<-aldex.ttest(Aldex_family_16sY2_F2F3.clr)

# Aldex_family_ITSY1_F1F2.t<-aldex.ttest(Aldex_family_ITSY1_F1F2.clr)
# Aldex_family_ITSY1_F1F3.t<-aldex.ttest(Aldex_family_ITSY1_F1F3.clr)
# Aldex_family_ITSY1_F2F3.t<-aldex.ttest(Aldex_family_ITSY1_F2F3.clr)
# 
# Aldex_family_ITSY2_F1F2.t<-aldex.ttest(Aldex_family_ITSY2_F1F2.clr)
# Aldex_family_ITSY2_F1F3.t<-aldex.ttest(Aldex_family_ITSY2_F1F3.clr)
# Aldex_family_ITSY2_F2F3.t<-aldex.ttest(Aldex_family_ITSY2_F2F3.clr)

#Merge data frames
#At family level
Aldex_family_16sY1_F1F2.all<-data.frame(Aldex_family_16sY1_F1F2.e,Aldex_family_16sY1_F1F2.t)
Aldex_family_16sY1_F1F3.all<-data.frame(Aldex_family_16sY1_F1F3.e,Aldex_family_16sY1_F1F3.t)
Aldex_family_16sY1_F2F3.all<-data.frame(Aldex_family_16sY1_F2F3.e,Aldex_family_16sY1_F2F3.t)

Aldex_family_16sY2_F1F2.all<-data.frame(Aldex_family_16sY2_F1F2.e,Aldex_family_16sY2_F1F2.t)
Aldex_family_16sY2_F1F3.all<-data.frame(Aldex_family_16sY2_F1F3.e,Aldex_family_16sY2_F1F3.t)
Aldex_family_16sY2_F2F3.all<-data.frame(Aldex_family_16sY2_F2F3.e,Aldex_family_16sY2_F2F3.t)

# Aldex_family_ITSY1_F1F2.all<-data.frame(Aldex_family_ITSY1_F1F2.e,Aldex_family_ITSY1_F1F2.t)
# Aldex_family_ITSY1_F1F3.all<-data.frame(Aldex_family_ITSY1_F1F3.e,Aldex_family_ITSY1_F1F3.t)
# Aldex_family_ITSY1_F2F3.all<-data.frame(Aldex_family_ITSY1_F2F3.e,Aldex_family_ITSY1_F2F3.t)
# 
# Aldex_family_ITSY2_F1F2.all<-data.frame(Aldex_family_ITSY2_F1F2.e,Aldex_family_ITSY2_F1F2.t)
# Aldex_family_ITSY2_F1F3.all<-data.frame(Aldex_family_ITSY2_F1F3.e,Aldex_family_ITSY2_F1F3.t)
# Aldex_family_ITSY2_F2F3.all<-data.frame(Aldex_family_ITSY2_F2F3.e,Aldex_family_ITSY2_F2F3.t)

#Determine which corrected values fall below a threshold
#At family level
Aldex_family_16sY1_F1F2.sig<-which(Aldex_family_16sY1_F1F2.all$wi.eBH <=0.05)
Aldex_family_16sY1_F1F3.sig<-which(Aldex_family_16sY1_F1F3.all$wi.eBH <=0.05)
Aldex_family_16sY1_F2F3.sig<-which(Aldex_family_16sY1_F2F3.all$wi.eBH <=0.05)

Aldex_family_16sY2_F1F2.sig<-which(Aldex_family_16sY2_F1F2.all$wi.eBH <=0.05)
Aldex_family_16sY2_F1F3.sig<-which(Aldex_family_16sY2_F1F3.all$wi.eBH <=0.05)
Aldex_family_16sY2_F2F3.sig<-which(Aldex_family_16sY2_F2F3.all$wi.eBH <=0.05)

# Aldex_family_ITSY1_F1F2.sig<-which(Aldex_family_ITSY1_F1F2.all$wi.eBH <=0.05)
# Aldex_family_ITSY1_F1F3.sig<-which(Aldex_family_ITSY1_F1F3.all$wi.eBH <=0.05)
# Aldex_family_ITSY1_F2F3.sig<-which(Aldex_family_ITSY1_F2F3.all$wi.eBH <=0.05)
# 
# Aldex_family_ITSY2_F1F2.sig<-which(Aldex_family_ITSY2_F1F2.all$wi.eBH <=0.05)
# Aldex_family_ITSY2_F1F3.sig<-which(Aldex_family_ITSY2_F1F3.all$wi.eBH <=0.05)
# Aldex_family_ITSY2_F2F3.sig<-which(Aldex_family_ITSY2_F2F3.all$wi.eBH <=0.05)

#Plot the results - Effect plots described by the documentation

png("Aldex -Bacteria by Facility - effect plot.png", width = 10, height = 10, units = 'in', res=600)

#16sY1
par(mar=c(4,6,3,4))
par(mfrow=c(4,3))
plot(Aldex_family_16sY1_F1F2.all$diff.win, Aldex_family_16sY1_F1F2.all$diff.btw, pch=19, cex=0.2,
     col='grey', xlab='within group dispersion', ylab='between group difference', main="Bacteria Y1 - F1 v F2")
points(Aldex_family_16sY1_F1F2.all$diff.win[Aldex_family_16sY1_F1F2.sig], 
       Aldex_family_16sY1_F1F2.all$diff.btw[Aldex_family_16sY1_F1F2.sig], pch=19, cex=0.2,col='red')#color significant ones in reg
abline(0,1, col='grey', lty=1)#Add approximate effect size lines
abline(0,-1, col='grey', lty=1)
abline(0,0.5, col='grey', lty=2)
abline(0,-0.5, col='grey', lty=2)

plot(Aldex_family_16sY1_F1F3.all$diff.win, Aldex_family_16sY1_F1F3.all$diff.btw, pch=19, cex=0.2,
     col='grey', xlab="within group dispersion", ylab='between group difference', main="Bacteria Y1 - F1 v F3")
points(Aldex_family_16sY1_F1F3.all$diff.win[Aldex_family_16sY1_F1F3.sig], 
       Aldex_family_16sY1_F1F3.all$diff.btw[Aldex_family_16sY1_F1F3.sig], pch=19, cex=0.2,col='red')#color significant ones in reg
abline(0,1, col='grey', lty=1)#Add approximate effect size lines
abline(0,-1, col='grey', lty=1)
abline(0,0.5, col='grey', lty=2)
abline(0,-0.5, col='grey', lty=2)

plot(Aldex_family_16sY1_F2F3.all$diff.win, Aldex_family_16sY1_F2F3.all$diff.btw, pch=19, cex=0.2,
     col='grey', xlab="within group dispersion", ylab='between group difference', main="Bacteria Y1 - F2 v F3")
points(Aldex_family_16sY1_F2F3.all$diff.win[Aldex_family_16sY1_F2F3.sig], 
       Aldex_family_16sY1_F2F3.all$diff.btw[Aldex_family_16sY1_F2F3.sig], pch=19, cex=0.2,col='red')#color significant ones in reg
abline(0,1, col='grey', lty=1)#Add approximate effect size lines
abline(0,-1, col='grey', lty=1)
abline(0,0.5, col='grey', lty=2)
abline(0,-0.5, col='grey', lty=2)

#16sY2
plot(Aldex_family_16sY2_F1F2.all$diff.win, Aldex_family_16sY2_F1F2.all$diff.btw, pch=19, cex=0.2,
     col='grey', xlab="within group dispersion", ylab='between group difference', main="Bacteria Y2 - F1 v F2")
points(Aldex_family_16sY2_F1F2.all$diff.win[Aldex_family_16sY2_F1F2.sig], 
       Aldex_family_16sY2_F1F2.all$diff.btw[Aldex_family_16sY2_F1F2.sig], pch=19, cex=0.2,col='red')#color significant ones in reg
abline(0,1, col='grey', lty=1)#Add approximate effect size lines
abline(0,-1, col='grey', lty=1)
abline(0,0.5, col='grey', lty=2)
abline(0,-0.5, col='grey', lty=2)

plot(Aldex_family_16sY2_F1F3.all$diff.win, Aldex_family_16sY2_F1F3.all$diff.btw, pch=19, cex=0.2,
     col='grey', xlab="within group dispersion", ylab='between group difference', main="Bacteria Y2 - F1 v F3")
points(Aldex_family_16sY2_F1F3.all$diff.win[Aldex_family_16sY2_F1F3.sig], 
       Aldex_family_16sY2_F1F3.all$diff.btw[Aldex_family_16sY2_F1F3.sig], pch=19, cex=0.2,col='red')#color significant ones in reg
abline(0,1, col='grey', lty=1)#Add approximate effect size lines
abline(0,-1, col='grey', lty=1)
abline(0,0.5, col='grey', lty=2)
abline(0,-0.5, col='grey', lty=2)

plot(Aldex_family_16sY2_F2F3.all$diff.win, Aldex_family_16sY2_F2F3.all$diff.btw, pch=19, cex=0.2,
     col='grey', xlab="within group dispersion", ylab='between group difference', main="Bacteria Y2 - F2 v F3")
points(Aldex_family_16sY2_F2F3.all$diff.win[Aldex_family_16sY2_F2F3.sig], 
       Aldex_family_16sY2_F2F3.all$diff.btw[Aldex_family_16sY2_F2F3.sig], pch=19, cex=0.2,col='red')#color significant ones in reg
abline(0,1, col='grey', lty=1)#Add approximate effect size lines
abline(0,-1, col='grey', lty=1)
abline(0,0.5, col='grey', lty=2)
abline(0,-0.5, col='grey', lty=2)

# #ITSY1
# plot(Aldex_family_ITSY1_F1F2.all$diff.win, Aldex_family_ITSY1_F1F2.all$diff.btw, pch=19, cex=0.2,
#      col='grey', xlab="within group dispersion", ylab='between group difference', main="Fungi Y1 - F1 v F2")
# points(Aldex_family_ITSY1_F1F2.all$diff.win[Aldex_family_ITSY1_F1F2.sig], 
#        Aldex_family_ITSY1_F1F2.all$diff.btw[Aldex_family_ITSY1_F1F2.sig], pch=19, cex=0.2,col='red')#color significant ones in reg
# abline(0,1, col='grey', lty=1)#Add approximate effect size lines
# abline(0,-1, col='grey', lty=1)
# abline(0,0.5, col='grey', lty=2)
# abline(0,-0.5, col='grey', lty=2)
# 
# plot(Aldex_family_ITSY1_F1F3.all$diff.win, Aldex_family_ITSY1_F1F3.all$diff.btw, pch=19, cex=0.2,
#      col='grey', xlab="within group dispersion", ylab='between group difference', main="Fungi Y1 - F1 v F3")
# points(Aldex_family_ITSY1_F1F3.all$diff.win[Aldex_family_ITSY1_F1F3.sig], 
#        Aldex_family_ITSY1_F1F3.all$diff.btw[Aldex_family_ITSY1_F1F3.sig], pch=19, cex=0.2,col='red')#color significant ones in reg
# abline(0,1, col='grey', lty=1)#Add approximate effect size lines
# abline(0,-1, col='grey', lty=1)
# abline(0,0.5, col='grey', lty=2)
# abline(0,-0.5, col='grey', lty=2)
# 
# plot(Aldex_family_ITSY1_F2F3.all$diff.win, Aldex_family_ITSY1_F2F3.all$diff.btw, pch=19, cex=0.2,
#      col='grey', xlab="within group dispersion", ylab='between group difference', main="Fungi Y1 - F2 v F3")
# points(Aldex_family_ITSY1_F2F3.all$diff.win[Aldex_family_ITSY1_F2F3.sig], 
#        Aldex_family_ITSY1_F2F3.all$diff.btw[Aldex_family_ITSY1_F2F3.sig], pch=19, cex=0.2,col='red')#color significant ones in reg
# abline(0,1, col='grey', lty=1)#Add approximate effect size lines
# abline(0,-1, col='grey', lty=1)
# abline(0,0.5, col='grey', lty=2)
# abline(0,-0.5, col='grey', lty=2)
# 
# #ITSY2
# plot(Aldex_family_ITSY2_F1F2.all$diff.win, Aldex_family_ITSY2_F1F2.all$diff.btw, pch=19, cex=0.2,
#      col='grey', xlab="within group dispersion", ylab='between group difference', main="Fungi Y2 - F1 v F2")
# points(Aldex_family_ITSY2_F1F2.all$diff.win[Aldex_family_ITSY2_F1F2.sig], 
#        Aldex_family_ITSY2_F1F2.all$diff.btw[Aldex_family_ITSY2_F1F2.sig], pch=19, cex=0.2,col='red')#color significant ones in reg
# abline(0,1, col='grey', lty=1)#Add approximate effect size lines
# abline(0,-1, col='grey', lty=1)
# abline(0,0.5, col='grey', lty=2)
# abline(0,-0.5, col='grey', lty=2)
# 
# plot(Aldex_family_ITSY2_F1F3.all$diff.win, Aldex_family_ITSY2_F1F3.all$diff.btw, pch=19, cex=0.2,
#      col='grey', xlab="within group dispersion", ylab='between group difference', main="Fungi Y2 - F1 v F3")
# points(Aldex_family_ITSY2_F1F3.all$diff.win[Aldex_family_ITSY2_F1F3.sig], 
#        Aldex_family_ITSY2_F1F3.all$diff.btw[Aldex_family_ITSY2_F1F3.sig], pch=19, cex=0.2,col='red')#color significant ones in reg
# abline(0,1, col='grey', lty=1)#Add approximate effect size lines
# abline(0,-1, col='grey', lty=1)
# abline(0,0.5, col='grey', lty=2)
# abline(0,-0.5, col='grey', lty=2)
# 
# plot(Aldex_family_ITSY2_F2F3.all$diff.win, Aldex_family_ITSY2_F2F3.all$diff.btw, pch=19, cex=0.2,
#      col='grey', xlab="within group dispersion", ylab='between group difference', main="Fungi Y2 - F2 v F3")
# points(Aldex_family_ITSY2_F2F3.all$diff.win[Aldex_family_ITSY2_F2F3.sig], 
#        Aldex_family_ITSY2_F2F3.all$diff.btw[Aldex_family_ITSY2_F2F3.sig], pch=19, cex=0.2,col='red')#color significant ones in reg
# abline(0,1, col='grey', lty=1)#Add approximate effect size lines
# abline(0,-1, col='grey', lty=1)
# abline(0,0.5, col='grey', lty=2)
# abline(0,-0.5, col='grey', lty=2)

dev.off()


#Plots of significant families relative abundances
#Extract significant OTU data from ALDEx2 output
# 16sY1 
Aldex_family_16sY1_F1F2.sig.row<-rownames(Aldex_family_16sY1_F1F2.all)[which(Aldex_family_16sY1_F1F2.all$wi.eBH <=0.05)] #select the rownames of significant families
Aldex_family_16sY1_F1F2.sig.table<-subset(Aldex_family_16sY1_F1F2.all, rownames(Aldex_family_16sY1_F1F2.all) %in% Aldex_family_16sY1_F1F2.sig.row) #Subset significant families
# Aldex_family_16sY1_F1F2.sig.otu<-as.data.frame(subset((t(otu.n0.acomp_16sY1)), rownames(t(otu.n0.acomp_16sY1)) %in% Aldex_family_16sY1_F1F2.sig.row)) #Subset table with relative abundance per sample, only for significant fammilies
Aldex_family_16sY1_F1F2.sig.taxon<-subset(taxon_16s, rownames(taxon_16s) %in% Aldex_family_16sY1_F1F2.sig.row) #Subset taxon table to get the taxonomy of significant families
Aldex_family_16sY1_F1F2.sig.table.all<-bind_cols(Aldex_family_16sY1_F1F2.sig.taxon, Aldex_family_16sY1_F1F2.sig.table) #combine tables
Aldex_family_16sY1_F1F2.sig.table.all$Comparison<-rep("F1 v F2", 36)
Aldex_family_16sY1_F1F2.sig.table.all$Facility<-ifelse(Aldex_family_16sY1_F1F2.sig.table.all$effect<0, "F2", "F1") #Add direction of significance


Aldex_family_16sY1_F1F3.sig.row<-rownames(Aldex_family_16sY1_F1F3.all)[which(Aldex_family_16sY1_F1F3.all$wi.eBH <=0.05)] #select the rownames of significant families
Aldex_family_16sY1_F1F3.sig.table<-subset(Aldex_family_16sY1_F1F3.all, rownames(Aldex_family_16sY1_F1F3.all) %in% Aldex_family_16sY1_F1F3.sig.row) #Subset significant families
# Aldex_family_16sY1_F1F3.sig.otu<-as.data.frame(subset((t(otu.n0.acomp_16sY1)), rownames(t(otu.n0.acomp_16sY1)) %in% Aldex_family_16sY1_F1F3.sig.row)) #Subset table with relative abundance per sample, only for significant fammilies
Aldex_family_16sY1_F1F3.sig.taxon<-subset(taxon_16s, rownames(taxon_16s) %in% Aldex_family_16sY1_F1F3.sig.row) #Subset taxon table to get the taxonomy of significant families
Aldex_family_16sY1_F1F3.sig.table.all<-bind_cols(Aldex_family_16sY1_F1F3.sig.taxon, Aldex_family_16sY1_F1F3.sig.table) #combine tables
Aldex_family_16sY1_F1F3.sig.table.all$Comparison<-rep("F1 v F3", 40)
Aldex_family_16sY1_F1F3.sig.table.all$Facility<-ifelse(Aldex_family_16sY1_F1F3.sig.table.all$effect<0, "F3", "F1") #Add direction of significance


Aldex_family_16sY1_F2F3.sig.row<-rownames(Aldex_family_16sY1_F2F3.all)[which(Aldex_family_16sY1_F2F3.all$wi.eBH <=0.05)] #select the rownames of significant families
Aldex_family_16sY1_F2F3.sig.table<-subset(Aldex_family_16sY1_F2F3.all, rownames(Aldex_family_16sY1_F2F3.all) %in% Aldex_family_16sY1_F2F3.sig.row) #Subset significant families
# Aldex_family_16sY1_F2F3.sig.otu<-as.data.frame(subset((t(otu.n0.acomp_16sY1)), rownames(t(otu.n0.acomp_16sY1)) %in% Aldex_family_16sY1_F2F3.sig.row)) #Subset table with relative abundance per sample, only for significant fammilies
Aldex_family_16sY1_F2F3.sig.taxon<-subset(taxon_16s, rownames(taxon_16s) %in% Aldex_family_16sY1_F2F3.sig.row) #Subset taxon table to get the taxonomy of significant families
Aldex_family_16sY1_F2F3.sig.table.all<-bind_cols(Aldex_family_16sY1_F2F3.sig.taxon, Aldex_family_16sY1_F2F3.sig.table) #combine tables
Aldex_family_16sY1_F2F3.sig.table.all$Comparison<-rep("F2 v F3", 67)
Aldex_family_16sY1_F2F3.sig.table.all$Facility<-ifelse(Aldex_family_16sY1_F2F3.sig.table.all$effect<0, "F3", "F2") #Add direction of significance

# 16sY2 
Aldex_family_16sY2_F1F2.sig.row<-rownames(Aldex_family_16sY2_F1F2.all)[which(Aldex_family_16sY2_F1F2.all$wi.eBH <=0.05)] #select the rownames of significant families
Aldex_family_16sY2_F1F2.sig.table<-subset(Aldex_family_16sY2_F1F2.all, rownames(Aldex_family_16sY2_F1F2.all) %in% Aldex_family_16sY2_F1F2.sig.row) #Subset significant families
# Aldex_family_16sY2_F1F2.sig.otu<-as.data.frame(subset((t(otu.n0.acomp_16sY2)), rownames(t(otu.n0.acomp_16sY2)) %in% Aldex_family_16sY2_F1F2.sig.row)) #Subset table with relative abundance per sample, only for significant fammilies
Aldex_family_16sY2_F1F2.sig.taxon<-subset(taxon_16s, rownames(taxon_16s) %in% Aldex_family_16sY2_F1F2.sig.row) #Subset taxon table to get the taxonomy of significant families
Aldex_family_16sY2_F1F2.sig.table.all<-bind_cols(Aldex_family_16sY2_F1F2.sig.taxon, Aldex_family_16sY2_F1F2.sig.table) #combine tables
Aldex_family_16sY2_F1F2.sig.table.all$Comparison<-rep("F1 v F2", 54)
Aldex_family_16sY2_F1F2.sig.table.all$Facility<-ifelse(Aldex_family_16sY2_F1F2.sig.table.all$effect<0, "F2", "F1") #Add direction of significance


Aldex_family_16sY2_F1F3.sig.row<-rownames(Aldex_family_16sY2_F1F3.all)[which(Aldex_family_16sY2_F1F3.all$wi.eBH <=0.05)] #select the rownames of significant families
Aldex_family_16sY2_F1F3.sig.table<-subset(Aldex_family_16sY2_F1F3.all, rownames(Aldex_family_16sY2_F1F3.all) %in% Aldex_family_16sY2_F1F3.sig.row) #Subset significant families
# Aldex_family_16sY2_F1F3.sig.otu<-as.data.frame(subset((t(otu.n0.acomp_16sY2)), rownames(t(otu.n0.acomp_16sY2)) %in% Aldex_family_16sY2_F1F3.sig.row)) #Subset table with relative abundance per sample, only for significant fammilies
Aldex_family_16sY2_F1F3.sig.taxon<-subset(taxon_16s, rownames(taxon_16s) %in% Aldex_family_16sY2_F1F3.sig.row) #Subset taxon table to get the taxonomy of significant families
Aldex_family_16sY2_F1F3.sig.table.all<-bind_cols(Aldex_family_16sY2_F1F3.sig.taxon, Aldex_family_16sY2_F1F3.sig.table) #combine tables
Aldex_family_16sY2_F1F3.sig.table.all$Comparison<-rep("F1 v F3",72)
Aldex_family_16sY2_F1F3.sig.table.all$Facility<-ifelse(Aldex_family_16sY2_F1F3.sig.table.all$effect<0, "F3", "F1") #Add direction of significance


Aldex_family_16sY2_F2F3.sig.row<-rownames(Aldex_family_16sY2_F2F3.all)[which(Aldex_family_16sY2_F2F3.all$wi.eBH <=0.05)] #select the rownames of significant families
Aldex_family_16sY2_F2F3.sig.table<-subset(Aldex_family_16sY2_F2F3.all, rownames(Aldex_family_16sY2_F2F3.all) %in% Aldex_family_16sY2_F2F3.sig.row) #Subset significant families
# Aldex_family_16sY2_F2F3.sig.otu<-as.data.frame(subset((t(otu.n0.acomp_16sY2)), rownames(t(otu.n0.acomp_16sY2)) %in% Aldex_family_16sY2_F2F3.sig.row)) #Subset table with relative abundance per sample, only for significant fammilies
Aldex_family_16sY2_F2F3.sig.taxon<-subset(taxon_16s, rownames(taxon_16s) %in% Aldex_family_16sY2_F2F3.sig.row) #Subset taxon table to get the taxonomy of significant families
Aldex_family_16sY2_F2F3.sig.table.all<-bind_cols(Aldex_family_16sY2_F2F3.sig.taxon, Aldex_family_16sY2_F2F3.sig.table) #combine tables
Aldex_family_16sY2_F2F3.sig.table.all$Comparison<-rep("F2 v F3", 33)
Aldex_family_16sY2_F2F3.sig.table.all$Facility<-ifelse(Aldex_family_16sY2_F2F3.sig.table.all$effect<0, "F3", "F2") #Add direction of significance



#Effect size plots
#16sY1
Effect_family_16sY1_F1F2<-ggplot(Aldex_family_16sY1_F1F2.sig.table.all, aes(x=effect, y=reorder(Family,effect), fill=Facility))+
  geom_bar(stat='identity', color='black')+
  theme(axis.title.y = element_blank(), axis.ticks=element_line(color='black'),
        axis.text=element_text(size=10, color='black')) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  xlab("Effect size") + 
  theme(panel.background = element_rect(fill=NA, color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background = element_rect(color='black', fill='white', size=1, linetype = 'solid'))+
  theme(strip.text = element_text(size=10, face='italic'),
        panel.border = element_rect(color="black", fill=NA))+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
  scale_fill_manual(values=c("#420A68","#BB3754"))+
  theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))+
  ggtitle("Bacteria Y1 - F1 v F2", subtitle = "Differential abundant families")

Effect_family_16sY1_F1F3<-ggplot(Aldex_family_16sY1_F1F3.sig.table.all, aes(x=effect, y=reorder(Family,effect), fill=Facility))+
  geom_bar(stat='identity', color='black')+
  theme(axis.title.y = element_blank(), axis.ticks=element_line(color='black'),
        axis.text=element_text(size=10, color='black')) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  xlab("Effect size") + 
  theme(panel.background = element_rect(fill=NA, color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background = element_rect(color='black', fill='white', size=1, linetype = 'solid'))+
  theme(strip.text = element_text(size=10, face='italic'),
        panel.border = element_rect(color="black", fill=NA))+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
  scale_fill_manual(values=c("#420A68","#FCA50A"))+
  theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))+
  ggtitle("Bacteria Y1 - F1 v F3", subtitle = "Differential abundant families")

Effect_family_16sY1_F2F3<-ggplot(Aldex_family_16sY1_F2F3.sig.table.all, aes(x=effect, y=reorder(Family,effect), fill=Facility))+
  geom_bar(stat='identity', color='black')+
  theme(axis.title.y = element_blank(), axis.ticks=element_line(color='black'),
        axis.text=element_text(size=10, color='black')) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  xlab("Effect size") + 
  theme(panel.background = element_rect(fill=NA, color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background = element_rect(color='black', fill='white', size=1, linetype = 'solid'))+
  theme(strip.text = element_text(size=10, face='italic'),
        panel.border = element_rect(color="black", fill=NA))+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
  scale_fill_manual(values=c("#BB3754","#FCA50A"))+
  theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))+
  ggtitle("Bacteria Y1 - F2 v F3", subtitle = "Differential abundant families")

#16sY2
Effect_family_16sY2_F1F2<-ggplot(Aldex_family_16sY2_F1F2.sig.table.all, aes(x=effect, y=reorder(Family,effect), fill=Facility))+
  geom_bar(stat='identity', color='black')+
  theme(axis.title.y = element_blank(), axis.ticks=element_line(color='black'),
        axis.text=element_text(size=10, color='black')) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  xlab("Effect size") + 
  theme(panel.background = element_rect(fill=NA, color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background = element_rect(color='black', fill='white', size=1, linetype = 'solid'))+
  theme(strip.text = element_text(size=10, face='italic'),
        panel.border = element_rect(color="black", fill=NA))+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
  scale_fill_manual(values=c("#420A68","#BB3754"))+
  theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))+
  ggtitle("Bacteria Y2 - F1 v F2", subtitle = "Differential abundant families")

Effect_family_16sY2_F1F3<-ggplot(Aldex_family_16sY2_F1F3.sig.table.all, aes(x=effect, y=reorder(Family,effect), fill=Facility))+
  geom_bar(stat='identity', color='black')+
  theme(axis.title.y = element_blank(), axis.ticks=element_line(color='black'),
        axis.text=element_text(size=10, color='black')) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  xlab("Effect size") + 
  theme(panel.background = element_rect(fill=NA, color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background = element_rect(color='black', fill='white', size=1, linetype = 'solid'))+
  theme(strip.text = element_text(size=10, face='italic'),
        panel.border = element_rect(color="black", fill=NA))+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
  scale_fill_manual(values=c("#420A68","#FCA50A"))+
  theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))+
  ggtitle("Bacteria Y2 - F1 v F3", subtitle = "Differential abundant families")

Effect_family_16sY2_F2F3<-ggplot(Aldex_family_16sY2_F2F3.sig.table.all, aes(x=effect, y=reorder(Family,effect), fill=Facility))+
  geom_bar(stat='identity', color='black')+
  theme(axis.title.y = element_blank(), axis.ticks=element_line(color='black'),
        axis.text=element_text(size=10, color='black')) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  xlab("Effect size") + 
  theme(panel.background = element_rect(fill=NA, color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background = element_rect(color='black', fill='white', size=1, linetype = 'solid'))+
  theme(strip.text = element_text(size=10, face='italic'),
        panel.border = element_rect(color="black", fill=NA))+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
  scale_fill_manual(values=c("#BB3754","#FCA50A"))+
  theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))+
  ggtitle("Bacteria Y2 - F2 v F3", subtitle = "Differential abundant families")



#Combine significant tables
Aldex_family_16sY1_facility<-bind_rows(Aldex_family_16sY1_F1F2.sig.table.all, Aldex_family_16sY1_F1F3.sig.table.all, Aldex_family_16sY1_F2F3.sig.table.all)


#Export in .csv format



Aldex_family_16sY2_F1F2.sig.row<-rownames(Aldex_family_16sY2_F1F2.all)[which(Aldex_family_16sY2_F1F2.all$wi.eBH <=0.05)] #select the rownames of significant families
Aldex_family_16sY2_F1F3.sig.row<-rownames(Aldex_family_16sY2_F1F3.all)[which(Aldex_family_16sY2_F1F3.all$wi.eBH <=0.05)] #select the rownames of significant families
Aldex_family_16sY2_F2F3.sig.row<-rownames(Aldex_family_16sY2_F2F3.all)[which(Aldex_family_16sY2_F2F3.all$wi.eBH <=0.05)] #select the rownames of significant families

# Aldex_family_ITSY1_F1F2.sig.row<-rownames(Aldex_family_ITSY1_F1F2.all)[which(Aldex_family_ITSY1_F1F2.all$wi.eBH <=0.05)] #select the rownames of significant families
# Aldex_family_ITSY1_F1F3.sig.row<-rownames(Aldex_family_ITSY1_F1F3.all)[which(Aldex_family_ITSY1_F1F3.all$wi.eBH <=0.05)] #select the rownames of significant families
# Aldex_family_ITSY1_F2F3.sig.row<-rownames(Aldex_family_ITSY1_F2F3.all)[which(Aldex_family_ITSY1_F2F3.all$wi.eBH <=0.05)] #select the rownames of significant families
# 
# Aldex_family_ITSY2_F1F2.sig.row<-rownames(Aldex_family_ITSY2_F1F2.all)[which(Aldex_family_ITSY2_F1F2.all$wi.eBH <=0.05)] #select the rownames of significant families
# Aldex_family_ITSY2_F1F3.sig.row<-rownames(Aldex_family_ITSY2_F1F3.all)[which(Aldex_family_ITSY2_F1F3.all$wi.eBH <=0.05)] #select the rownames of significant families
# Aldex_family_ITSY2_F2F3.sig.row<-rownames(Aldex_family_ITSY2_F2F3.all)[which(Aldex_family_ITSY2_F2F3.all$wi.eBH <=0.05)] #select the rownames of significant families

#Combine list of significant OTU names and remove duplicates
Aldex_family_16sY1.sig.row<-unique(c(Aldex_family_16sY1_F1F2.sig.row, Aldex_family_16sY1_F1F3.sig.row, Aldex_family_16sY1_F2F3.sig.row))
Aldex_family_16sY2.sig.row<-unique(c(Aldex_family_16sY2_F1F2.sig.row, Aldex_family_16sY2_F1F3.sig.row, Aldex_family_16sY2_F2F3.sig.row))
# Aldex_family_ITSY1.sig.row<-unique(c(Aldex_family_ITSY1_F1F2.sig.row, Aldex_family_ITSY1_F1F3.sig.row, Aldex_family_ITSY1_F2F3.sig.row))
# Aldex_family_ITSY2.sig.row<-unique(c(Aldex_family_ITSY2_F1F2.sig.row, Aldex_family_ITSY2_F1F3.sig.row, Aldex_family_ITSY2_F2F3.sig.row))

#Subset table with relative abundance to keep only differential abundant families
SigFamilies_16sY1<-as.data.frame(subset(family_16sY1_, family_16sY1_$OTU %in% Aldex_family_16sY1.sig.row)) 
SigFamilies_16sY2<-as.data.frame(subset(family_16sY2_, family_16sY2_$OTU %in% Aldex_family_16sY2.sig.row)) 
# SigFamilies_ITSY1<-as.data.frame(subset(family_ITSY1, family_ITSY1$OTU %in% Aldex_family_ITSY1.sig.row)) 
# SigFamilies_ITSY2<-as.data.frame(subset(family_ITSY2, family_ITSY2$OTU %in% Aldex_family_ITSY2.sig.row)) 

#Calculate mean, sd and se for each Family by Facility
SigFamilies_16sY1_summarystat <- describeBy(SigFamilies_16sY1$Abundance, list(SigFamilies_16sY1$Family,SigFamilies_16sY1$Facility), mat = TRUE)
SigFamilies_16sY2_summarystat <- describeBy(SigFamilies_16sY2$Abundance, list(SigFamilies_16sY2$Family,SigFamilies_16sY2$Facility), mat = TRUE)
# SigFamilies_ITSY1_summarystat <- describeBy(SigFamilies_ITSY1$Abundance, list(SigFamilies_ITSY1$Family,SigFamilies_ITSY1$Facility), mat = TRUE)
# SigFamilies_ITSY2_summarystat <- describeBy(SigFamilies_ITSY2$Abundance, list(SigFamilies_ITSY2$Family,SigFamilies_ITSY2$Facility), mat = TRUE)


#Plot Significant families barplots of mean relative abundance by facility

SigFamilies_16sY1_barplot<-ggplot(SigFamilies_16sY1_summarystat, aes(x=group2, y=mean, fill=group2)) + 
  geom_bar(stat="identity", color="black", position=position_dodge()) + 
  facet_wrap(~group1, scales='free_y')+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(.9))+
  theme(axis.title.x = element_blank(), axis.ticks=element_line(color='black'),
        axis.text=element_text(size=12, color='black')) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  ylab("Relative Abundance (%)") + 
  theme(panel.background = element_rect(fill=NA, color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background = element_rect(color='black', fill='white', size=1, linetype = 'solid'))+
  theme(strip.text = element_text(size=10, face='italic'),
        panel.border = element_rect(color="black", fill=NA))+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
  scale_fill_viridis_d(begin = 0.2, end = 0.8, option='inferno')+
  theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))+
  ggtitle("Bacteria - Y1", subtitle = "Differential abundant families")
SigFamilies_16sY1_barplot
ggsave("Barplot_SigFam_16sY1.png", plot=SigFamilies_16sY1_barplot, device="png", width=24, height=20, units="in", dpi=600)

SigFamilies_16sY2_barplot<-ggplot(SigFamilies_16sY2_summarystat, aes(x=group2, y=mean, fill=group2)) + 
  geom_bar(stat="identity", color="black", position=position_dodge()) + 
  facet_wrap(~group1, scales='free_y')+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(.9))+
  theme(axis.title.x = element_blank(), axis.ticks=element_line(color='black'),
        axis.text=element_text(size=12, color='black')) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  ylab("Relative Abundance (%)") + 
  theme(panel.background = element_rect(fill=NA, color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background = element_rect(color='black', fill='white', size=1, linetype = 'solid'))+
  theme(strip.text = element_text(size=10, face='italic'),
        panel.border = element_rect(color="black", fill=NA))+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
  scale_fill_viridis_d(begin = 0.2, end = 0.8, option='inferno')+
  theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))+
  ggtitle("Bacteria - Y2", subtitle = "Differential abundant families")

ggsave("Barplot_SigFam_16sY2.png", plot=SigFamilies_16sY2_barplot, device="png", width=24, height=20, units="in", dpi=600)

# SigFamilies_ITSY1_barplot<-ggplot(SigFamilies_ITSY1_summarystat, aes(x=group2, y=mean, fill=group2)) + 
#   geom_bar(stat="identity", color="black", position=position_dodge()) + 
#   facet_wrap(~group1, scales='free_y')+
#   geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(.9))+
#   theme(axis.title.x = element_blank(), axis.ticks=element_line(color='black'),
#         axis.text=element_text(size=12, color='black')) + 
#   guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
#   ylab("Relative Abundance (%)") + 
#   theme(panel.background = element_rect(fill=NA, color =NA),
#         plot.background = element_rect(fill="transparent", color =NA)) +
#   theme(strip.background = element_rect(color='black', fill='white', size=1, linetype = 'solid'))+
#   theme(strip.text = element_text(size=10, face='italic'),
#         panel.border = element_rect(color="black", fill=NA))+
#   theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
#   scale_fill_viridis_d(begin = 0.2, end = 0.8, option='inferno')+
#   theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))+
#   ggtitle("Fungi - Y1", subtitle = "Differential abundant families")
# 
# ggsave("Barplot_SigFam_ITSY1.png", plot=SigFamilies_ITSY1_barplot, device="png", width=24, height=20, units="in", dpi=600)
# 
# SigFamilies_ITSY2_barplot<-ggplot(SigFamilies_ITSY2_summarystat, aes(x=group2, y=mean, fill=group2)) + 
#   geom_bar(stat="identity", color="black", position=position_dodge()) + 
#   facet_wrap(~group1, scales='free_y')+
#   geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(.9))+
#   theme(axis.title.x = element_blank(), axis.ticks=element_line(color='black'),
#         axis.text=element_text(size=12, color='black')) + 
#   guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
#   ylab("Relative Abundance (%)") + 
#   theme(panel.background = element_rect(fill=NA, color =NA),
#         plot.background = element_rect(fill="transparent", color =NA)) +
#   theme(strip.background = element_rect(color='black', fill='white', size=1, linetype = 'solid'))+
#   theme(strip.text = element_text(size=10, face='italic'),
#         panel.border = element_rect(color="black", fill=NA))+
#   theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
#   scale_fill_viridis_d(begin = 0.2, end = 0.8, option='inferno')+
#   theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))+
#   ggtitle("Fungi - Y2", subtitle = "Differential abundant families")
# 
# ggsave("Barplot_SigFam_ITSY2.png", plot=SigFamilies_ITSY2_barplot, device="png", width=24, height=20, units="in", dpi=600)


#Plots for differential abundant families above 5% relative abundance
Sigfamily_16sY1_over5<-c('Arcobacteraceae','Azospirillaceae','Caulobacteraceae','Chitinophagaceae','Enterobacteriaceae',
                         'Gammaproteobacteria_unclassified','Moraxellaceae','Pseudomonadaceae','Rhizobiaceae','Rhodobacteraceae','Sphingobacteriaceae',
                         'Sphingomonadaceae','Xanthomonadaceae')

Sigfamily_16sY2_over5<-c('Arcobacteraceae','Azospirillaceae','Beijerinckiaceae','Microbacteriaceae',
                         'Micrococcaceae','Nocardiaceae','Pseudomonadaceae','Rhizobiaceae')

# Sigfamily_ITSY1_over5<-c('Aspergillaceae','Aureobasidiaceae','Bulleribasidiaceae','Cucurbitariaceae','Didymellaceae','Didymosphaeriaceae','Dipodascaceae','Helotiales_fam_Incertae_sedis',
#                          'Mycosphaerellaceae','Pleosporaceae','Pleosporales_unclassified','Trichosporonaceae','Ascomycota_unclassified')
# 
# 
# Sigfamily_ITSY2_over5<-c('Ascomycota_unclassified','Aspergillaceae','Aureobasidiaceae','Bulleribasidiaceae',
#                          'Cucurbitariaceae','Didymosphaeriaceae','Dipodascaceae','Mucoraceae','Pichiaceae','Saccharomycetales_fam_Incertae_sedis',
#                          'Sporidiobolaceae','Trichomeriaceae','Trichosporonaceae','Pleosporales_unclassified')

Sigfamily_16sY1_over5abund <- filter(SigFamilies_16sY1, Family %in% Sigfamily_16sY1_over5)
Sigfamily_16sY2_over5abund <- filter(SigFamilies_16sY2, Family %in% Sigfamily_16sY2_over5)
# Sigfamily_ITSY1_over5abund <- filter(SigFamilies_ITSY1, Family %in% Sigfamily_ITSY1_over5)
# Sigfamily_ITSY2_over5abund <- filter(SigFamilies_ITSY2, Family %in% Sigfamily_ITSY2_over5)



#Calculate mean, sd and se for each Family by Facility
SigFamilies_16sY1_over5_summarystat <- describeBy(Sigfamily_16sY1_over5abund$Abundance, list(Sigfamily_16sY1_over5abund$Family,Sigfamily_16sY1_over5abund$Facility), mat = TRUE)
SigFamilies_16sY2_over5_summarystat <- describeBy(Sigfamily_16sY2_over5abund$Abundance, list(Sigfamily_16sY2_over5abund$Family,Sigfamily_16sY2_over5abund$Facility), mat = TRUE)
# SigFamilies_ITSY1_over5_summarystat <- describeBy(Sigfamily_ITSY1_over5abund$Abundance, list(Sigfamily_ITSY1_over5abund$Family,Sigfamily_ITSY1_over5abund$Facility), mat = TRUE)
# SigFamilies_ITSY2_over5_summarystat <- describeBy(Sigfamily_ITSY2_over5abund$Abundance, list(Sigfamily_ITSY2_over5abund$Family,Sigfamily_ITSY2_over5abund$Facility), mat = TRUE)



#Plot Significant families barplots of mean relative abundance by facility

SigFamilies_16sY1_over5_barplot<-ggplot(SigFamilies_16sY1_over5_summarystat, aes(x=group2, y=mean, fill=group2)) + 
  geom_bar(stat="identity", color="black", position=position_dodge()) + 
  facet_wrap(~group1, scales='free_y')+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(.9))+
  theme(axis.title.x = element_blank(), axis.ticks=element_line(color='black'),
        axis.text=element_text(size=12, color='black')) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  ylab("Relative Abundance (%)") + 
  theme(panel.background = element_rect(fill=NA, color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background = element_rect(color='black', fill='white', size=1, linetype = 'solid'))+
  theme(strip.text = element_text(size=10, face='italic'),
        panel.border = element_rect(color="black", fill=NA))+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
  scale_fill_viridis_d(begin = 0.2, end = 0.8, option='inferno')+
  theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))+
  guides(fill=guide_legend(title="Facility"))+
  ggtitle("Bacteria - Y1", subtitle = "Differential abundant families above 5% relative abundace")
ggsave("Barplot_SigFam_over5_16sY1.png", plot=SigFamilies_16sY1_over5_barplot, device="png", width=10, height=8, units="in", dpi=600)

SigFamilies_16sY2_over5_barplot<-ggplot(SigFamilies_16sY2_over5_summarystat, aes(x=group2, y=mean, fill=group2)) + 
  geom_bar(stat="identity", color="black", position=position_dodge()) + 
  facet_wrap(~group1, scales='free_y')+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(.9))+
  theme(axis.title.x = element_blank(), axis.ticks=element_line(color='black'),
        axis.text=element_text(size=12, color='black')) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  ylab("Relative Abundance (%)") + 
  theme(panel.background = element_rect(fill=NA, color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background = element_rect(color='black', fill='white', size=1, linetype = 'solid'))+
  theme(strip.text = element_text(size=10, face='italic'),
        panel.border = element_rect(color="black", fill=NA))+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
  scale_fill_viridis_d(begin = 0.2, end = 0.8, option='inferno')+
  theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))+
  guides(fill=guide_legend(title="Facility"))+
  ggtitle("Bacteria - Y2", subtitle = "Differential abundant families above 5% relative abundace")
ggsave("Barplot_SigFam_over5_16sY2.png", plot=SigFamilies_16sY2_over5_barplot, device="png", width=10, height=8, units="in", dpi=600)

# SigFamilies_ITSY1_over5_barplot<-ggplot(SigFamilies_ITSY1_over5_summarystat, aes(x=group2, y=mean, fill=group2)) + 
#   geom_bar(stat="identity", color="black", position=position_dodge()) + 
#   facet_wrap(~group1, scales='free_y')+
#   geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(.9))+
#   theme(axis.title.x = element_blank(), axis.ticks=element_line(color='black'),
#         axis.text=element_text(size=12, color='black')) + 
#   guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
#   ylab("Relative Abundance (%)") + 
#   theme(panel.background = element_rect(fill=NA, color =NA),
#         plot.background = element_rect(fill="transparent", color =NA)) +
#   theme(strip.background = element_rect(color='black', fill='white', size=1, linetype = 'solid'))+
#   theme(strip.text = element_text(size=10, face='italic'),
#         panel.border = element_rect(color="black", fill=NA))+
#   theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
#   scale_fill_viridis_d(begin = 0.2, end = 0.8, option='inferno')+
#   theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))+
#   guides(fill=guide_legend(title="Facility"))+
#   ggtitle("Fungi - Y1", subtitle = "Differential abundant families above 5% relative abundace")
# ggsave("Barplot_SigFam_over5_ITSY1.png", plot=SigFamilies_ITSY1_over5_barplot, device="png", width=10, height=8, units="in", dpi=600)
# 
# SigFamilies_ITSY2_over5_barplot<-ggplot(SigFamilies_ITSY2_over5_summarystat, aes(x=group2, y=mean, fill=group2)) + 
#   geom_bar(stat="identity", color="black", position=position_dodge()) + 
#   facet_wrap(~group1, scales='free_y')+
#   geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(.9))+
#   theme(axis.title.x = element_blank(), axis.ticks=element_line(color='black'),
#         axis.text=element_text(size=12, color='black')) + 
#   guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
#   ylab("Relative Abundance (%)") + 
#   theme(panel.background = element_rect(fill=NA, color =NA),
#         plot.background = element_rect(fill="transparent", color =NA)) +
#   theme(strip.background = element_rect(color='black', fill='white', size=1, linetype = 'solid'))+
#   theme(strip.text = element_text(size=10, face='italic'),
#         panel.border = element_rect(color="black", fill=NA))+
#   theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
#   scale_fill_viridis_d(begin = 0.2, end = 0.8, option='inferno')+
#   theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))+
#   guides(fill=guide_legend(title="Facility"))+
#   ggtitle("Fungi - Y2", subtitle = "Differential abundant families above 5% relative abundace")
# ggsave("Barplot_SigFam_over5_ITSY2.png", plot=SigFamilies_ITSY2_over5_barplot, device="png", width=10, height=8, units="in", dpi=600)
# 


#### Compositional analysis by facility throughout seasons at the OTU level ####

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

#Remove otus with 0 counts in all samples
otus_16sF1<-otus_16sF1[ which(rowSums(otus_16sF1)>0),]
otus_16sF2<-otus_16sF2[ which(rowSums(otus_16sF2)>0),]
otus_16sF3<-otus_16sF3[ which(rowSums(otus_16sF3)>0),]

otus_ITSF1<-otus_ITSF1[ which(rowSums(otus_ITSF1)>0),]
otus_ITSF2<-otus_ITSF2[ which(rowSums(otus_ITSF2)>0),]
otus_ITSF3<-otus_ITSF3[ which(rowSums(otus_ITSF3)>0),]

#Subset metadata
metadata_16sF1<-subset(metadata_16s, Facility=="F1")
metadata_16sF2<-subset(metadata_16s, Facility=="F2")
metadata_16sF3<-subset(metadata_16s, Facility=="F3")

metadata_ITSF1<-subset(metadata_ITS, Facility=="F1")
metadata_ITSF2<-subset(metadata_ITS, Facility=="F2")
metadata_ITSF3<-subset(metadata_ITS, Facility=="F3")


#Step 1: Convert OTU table to appropriate format
#Following step requires samples on rows and OTUs in columns
head(t(otus_16sF1)) 
head(t(otus_16sF2)) 
head(t(otus_16sF3)) 

head(t(otus_ITSF1)) 
head(t(otus_ITSF2)) 
head(t(otus_ITSF3)) 

#Step 2: Replace zero values before clr transformation. 
#Uses Count Zero Multiplicative method to replace zeros and outputs pseudo-counts
otu.n0_16sF1<-t(cmultRepl(t(otus_16sF1), label=0, method="CZM", output="p-counts")) #214085   corrected values
otu.n0_16sF2<-t(cmultRepl(t(otus_16sF2), label=0, method="CZM", output="p-counts")) #136091 
otu.n0_16sF3<-t(cmultRepl(t(otus_16sF3), label=0, method="CZM", output="p-counts")) #360590

otu.n0_ITSF1<-t(cmultRepl(t(otus_ITSF1), label=0, method="CZM", output="p-counts")) #140754  corrected values
otu.n0_ITSF2<-t(cmultRepl(t(otus_ITSF2), label=0, method="CZM", output="p-counts")) #106025
otu.n0_ITSF3<-t(cmultRepl(t(otus_ITSF3), label=0, method="CZM", output="p-counts")) #140912

#Check output table - needs to have samples in columns and OTUs in rows
head(otu.n0_16sF1) 
head(otu.n0_16sF2) 
head(otu.n0_16sF3) 

head(otu.n0_ITSF1) 
head(otu.n0_ITSF2) 
head(otu.n0_ITSF3) 

#Step 3: Convert data to proportions
otu.n0_16sF1_prop<-apply(otu.n0_16sF1, 2, function(x) {x/sum(x)})
otu.n0_16sF2_prop<-apply(otu.n0_16sF2, 2, function(x) {x/sum(x)})
otu.n0_16sF3_prop<-apply(otu.n0_16sF3, 2, function(x) {x/sum(x)})

otu.n0_ITSF1_prop<-apply(otu.n0_ITSF1, 2, function(x) {x/sum(x)})
otu.n0_ITSF2_prop<-apply(otu.n0_ITSF2, 2, function(x) {x/sum(x)})
otu.n0_ITSF3_prop<-apply(otu.n0_ITSF3, 2, function(x) {x/sum(x)})

head(otu.n0_16sF1_prop)
head(otu.n0_16sF2_prop)
head(otu.n0_16sF3_prop)

head(otu.n0_ITSF1_prop)
head(otu.n0_ITSF2_prop)
head(otu.n0_ITSF3_prop)

#Step 4: Perform abundance and sample filtering and deal sparsity
#Filter the data to remove all taxa with less than 0.00001% abundance in any sample
otu.n0_16sF1_prop_f<-otu.n0_16sF1[apply(otu.n0_16sF1_prop, 1, min) > 0.00000001, ]
otu.n0_16sF2_prop_f<-otu.n0_16sF2[apply(otu.n0_16sF2_prop, 1, min) > 0.00000001, ]
otu.n0_16sF3_prop_f<-otu.n0_16sF3[apply(otu.n0_16sF3_prop, 1, min) > 0.00000001, ]

otu.n0_ITSF1_prop_f<-otu.n0_ITSF1[apply(otu.n0_ITSF1_prop, 1, min) > 0.00000001, ]
otu.n0_ITSF2_prop_f<-otu.n0_ITSF2[apply(otu.n0_ITSF2_prop, 1, min) > 0.00000001, ]
otu.n0_ITSF3_prop_f<-otu.n0_ITSF3[apply(otu.n0_ITSF3_prop, 1, min) > 0.00000001, ]

#Check that samples are on columns and OTUs in rows
head(otu.n0_16sF1_prop_f) 
head(otu.n0_16sF2_prop_f) 
head(otu.n0_16sF3_prop_f) 

head(otu.n0_ITSF1_prop_f) 
head(otu.n0_ITSF2_prop_f) 
head(otu.n0_ITSF3_prop_f)  

#Step 5: perform CLR transformation
otu.n0.clr_16sF1<-t(apply(otu.n0_16sF1_prop_f, 2, function(x){log(x)-mean(log(x))}))
otu.n0.clr_16sF2<-t(apply(otu.n0_16sF2_prop_f, 2, function(x){log(x)-mean(log(x))}))
otu.n0.clr_16sF3<-t(apply(otu.n0_16sF3_prop_f, 2, function(x){log(x)-mean(log(x))}))

otu.n0.clr_ITSF1<-t(apply(otu.n0_ITSF1_prop_f, 2, function(x){log(x)-mean(log(x))}))
otu.n0.clr_ITSF2<-t(apply(otu.n0_ITSF2_prop_f, 2, function(x){log(x)-mean(log(x))}))
otu.n0.clr_ITSF3<-t(apply(otu.n0_ITSF3_prop_f, 2, function(x){log(x)-mean(log(x))}))

#Check output table. Samples should be in rows and OTUs in columns
head(otu.n0.clr_16sF1) 
head(otu.n0.clr_16sF2) 
head(otu.n0.clr_16sF3) 

head(otu.n0.clr_ITSF1) 
head(otu.n0.clr_ITSF2) 
head(otu.n0.clr_ITSF3)

#Step 6: Perform Singular Value Decomposition (PCA)
pc.clr_16sF1<-prcomp(otu.n0.clr_16sF1)
pc.clr_16sF2<-prcomp(otu.n0.clr_16sF2)
pc.clr_16sF3<-prcomp(otu.n0.clr_16sF3)

pc.clr_ITSF1<-prcomp(otu.n0.clr_ITSF1)
pc.clr_ITSF2<-prcomp(otu.n0.clr_ITSF2)
pc.clr_ITSF3<-prcomp(otu.n0.clr_ITSF3)

png("Screeplot - PCA by facility .png", width = 1000, height = 500, units = 'px')
par(mar=c(2,2,2,2))
par(mfrow=c(2,3))
screeplot(pc.clr_16sF1, type='lines', main="Bacteria F1")
screeplot(pc.clr_16sF2, type='lines', main="Bacteria F2")
screeplot(pc.clr_16sF3, type='lines', main="Bacteria F3")
screeplot(pc.clr_ITSF1, type='lines', main="Fungi F1")
screeplot(pc.clr_ITSF2, type='lines', main="Fungi F2")
screeplot(pc.clr_ITSF3, type='lines', main="Fungi F3")
dev.off()


#Calculate total variance of the data
mvar.clr_16sF1<-mvar(otu.n0.clr_16sF1)
mvar.clr_16sF2<-mvar(otu.n0.clr_16sF2)
mvar.clr_16sF3<-mvar(otu.n0.clr_16sF3)

mvar.clr_ITSF1<-mvar(otu.n0.clr_ITSF1)
mvar.clr_ITSF2<-mvar(otu.n0.clr_ITSF2)
mvar.clr_ITSF3<-mvar(otu.n0.clr_ITSF3)


#Display results - 16sF1
row_16sF1<-rownames(otu.n0.clr_16sF1) #Make vector with sample names
pc_out_16sF1<-as.data.frame(pc.clr_16sF1$x[,1:2]) #Get PC1 and PC2 
pc_out_meta_16sF1<-as.data.frame(bind_cols(pc_out_16sF1,metadata_16sF1)) #Add metadata information
row.names(pc_out_meta_16sF1)<-row_16sF1 #Add rownames to dataframe
pc_out_meta_16sF1$Facility<-as.factor(pc_out_meta_16sF1$Facility)
pc_out_meta_16sF1$Year<-as.factor(pc_out_meta_16sF1$Year)
pc_out_meta_16sF1$Week<-as.factor(pc_out_meta_16sF1$Week)
pc_out_meta_16sF1$Month<-as.factor(pc_out_meta_16sF1$Month)
pc_out_meta_16sF1$Collection<-as.factor(pc_out_meta_16sF1$Collection)

# Make PCA plot - First 2 axis- color by facility and shape by L..mono
PCA_16sF1_date <- ggplot(pc_out_meta_16sF1, aes(x=PC1,y=PC2, color=Month, shape=Year))+
  geom_point(size=3)+
  theme(legend.text=element_text(size=13,color='black'), legend.title= element_text(size=13, face='italic')) +
  theme(legend.position = 'right')+
  theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  theme(axis.title = element_text(size=13,color='black')) +
  theme(panel.background = element_rect(fill='transparent', color = NA),
        plot.background = element_rect(fill = 'transparent',color = NA), 
        panel.border = element_rect(color='black',fill = NA,size=1))+
  scale_x_continuous(name = paste("PC1: ", round(pc.clr_16sF1$sdev[1]^2/mvar.clr_16sF1*100, digits=1), "%", sep="")) +
  scale_y_continuous(name = paste("PC2: ", round(pc.clr_16sF1$sdev[2]^2/mvar.clr_16sF1*100, digits=1), "%", sep="")) +
  ggtitle("Bacteria - F1", subtitle = "PCA by time")+theme(plot.title = element_text(hjust = 0.5, face = 'bold'))+
  scale_color_viridis_d(begin = 0.1, end = 0.9, option='viridis')
PCA_16sF1_date

#Display results - 16sF2
row_16sF2<-rownames(otu.n0.clr_16sF2) #Make vector with sample names
pc_out_16sF2<-as.data.frame(pc.clr_16sF2$x[,1:2]) #Get PC1 and PC2 
pc_out_meta_16sF2<-as.data.frame(bind_cols(pc_out_16sF2,metadata_16sF2)) #Add metadata information
row.names(pc_out_meta_16sF2)<-row_16sF2 #Add rownames to dataframe
pc_out_meta_16sF2$Facility<-as.factor(pc_out_meta_16sF2$Facility)
pc_out_meta_16sF2$Year<-as.factor(pc_out_meta_16sF2$Year)
pc_out_meta_16sF2$Week<-as.factor(pc_out_meta_16sF2$Week)
pc_out_meta_16sF2$Month<-as.factor(pc_out_meta_16sF2$Month)
pc_out_meta_16sF2$Collection<-as.factor(pc_out_meta_16sF2$Collection)

# Make PCA plot - First 2 axis- color by Year and month
PCA_16sF2_date <- ggplot(pc_out_meta_16sF2, aes(x=PC1,y=PC2, color=Month, shape=Year))+
  geom_point(size=3)+
  theme(legend.text=element_text(size=13,color='black'), legend.title= element_text(size=13, face='italic')) +
  theme(legend.position = 'right')+
  theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  theme(axis.title = element_text(size=13,color='black')) +
  theme(panel.background = element_rect(fill='transparent', color = NA),
        plot.background = element_rect(fill = 'transparent',color = NA), 
        panel.border = element_rect(color='black',fill = NA,size=1))+
  scale_x_continuous(name = paste("PC1: ", round(pc.clr_16sF2$sdev[1]^2/mvar.clr_16sF2*100, digits=1), "%", sep="")) +
  scale_y_continuous(name = paste("PC2: ", round(pc.clr_16sF2$sdev[2]^2/mvar.clr_16sF2*100, digits=1), "%", sep="")) +
  ggtitle("Bacteria - F2", subtitle = "PCA by time")+theme(plot.title = element_text(hjust = 0.5, face = 'bold'))+
  scale_color_viridis_d(begin = 0.1, end = 0.9, option='viridis')
PCA_16sF2_date

#Display results - 16sF3
row_16sF3<-rownames(otu.n0.clr_16sF3) #Make vector with sample names
pc_out_16sF3<-as.data.frame(pc.clr_16sF3$x[,1:2]) #Get PC1 and PC2 
pc_out_meta_16sF3<-as.data.frame(bind_cols(pc_out_16sF3,metadata_16sF3)) #Add metadata information
row.names(pc_out_meta_16sF3)<-row_16sF3 #Add rownames to dataframe
pc_out_meta_16sF3$Facility<-as.factor(pc_out_meta_16sF3$Facility)
pc_out_meta_16sF3$Year<-as.factor(pc_out_meta_16sF3$Year)
pc_out_meta_16sF3$Week<-as.factor(pc_out_meta_16sF3$Week)
pc_out_meta_16sF3$Month<-as.factor(pc_out_meta_16sF3$Month)
pc_out_meta_16sF3$Collection<-as.factor(pc_out_meta_16sF3$Collection)

# Make PCA plot - First 2 axis- color by Year and month
PCA_16sF3_date <- ggplot(pc_out_meta_16sF3, aes(x=PC1,y=PC2, color=Month, shape=Year))+
  geom_point(size=3)+
  theme(legend.text=element_text(size=13,color='black'), legend.title= element_text(size=13, face='italic')) +
  theme(legend.position = 'right')+
  theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  theme(axis.title = element_text(size=13,color='black')) +
  theme(panel.background = element_rect(fill='transparent', color = NA),
        plot.background = element_rect(fill = 'transparent',color = NA), 
        panel.border = element_rect(color='black',fill = NA,size=1))+
  scale_x_continuous(name = paste("PC1: ", round(pc.clr_16sF3$sdev[1]^2/mvar.clr_16sF3*100, digits=1), "%", sep="")) +
  scale_y_continuous(name = paste("PC2: ", round(pc.clr_16sF3$sdev[2]^2/mvar.clr_16sF3*100, digits=1), "%", sep="")) +
  ggtitle("Bacteria - F3", subtitle = "PCA by time")+theme(plot.title = element_text(hjust = 0.5, face = 'bold'))+
  scale_color_viridis_d(begin = 0.1, end = 0.9, option='viridis')
PCA_16sF3_date

#Display results - ITSF1
row_ITSF1<-rownames(otu.n0.clr_ITSF1) #Make vector with sample names
pc_out_ITSF1<-as.data.frame(pc.clr_ITSF1$x[,1:2]) #Get PC1 and PC2 
pc_out_meta_ITSF1<-as.data.frame(bind_cols(pc_out_ITSF1,metadata_ITSF1)) #Add metadata information
row.names(pc_out_meta_ITSF1)<-row_ITSF1 #Add rownames to dataframe
pc_out_meta_ITSF1$Facility<-as.factor(pc_out_meta_ITSF1$Facility)
pc_out_meta_ITSF1$Year<-as.factor(pc_out_meta_ITSF1$Year)
pc_out_meta_ITSF1$Week<-as.factor(pc_out_meta_ITSF1$Week)
pc_out_meta_ITSF1$Month<-as.factor(pc_out_meta_ITSF1$Month)
pc_out_meta_ITSF1$Collection<-as.factor(pc_out_meta_ITSF1$Collection)

# Make PCA plot - First 2 axis- color by facility and shape by L..mono
PCA_ITSF1_date <- ggplot(pc_out_meta_ITSF1, aes(x=PC1,y=PC2, color=Month, shape=Year))+
  geom_point(size=3)+
  theme(legend.text=element_text(size=13,color='black'), legend.title= element_text(size=13, face='italic')) +
  theme(legend.position = 'right')+
  theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  theme(axis.title = element_text(size=13,color='black')) +
  theme(panel.background = element_rect(fill='transparent', color = NA),
        plot.background = element_rect(fill = 'transparent',color = NA), 
        panel.border = element_rect(color='black',fill = NA,size=1))+
  scale_x_continuous(name = paste("PC1: ", round(pc.clr_ITSF1$sdev[1]^2/mvar.clr_ITSF1*100, digits=1), "%", sep="")) +
  scale_y_continuous(name = paste("PC2: ", round(pc.clr_ITSF1$sdev[2]^2/mvar.clr_ITSF1*100, digits=1), "%", sep="")) +
  ggtitle("Fungi - F1", subtitle = "PCA by time")+theme(plot.title = element_text(hjust = 0.5, face = 'bold'))+
  scale_color_viridis_d(begin = 0.1, end = 0.9, option='viridis')
PCA_ITSF1_date

#Display results - ITSF2
row_ITSF2<-rownames(otu.n0.clr_ITSF2) #Make vector with sample names
pc_out_ITSF2<-as.data.frame(pc.clr_ITSF2$x[,1:2]) #Get PC1 and PC2 
pc_out_meta_ITSF2<-as.data.frame(bind_cols(pc_out_ITSF2,metadata_ITSF2)) #Add metadata information
row.names(pc_out_meta_ITSF2)<-row_ITSF2 #Add rownames to dataframe
pc_out_meta_ITSF2$Facility<-as.factor(pc_out_meta_ITSF2$Facility)
pc_out_meta_ITSF2$Year<-as.factor(pc_out_meta_ITSF2$Year)
pc_out_meta_ITSF2$Week<-as.factor(pc_out_meta_ITSF2$Week)
pc_out_meta_ITSF2$Month<-as.factor(pc_out_meta_ITSF2$Month)
pc_out_meta_ITSF2$Collection<-as.factor(pc_out_meta_ITSF2$Collection)

# Make PCA plot - First 2 axis- color by Year and month
PCA_ITSF2_date <- ggplot(pc_out_meta_ITSF2, aes(x=PC1,y=PC2, color=Month, shape=Year))+
  geom_point(size=3)+
  theme(legend.text=element_text(size=13,color='black'), legend.title= element_text(size=13, face='italic')) +
  theme(legend.position = 'right')+
  theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  theme(axis.title = element_text(size=13,color='black')) +
  theme(panel.background = element_rect(fill='transparent', color = NA),
        plot.background = element_rect(fill = 'transparent',color = NA), 
        panel.border = element_rect(color='black',fill = NA,size=1))+
  scale_x_continuous(name = paste("PC1: ", round(pc.clr_ITSF2$sdev[1]^2/mvar.clr_ITSF2*100, digits=1), "%", sep="")) +
  scale_y_continuous(name = paste("PC2: ", round(pc.clr_ITSF2$sdev[2]^2/mvar.clr_ITSF2*100, digits=1), "%", sep="")) +
  ggtitle("Fungi - F2", subtitle = "PCA by time")+theme(plot.title = element_text(hjust = 0.5, face = 'bold'))+
  scale_color_viridis_d(begin = 0.1, end = 0.9, option='viridis')
PCA_ITSF2_date

#Display results - ITSF3
row_ITSF3<-rownames(otu.n0.clr_ITSF3) #Make vector with sample names
pc_out_ITSF3<-as.data.frame(pc.clr_ITSF3$x[,1:2]) #Get PC1 and PC2 
pc_out_meta_ITSF3<-as.data.frame(bind_cols(pc_out_ITSF3,metadata_ITSF3)) #Add metadata information
row.names(pc_out_meta_ITSF3)<-row_ITSF3 #Add rownames to dataframe
pc_out_meta_ITSF3$Facility<-as.factor(pc_out_meta_ITSF3$Facility)
pc_out_meta_ITSF3$Year<-as.factor(pc_out_meta_ITSF3$Year)
pc_out_meta_ITSF3$Week<-as.factor(pc_out_meta_ITSF3$Week)
pc_out_meta_ITSF3$Month<-as.factor(pc_out_meta_ITSF3$Month)
pc_out_meta_ITSF3$Collection<-as.factor(pc_out_meta_ITSF3$Collection)

# Make PCA plot - First 2 axis- color by Year and month
PCA_ITSF3_date <- ggplot(pc_out_meta_ITSF3, aes(x=PC1,y=PC2, color=Month, shape=Year))+
  geom_point(size=3)+
  theme(legend.text=element_text(size=13,color='black'), legend.title= element_text(size=13, face='italic')) +
  theme(legend.position = 'right')+
  theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  theme(axis.title = element_text(size=13,color='black')) +
  theme(panel.background = element_rect(fill='transparent', color = NA),
        plot.background = element_rect(fill = 'transparent',color = NA), 
        panel.border = element_rect(color='black',fill = NA,size=1))+
  scale_x_continuous(name = paste("PC1: ", round(pc.clr_ITSF3$sdev[1]^2/mvar.clr_ITSF3*100, digits=1), "%", sep="")) +
  scale_y_continuous(name = paste("PC2: ", round(pc.clr_ITSF3$sdev[2]^2/mvar.clr_ITSF3*100, digits=1), "%", sep="")) +
  ggtitle("Fungi - F3", subtitle = "PCA by time")+theme(plot.title = element_text(hjust = 0.5, face = 'bold'))+
  scale_color_viridis_d(begin = 0.1, end = 0.9, option='viridis')
PCA_ITSF3_date



PCA_all_date = plot_grid(PCA_16sF1_date, PCA_16sF2_date, PCA_16sF3_date,
                         PCA_ITSF1_date, PCA_ITSF2_date, PCA_ITSF3_date,
                         ncol=3, nrow=2)
PCA_all_date


ggsave("PC_all_date.png", plot =PCA_all_date, device="png", width=16, height=8, units="in",dpi=600)


#Beta diversity over time
#Generate the distance matrix
dist_16sF1<-dist(otu.n0.clr_16sF1, method='euclidean')
dist_16sF2<-dist(otu.n0.clr_16sF2, method='euclidean')
dist_16sF3<-dist(otu.n0.clr_16sF3, method='euclidean')

dist_ITSF1<-dist(otu.n0.clr_ITSF1, method='euclidean')
dist_ITSF2<-dist(otu.n0.clr_ITSF2, method='euclidean')
dist_ITSF3<-dist(otu.n0.clr_ITSF3, method='euclidean')
# 
# ## Hierarchical Clustering #
# #Generate the distance matrix
# dist_16sY1_<-dist(otu.n0.clr_16sY1_, method='euclidean') #Distance matrix based on Aitchinson simplex
# # dist_16sY1_v138<-dist(otu.n0.clr_16sY1_v138, method='euclidean')
# dist_ITSY1<-dist(otu.n0.clr_ITSY1, method='euclidean')
# 
# dist_16sY2_<-dist(otu.n0.clr_16sY2_, method='euclidean') #Distance matrix based on Aitchinson simplex
# # dist_16sY2_v138<-dist(otu.n0.clr_16sY2_v138, method='euclidean')
# dist_ITSY2<-dist(otu.n0.clr_ITSY2, method='euclidean')
# 
# #Cluster the data
# hc_16sY1_<-hclust(dist_16sY1_, method = 'ward.D2')
# # hc_16sY1_v138<-hclust(dist_16sY1_v138, method = 'ward.D2')
# hc_ITSY1<-hclust(dist_ITSY1, method = 'ward.D2')
# 
# hc_16sY2_<-hclust(dist_16sY2_, method = 'ward.D2')
# # hc_16sY2_v138<-hclust(dist_16sY2_v138, method = 'ward.D2')
# hc_ITSY2<-hclust(dist_ITSY2, method = 'ward.D2')
# 
# 
# #Plot the dendrogram
# metadata_16sY1$Facility<-as.factor(metadata_16sY1$Facility) #Make Facility a factor
# metadata_16sY1$L..monocytogenes<-as.factor(metadata_16sY1$L..monocytogenes) #make L.mono a factor
# metadata_16sY1$Week<-as.factor(metadata_16sY1$Week) #make Year a factor
# 
# metadata_16sY2$Facility<-as.factor(metadata_16sY2$Facility) #Make Facility a factor
# metadata_16sY2$L..monocytogenes<-as.factor(metadata_16sY2$L..monocytogenes) #make L.mono a factor
# metadata_16sY2$Week<-as.factor(metadata_16sY2$Week) #make Year a factor
# 
# 
# metadata_ITSY1$Facility<-as.factor(metadata_ITSY1$Facility) #Make Facility a factor
# metadata_ITSY1$L..monocytogenes<-as.factor(metadata_ITSY1$L..monocytogenes) #make L.mono a factor
# metadata_ITSY1$Week<-as.factor(metadata_ITSY1$Week) #make Year a factor
# 
# metadata_ITSY2$Facility<-as.factor(metadata_ITSY2$Facility) #Make Facility a factor
# metadata_ITSY2$L..monocytogenes<-as.factor(metadata_ITSY2$L..monocytogenes) #make L.mono a factor
# metadata_ITSY2$Week<-as.factor(metadata_ITSY2$Week) #make Year a factor
# 
# 
# #Make vector with colors y facility and order them as in dendrogram
# cols_fac <- c('#420A68','#BB3754','#FCA50A') #colors for facilities
# col_16sY1__fac <- (cols_fac[metadata_16sY1$Facility])[order.dendrogram(as.dendrogram(hc_16sY1_))]
# # col_16sY1_v138_fac <- (cols_fac[metadata_16sY1$Facility])[order.dendrogram(as.dendrogram(hc_16sY1_v138))]
# col_ITSY1_fac <- (cols_fac[metadata_ITSY1$Facility])[order.dendrogram(as.dendrogram(hc_ITSY1))]
# 
# col_16sY2__fac <- (cols_fac[metadata_16sY2$Facility])[order.dendrogram(as.dendrogram(hc_16sY2_))]
# # col_16sY2_v138_fac <- (cols_fac[metadata_16sY2$Facility])[order.dendrogram(as.dendrogram(hc_16sY2_v138))]
# col_ITSY2_fac <- (cols_fac[metadata_ITSY2$Facility])[order.dendrogram(as.dendrogram(hc_ITSY2))]
# 
# #Make vector with colors L. monocytogenes and order them as in dendrogram
# cols_Lm <- c('#F8CD9C','#EA7580') #colors for year
# col_16sY1__Lm <- (cols_Lm[metadata_16sY1$L..monocytogenes])[order.dendrogram(as.dendrogram(hc_16sY1_))]
# # col_16sY1_v138_Lm <- (cols_Lm[metadata_16sY1$L..monocytogenes])[order.dendrogram(as.dendrogram(hc_16sY1_v138))]
# col_ITSY1_Lm <- (cols_Lm[metadata_ITSY1$L..monocytogenes])[order.dendrogram(as.dendrogram(hc_ITSY1))]
# 
# col_16sY2__Lm <- (cols_Lm[metadata_16sY2$L..monocytogenes])[order.dendrogram(as.dendrogram(hc_16sY2_))]
# # col_16sY2_v138_Lm <- (cols_Lm[metadata_16sY2$L..monocytogenes])[order.dendrogram(as.dendrogram(hc_16sY2_v138))]
# col_ITSY2_Lm <- (cols_Lm[metadata_ITSY2$L..monocytogenes])[order.dendrogram(as.dendrogram(hc_ITSY2))]
# 
# # #Make vector with colors Week and order them as in dendrogram
# # cols_week <- c("#482576","#443A83","#3D4D8A","#35608D","#2D718E","#27808E","#21908C","#1FA188","#29AF7F","#43BF71","#66CB5D","#8FD744","#BBDF27") #colors for facilities
# # col_16sY1__week <- (cols_week[metadata_16sY1$Week])[order.dendrogram(as.dendrogram(hc_16sY1_))]
# # # col_16sY1_v138_week <- (cols_week[metadata_16sY1$Week])[order.dendrogram(as.dendrogram(hc_16sY1_v138))]
# # col_ITS_week <- (cols_week[metadata_ITS$Week])[order.dendrogram(as.dendrogram(hc_ITS))]
# 
# # col_16sY2__week <- (cols_week[metadata_16sY2$Week])[order.dendrogram(as.dendrogram(hc_16sY2_))]
# # col_16sY2_v138_week <- (cols_week[metadata_16sY2$Week])[order.dendrogram(as.dendrogram(hc_16sY2_v138))]
# # # col_ITS_week <- (cols_week[metadata_ITS$Week])[order.dendrogram(as.dendrogram(hc_ITS))]
# # 
# 
# #Make dendrogram object
# dendro_16sY1_<- as.dendrogram(hc_16sY1_) %>%
#   set("labels", "")
# dendro_16sY2_<- as.dendrogram(hc_16sY2_) %>%
#   set("labels", "")
# dendro_ITSY1<- as.dendrogram(hc_ITSY1) %>%
#   set("labels", "")
# dendro_ITSY2<- as.dendrogram(hc_ITSY2) %>%
#   set("labels", "")
# 
# 
# #plot and save
# png("Dendrogram - 16sY1_.png", width = 8, height = 4, units = 'in', res=600)
# par(mar=c(3,6,3,4), xpd=T)
# plot(dendro_16sY1_, horiz=FALSE, main= "Bacteria Y1 - ", axes=FALSE)
# colored_bars(colors = cbind(col_16sY1__Lm,col_16sY1__fac), dend=dendro_16sY1_, rowLabels = c('L. monocytogenes','Facility'), 
#              sort_by_labels_order = FALSE)
# legend("topright", inset = c(-0.1,-0.1), legend = c('-','+','F1','F2','F3'), 
#        fill = c('#F8CD9C','#EA7580','#420A68','#BB3754','#FCA50A'), bty='n')
# dev.off() 
# 
# png("Dendrogram - 16sY2_.png", width = 8, height = 4, units = 'in', res=600)
# par(mar=c(3,6,3,4), xpd=T)
# plot(dendro_16sY2_, horiz=FALSE, main= "Bacteria Y2 - ", axes=FALSE)
# colored_bars(colors = cbind(col_16sY2__Lm,col_16sY2__fac), dend=dendro_16sY2_, rowLabels = c('L. monocytogenes','Facility'), 
#              sort_by_labels_order = FALSE)
# legend("topright", inset = c(-0.1,-0.1), legend = c('-','+','F1','F2','F3'), 
#        fill = c('#F8CD9C','#EA7580','#420A68','#BB3754','#FCA50A'), bty='n')
# dev.off()  
# 
# png("Dendrogram - ITSY1.png", width = 8, height = 4, units = 'in', res=600)
# par(mar=c(3,6,3,4), xpd=T)
# plot(dendro_ITSY1, horiz=FALSE, main= "Fungi Y1", axes=FALSE)
# colored_bars(colors = cbind(col_ITSY1_Lm,col_ITSY1_fac), dend=dendro_ITSY1, rowLabels = c('L. monocytogenes','Facility'), 
#              sort_by_labels_order = FALSE)
# legend("topright", inset = c(-0.1,-0.1), legend = c('-','+','F1','F2','F3'), 
#        fill = c('#F8CD9C','#EA7580','#420A68','#BB3754','#FCA50A'), bty='n')
# dev.off() 
# 
# png("Dendrogram - ITSY2.png", width = 8, height = 4, units = 'in', res=600)
# par(mar=c(3,6,3,4), xpd=T)
# plot(dendro_ITSY2, horiz=FALSE, main= "Fungi Y2 ", axes=FALSE)
# colored_bars(colors = cbind(col_ITSY2_Lm,col_ITSY2_fac), dend=dendro_ITSY2, rowLabels = c('L. monocytogenes','Facility'), 
#              sort_by_labels_order = FALSE)
# legend("topright", inset = c(-0.1,-0.1), legend = c('-','+','F1','F2','F3'), 
#        fill = c('#F8CD9C','#EA7580','#420A68','#BB3754','#FCA50A'), bty='n')
# dev.off()  
# 




# PERMANOVA #
#Using Aitchinson distances

#Two way anova by facility and Lmono presence
#16s-

permanova_16sF1<-pairwise.adonis2(dist_16sF1~Year:Month+Year+Month, data=metadata_16sF1, perm = 999, p.adjust.m = 'bonferroni')
p.adjust(permanova_16sF1$`Y1_vs_Y2`$`Pr(>F)`, method = 'bonferroni')

permanova_16sF2<-pairwise.adonis2(dist_16sF2~Year:Month+Year+Month, data=metadata_16sF2, perm = 999, p.adjust.m = 'bonferroni')
p.adjust(permanova_16sF2$`Y1_vs_Y2`$`Pr(>F)`, method = 'bonferroni')

permanova_16sF3<-pairwise.adonis2(dist_16sF3~Year:Month+Year+Month, data=metadata_16sF3, perm = 999, p.adjust.m = 'bonferroni')
p.adjust(permanova_16sF3$`Y1_vs_Y2`$`Pr(>F)`, method = 'bonferroni')


#ITS
permanova_ITSF1<-pairwise.adonis2(dist_ITSF1~Year:Month+Year+Month, data=metadata_ITSF1, perm = 999, p.adjust.m = 'bonferroni')
p.adjust(permanova_ITSF1$`Y1_vs_Y2`$`Pr(>F)`, method = 'bonferroni')

permanova_ITSF2<-pairwise.adonis2(dist_ITSF2~Year:Month+Year+Month, data=metadata_ITSF2, perm = 999, p.adjust.m = 'bonferroni')
p.adjust(permanova_ITSF2$`Y1_vs_Y2`$`Pr(>F)`, method = 'bonferroni')

permanova_ITSF3<-pairwise.adonis2(dist_ITSF3~Year:Month+Year+Month, data=metadata_ITSF3, perm = 999, p.adjust.m = 'bonferroni')
p.adjust(permanova_ITSF3$`Y1_vs_Y2`$`Pr(>F)`, method = 'bonferroni')




#### Random forest ####
#Adapted from Taejung: https://github.com/tuc289/SurfaceWaterMicrobiome/blob/master/Year2/year2_GSD_poster.R
#Determine if there are taxa (OTU) that can predict the presence of Lm
# Percentage positive Lm samples: both years 65.5%, Y1 56.4% and Y2 75.7%

#Using OTU tables combining Y1 and Y2 for bacteria and fungi
physeq <- phyloseq(otu_table(biom),tax_table(biom),  sample_data)
library(ape)
TREE= rtree(ntaxa(physeq), rooted=TRUE, tip.label = taxa_names(physeq))
physeq = phyloseq(otu_table(biom),tax_table(biom),  sample_data,TREE)
physeq
physeq = transform_sample_counts(physeq, function(x) x/sum(x))
physeq_genus <- taxa_level(physeq, "Rank6")
physeq_family <- taxa_level(physeq, "Rank5")

devtools::install_github("vmikk/metagMisc")
library(metagMisc)

phyloseq_prevalence_plot(physeq_family)

# Identify bacterial families associated with the presence of Salmonella in sediment fractions.
model_sal <- as.data.frame(cbind(otu_table(physeq_family), sample_data(physeq_family)[,11]))
model_sal$salmonlla <- as.factor(model_sal$salmonlla)

# Load Dataset
x <- model_sal[,1:ncol(model_sal)-1]
y <- model_sal[,ncol(model_sal)] #Here, you can identify which colume will be used as classifier

dataset_rf <- as.data.frame(cbind(x,y))

# Random Forest 

library(caret)
library(randomForest)
library(MLeval)

bestMtry <- tuneRF(x,y, mtryStart = 13,stepFactor = 1.5, improve = 1e-5, ntree = 10001)

## mtry set as 67 ##
control <- trainControl(method = 'repeatedcv',
                        number = 10,
                        repeats = 10,
                        search = 'grid', classProbs = TRUE, savePredictions = TRUE)
#create tunegrid
tunegrid <- expand.grid(.mtry = 100)
modellist <- list()

#train with different ntree parameters
for (ntree in c(1001,1501,2001,2501)){
  set.seed(123)
  fit <- train(y~.,
               data = dataset_rf,
               method = 'rf',
               metric = 'Accuracy',
               tuneGrid = tunegrid,
               trControl = control,
               ntree = ntree, )
  key <- toString(ntree)
  modellist[[key]] <- fit
}

#Compare results
results <- resamples(modellist)
summary(results)
dotplot(results)

## ntree set as 1501 ##
set.seed(375)

rf_default <- train(y~., 
                    data=dataset_rf, 
                    method='rf', 
                    metric='Accuracy', 
                    tuneGrid=tunegrid, 
                    trControl=control, ntree = 1501)
print(rf_default)

rf_default$resample
varImp(rf_default, scale = T)


tiff('rf.tiff', units="in", width=5, height=4, res=600, compression = 'lzw')
res <- evalm(rf_default)
dev.off()








#### DIFFERENTIAL ABUNDANCE BY SEASON ####
#Modify OTU tables to replace the row names from OTU identifier to Family name
#Add Year info to Metadata and merge
#Merge OTU tables for both years

#Run Aldex2 by Year


#### NETWORK ANALYSIS ####















# #At the family level
# 
# ### Collapse OTU table to family level
# # Make phyloseq object
# OTU_16s <- otu_table(otus_16s, taxa_are_rows = TRUE)
# phyloseq_16s = phyloseq(OTU_16s, TAX_16s, META_16s)
# TREE_16s = rtree(ntaxa(phyloseq_16s), rooted=TRUE, tip.label = taxa_names(phyloseq_16s))  
# phyloseq16s <- phyloseq(OTU_16s, TAX_16s, TREE_16s, META_16s)
# 
# OTU_16sv138 <- otu_table(otus_16sv138, taxa_are_rows = TRUE)
# phyloseq_16sv138 = phyloseq(OTU_16sv138, TAX_16sv138, META_16s)
# TREE_16sv138 = rtree(ntaxa(phyloseq_16sv138), rooted=TRUE, tip.label = taxa_names(phyloseq_16sv138))  
# phyloseq16sv138 <- phyloseq(OTU_16sv138, TAX_16sv138, TREE_16sv138, META_16s)
# 
# # OTU_ITSY1 <- otu_table(otus_ITSY1, taxa_are_rows = TRUE)
# # phyloseq_ITSY1 = phyloseq(OTU_ITSY1, TAX_ITSY1, META_ITSY1)
# # TREE_ITSY1 = rtree(ntaxa(phyloseq_ITSY1), rooted=TRUE, tip.label = taxa_names(phyloseq_ITSY1))  
# # phyloseqITSY1 <- phyloseq(OTU_ITSY1, TAX_ITSY1, TREE_ITSY1, META_ITSY1)
# # 
# # OTU_ITSY2 <- otu_table(otus_ITSY2, taxa_are_rows = TRUE)
# # phyloseq_ITSY2 = phyloseq(OTU_ITSY2, TAX_ITSY2, META_ITSY2)
# # TREE_ITSY2 = rtree(ntaxa(phyloseq_ITSY2), rooted=TRUE, tip.label = taxa_names(phyloseq_ITSY2))  
# # phyloseqITSY2 <- phyloseq(OTU_ITSY2, TAX_ITSY2, TREE_ITSY2, META_ITSY2)
# 
# #Collapse data to family level
# phyloseq16s_family<-tax_glom(phyloseq16s,taxrank = "Family")
# phyloseq16sv138_family<-tax_glom(phyloseq16sv138,taxrank = "Family")
# # phyloseqITSY1_family<-tax_glom(phyloseqITSY1,taxrank = "Family")
# # phyloseqITSY2_family<-tax_glom(phyloseqITSY2,taxrank = "Family")
# 
# #Extract OTU level at the family level from Phyloseq object
# otus_16s_family<-as.data.frame(otu_table(phyloseq16s_family))
# otus_16sv138_family<-as.data.frame(otu_table(phyloseq16sv138_family))
# # otus_ITSY1_family<-as.data.frame(otu_table(phyloseqITSY1_family))
# # otus_ITSY2_family<-as.data.frame(otu_table(phyloseqITSY2_family))
# 
# 
# 
# #Step 1: Convert OTU table to appropriate format
# #Following step requires samples on rows and OTUs in columns
# head(t(otus_16s_family)) #check that data in the correct format: samples on rows and OTUs in columns
# head(t(otus_16sv138_family))
# # head(t(otus_ITSY1_family))
# # head(t(otus_ITSY2_family))
# 
# #Step 2: Replace zero values before clr transformation. 
# #Uses Count Zero Multiplicative method to replace zeros and outputs pseudo-counts
# otu.n0_16s<-t(cmultRepl(t(otus_16s_family), label=0, method="CZM", output="p-counts")) #49, 593  corrected values
# otu.n0_16sv138<-t(cmultRepl(t(otus_16sv138_family), label=0, method="CZM", output="p-counts")) #52517  corrected values
# # otu.n0_ITSY1<-t(cmultRepl(t(otus_ITSY1_family), label=0, method="CZM", output="p-counts")) #6,789 corrected values
# # otu.n0_ITSY2<-t(cmultRepl(t(otus_ITSY2_family), label=0, method="CZM", output="p-counts")) #5,730 corrected values
# 
# head(otu.n0_16s) #output table needs to have samples in columns and OTUs in rows
# head(otu.n0_16sv138)
# # head(otu.n0_ITSY1)
# # head(otu.n0_ITSY2)
# 
# #Step 3: Convert data to proportions
# otu.n0_16s_prop<-apply(otu.n0_16s, 2, function(x) {x/sum(x)})
# otu.n0_16sv138_prop<-apply(otu.n0_16sv138, 2, function(x) {x/sum(x)})
# # otu.n0_ITSY1_prop<-apply(otu.n0_ITSY1, 2, function(x) {x/sum(x)})
# # otu.n0_ITSY2_prop<-apply(otu.n0_ITSY2, 2, function(x) {x/sum(x)})
# 
# head(otu.n0_16s_prop)
# head(otu.n0_16sv138_prop)
# # head(otu.n0_ITSY1_prop)
# # head(otu.n0_ITSY2_prop)
# 
# #Step 4: Perform abundance and sample filtering and deal sparsity
# #Filter the data to remove all taxa with less than 0.00001% abundance in any sample
# 
# otu.n0_16s_prop_f<-otu.n0_16s[apply(otu.n0_16s_prop, 1, min) > 0.0000001, ]
# otu.n0_16sv138_prop_f<-otu.n0_16sv138[apply(otu.n0_16sv138_prop, 1, min) > 0.0000001, ]
# # otu.n0_ITSY1_prop_f<-otu.n0_ITSY1[apply(otu.n0_ITSY1_prop, 1, min) > 0.0000001, ]
# # otu.n0_ITSY2_prop_f<-otu.n0_ITSY2[apply(otu.n0_ITSY2_prop, 1, min) > 0.0000001, ]
# 
# head(otu.n0_16s_prop_f) #Check that samples are on colums and OTUs in rows
# head(otu.n0_16sv138_prop_f)
# # head(otu.n0_ITSY1_prop_f)
# # head(otu.n0_ITSY2_prop_f)
# 
# #Step 5: perform CLR transformation
# otu.n0.clr_16s<-t(apply(otu.n0_16s_prop_f, 2, function(x){log(x)-mean(log(x))}))
# otu.n0.clr_16sv138<-t(apply(otu.n0_16sv138_prop_f, 2, function(x){log(x)-mean(log(x))}))
# # otu.n0.clr_ITSY1<-t(apply(otu.n0_ITSY1_prop_f, 2, function(x){log(x)-mean(log(x))}))
# # otu.n0.clr_ITSY2<-t(apply(otu.n0_ITSY2_prop_f, 2, function(x){log(x)-mean(log(x))}))
# 
# head(otu.n0.clr_16s) #Check output table. Samples should be in rows and OTUs in columns
# head(otu.n0.clr_16sv138)
# # head(otu.n0.clr_ITSY1)
# # head(otu.n0.clr_ITSY2)
# 
# #Step 6: Perform Singular Value Decomposition (PCA)
# pc.clr_16s<-prcomp(otu.n0.clr_16s)
# pc.clr_16sv138<-prcomp(otu.n0.clr_16sv138)
# # pc.clr_ITSY1<-prcomp(otu.n0.clr_ITSY1)
# # pc.clr_ITSY2<-prcomp(otu.n0.clr_ITSY2)
# 
# png("Screeplot - PCA - by family.png", width = 5, height = 5, units = 'in', dpi=600)
# par(mar=c(5,5,2,2))
# par(mfrow=c(2,2))
# screeplot(pc.clr_16s, type='lines', main="Bacteria - ")
# screeplot(pc.clr_16sv138, type='lines', main="Bacteria - v138")
# screeplot(pc.clr_ITSY1, type='lines', main="Fungi - Y1")
# screeplot(pc.clr_ITSY2, type='lines', main="Fungi - Y2")
# dev.off()
# 
# #Calculate total variance of the data
# mvar.clr_16s<-mvar(otu.n0.clr_16s)
# mvar.clr_16s<-mvar(otu.n0.clr_16sv138)
# # mvar.clr_ITSY1<-mvar(otu.n0.clr_ITSY1)
# # mvar.clr_ITSY2<-mvar(otu.n0.clr_ITSY2)
# 
# 
# 
# #Display results - 16s v138
# row_16sv138<-rownames(otu.n0.clr_16sv138) #Make vector with sample names
# pc_out_16sv138<-as.data.frame(pc.clr_16sv138$x[,1:3]) #Get PC1 and PC2 
# pc_out_meta_16sv138<-as.data.frame(bind_cols(pc_out_16sv138,metadata_16s)) #Add metadata information
# row.names(pc_out_meta_16sv138)<-row_16sv138 #Add rownames to dataframe
# pc_out_meta_16sv138$Facility<-as.factor(pc_out_meta_16sv138$Facility)
# pc_out_meta_16sv138$L..monocytogenes<-as.factor(pc_out_meta_16sv138$L..monocytogenes)
# 
# # Make PCA plot - First 2 axis- color by facility and shape by L..mono
# PCA_16sv138_Lmono <- ggplot(pc_out_meta_16sv138, aes(x=PC1,y=PC2, color=Facility, shape=L..monocytogenes))+
#   geom_point(size=3)+facet_grid(.~Facility)+
#   theme(legend.text=element_text(size=13,color='black'), legend.title= element_text(size=13, face='italic')) +
#   theme(legend.position = 'bottom')+
#   theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
#   theme(axis.text = element_blank(), axis.ticks = element_blank()) +
#   theme(axis.title = element_text(size=13,color='black')) +
#   theme(panel.background = element_rect(fill='transparent', color = NA),
#         plot.background = element_rect(fill = 'transparent',color = NA), 
#         panel.border = element_rect(color='black',fill = NA,size=1))+
#   scale_x_continuous(name = paste("PC1: ", round(pc.clr_16sv138$sdev[1]^2/mvar.clr_16sv138*100, digits=1), "%", sep="")) +
#   scale_y_continuous(name = paste("PC2: ", round(pc.clr_16sv138$sdev[2]^2/mvar.clr_16sv138*100, digits=1), "%", sep="")) +
#   ggtitle("Bacteria - v138")+theme(plot.title = element_text(hjust = 0.5, face = 'bold'))+
#   scale_color_viridis_d(begin = 0.2, end = 0.8, option='inferno')
# PCA_16sv138_Lmono
# 
# # Make PCA plot - First 2 axis- color by facility and shape by Year
# PCA_16sv138_FacY <- ggplot(pc_out_meta_16sv138, aes(x=PC1,y=PC2, color=Facility, shape=Year))+
#   geom_point(size=3)+facet_grid(.~Facility)+
#   theme(legend.text=element_text(size=13,color='black'), legend.title= element_text(size=13, face='italic')) +
#   theme(legend.position = 'bottom')+
#   theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
#   theme(axis.text = element_blank(), axis.ticks = element_blank()) +
#   theme(axis.title = element_text(size=13,color='black')) +
#   theme(panel.background = element_rect(fill='transparent', color = NA),
#         plot.background = element_rect(fill = 'transparent',color = NA), 
#         panel.border = element_rect(color='black',fill = NA,size=1))+
#   scale_x_continuous(name = paste("PC1: ", round(pc.clr_16sv138$sdev[1]^2/mvar.clr_16sv138*100, digits=1), "%", sep="")) +
#   scale_y_continuous(name = paste("PC2: ", round(pc.clr_16sv138$sdev[2]^2/mvar.clr_16sv138*100, digits=1), "%", sep="")) +
#   ggtitle("Bacteria - v138")+theme(plot.title = element_text(hjust = 0.5, face = 'bold'))+
#   scale_color_viridis_d(begin = 0.2, end = 0.8, option='inferno')
# PCA_16sv138_FacY
# 
# # Make PCA plot - First 2 axis- color by week and shape by FacY
# PCA_16sv138_FacY <- ggplot(pc_out_meta_16sv138, aes(x=PC1,y=PC2, color=Month, shape=Year))+
#   geom_point(size=3)+facet_grid(.~Facility)+
#   theme(legend.text=element_text(size=13,color='black'), legend.title= element_text(size=13, face='italic')) +
#   theme(legend.position = 'bottom')+
#   theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
#   theme(axis.text = element_blank(), axis.ticks = element_blank()) +
#   theme(axis.title = element_text(size=13,color='black')) +
#   theme(panel.background = element_rect(fill='transparent', color = NA),
#         plot.background = element_rect(fill = 'transparent',color = NA), 
#         panel.border = element_rect(color='black',fill = NA,size=1))+
#   scale_x_continuous(name = paste("PC1: ", round(pc.clr_16sv138$sdev[1]^2/mvar.clr_16sv138*100, digits=1), "%", sep="")) +
#   scale_y_continuous(name = paste("PC2: ", round(pc.clr_16sv138$sdev[2]^2/mvar.clr_16sv138*100, digits=1), "%", sep="")) +
#   ggtitle("Bacteria - v138")+theme(plot.title = element_text(hjust = 0.5, face = 'bold'))+
#   scale_color_viridis_d(begin = 0.2, end = 0.8, option='inferno')
# PCA_16sv138_FacY
# 
# #Combine plots
# PCA_all = plot_grid(PCA +theme(legend.position = "none"),
#                             bubbleplot_16sY2+theme(legend.position = "none"), 
#                             bubbleplot_ITSY1+theme(legend.position = "none"), 
#                             bubbleplot_ITSY2 +theme(legend.position = "none"),
#                             ncol=2, nrow=2, labels = c("A","B","C","D"), label_size = 20, vjust = 2, hjust = -1.5)
# Bubbleplots_all
# 
# ggsave("PC_16sv138.png", plot =PCA_16sv138, device="png", width=6, height=6, units="in",dpi=600)


# #Plot with third dimension as size of point. Seems to crowded to me
# PCA_16sY1_three <- ggplot(pc_out_meta_16sY1, aes(x=PC1,y=PC2, color=Facility, shape=L..monocytogenes))+
#   geom_point(aes(size=PC3))+
#   theme(legend.text=element_text(size=13,color='black'), legend.title= element_text(size=13, face='italic')) +
#   theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
#   theme(axis.title = element_text(size=13,color='black')) +
#   theme(panel.background = element_rect(fill='transparent', color = NA),
#         plot.background = element_rect(fill = 'transparent',color = NA), 
#         panel.border = element_rect(color='black',fill = NA,size=1))+
#   scale_x_continuous(name = paste("PC1: ", round(pc.clr_16sY1$sdev[1]^2/mvar.clr_16sY1*100, digits=1), "%", sep="")) +
#   scale_y_continuous(name = paste("PC2: ", round(pc.clr_16sY1$sdev[2]^2/mvar.clr_16sY1*100, digits=1), "%", sep="")) +
#   ggtitle("Bacteria - Y1")+theme(plot.title = element_text(hjust = 0.5, face = 'bold'))+
#   scale_color_viridis_d(begin = 0.2, end = 0.8, option='inferno')
# PCA_16sY1_three
# ggsave("PC_16sY1_three.png", plot =PCA_16sY1_three, device="png", width=8, height=6, units="in",dpi=600)
# 
# 
# 
# #Display results - ITSY1
# row_ITSY1<-rownames(otu.n0.clr_ITSY1) #Make vector with sample names
# pc_out_ITSY1<-as.data.frame(pc.clr_ITSY1$x[,1:3]) #Get PC1 and PC2 
# pc_out_meta_ITSY1<-as.data.frame(bind_cols(pc_out_ITSY1,metadata_ITSY1)) #Add metadata information
# row.names(pc_out_meta_ITSY1)<-row_ITSY1 #Add rownames to dataframe
# pc_out_meta_ITSY1$Facility<-as.factor(pc_out_meta_ITSY1$Facility)
# pc_out_meta_ITSY1$L..monocytogenes<-as.factor(pc_out_meta_ITSY1$L..monocytogenes)
# 
# 
# PCA_ITSY1 <- ggplot(pc_out_meta_ITSY1, aes(x=PC1,y=PC2, color=Facility, shape=L..monocytogenes))+
#   geom_point(size=3)+
#   theme(legend.text=element_text(size=13,color='black'), legend.title= element_text(size=13, face='italic')) +
#   theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
#   theme(legend.position = 'bottom')+
#   theme(axis.text = element_blank(), axis.ticks = element_blank()) +
#   theme(axis.title = element_text(size=13,color='black')) +
#   theme(panel.background = element_rect(fill='transparent', color = NA),
#         plot.background = element_rect(fill = 'transparent',color = NA), 
#         panel.border = element_rect(color='black',fill = NA,size=1))+
#   scale_x_continuous(name = paste("PC1: ", round(pc.clr_ITSY1$sdev[1]^2/mvar.clr_ITSY1*100, digits=1), "%", sep="")) +
#   scale_y_continuous(name = paste("PC2: ", round(pc.clr_ITSY1$sdev[2]^2/mvar.clr_ITSY1*100, digits=1), "%", sep=""))+
#   ggtitle("Fungi - Y1")+theme(plot.title = element_text(hjust = 0.5, face = 'bold'))+
#   scale_color_viridis_d(begin = 0.2, end = 0.8, option='inferno')
# PCA_ITSY1
# ggsave("PC_ITSY1_Family.png", plot =PCA_ITSY1, device="png", width=6, height=6, units="in",dpi=600)
# 
# #Display results - 16sY2
# row_16sY2<-rownames(otu.n0.clr_16sY2) #Make vector with sample names
# pc_out_16sY2<-as.data.frame(pc.clr_16sY2$x[,1:32]) #Get PC1 and PC2 
# pc_out_meta_16sY2<-as.data.frame(bind_cols(pc_out_16sY2,metadata_16sY2)) #Add metadata information
# row.names(pc_out_meta_16sY2)<-row_16sY2 #Add rownames to dataframe
# pc_out_meta_16sY2$Facility<-as.factor(pc_out_meta_16sY2$Facility)
# pc_out_meta_16sY2$L..monocytogenes<-as.factor(pc_out_meta_16sY2$L..monocytogenes)
# 
# # Make PCA plot (Adapted from reference to use ggplot2 features)
# PCA_16sY2 <- ggplot(pc_out_meta_16sY2, aes(x=PC1,y=PC2, color=Facility, shape=L..monocytogenes))+
#   geom_point(size=3)+
#   theme(legend.text=element_text(size=13,color='black'), legend.title= element_text(size=13, face='italic')) +
#   theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
#   theme(axis.text = element_blank(), axis.ticks = element_blank()) +
#   theme(axis.title = element_text(size=13,color='black')) +
#   theme(legend.position = 'bottom')+
#   theme(panel.background = element_rect(fill='transparent', color = NA),
#         plot.background = element_rect(fill = 'transparent',color = NA), 
#         panel.border = element_rect(color='black',fill = NA,size=1))+
#   scale_x_continuous(name = paste("PC1: ", round(pc.clr_16sY2$sdev[1]^2/mvar.clr_16sY2*100, digits=1), "%", sep="")) +
#   scale_y_continuous(name = paste("PC2: ", round(pc.clr_16sY2$sdev[2]^2/mvar.clr_16sY2*100, digits=1), "%", sep="")) +
#   ggtitle("Bacteria - Y2")+theme(plot.title = element_text(hjust = 0.5, face = 'bold'))+
#   scale_color_viridis_d(begin = 0.2, end = 0.8, option='inferno')
# PCA_16sY2
# ggsave("PC_16sY2_Famiily.png", plot =PCA_16sY2, device="png", width=6, height=6, units="in",dpi=600)
# 
# 
# #Display results - ITSY2
# row_ITSY2<-rownames(otu.n0.clr_ITSY2) #Make vector with sample names
# pc_out_ITSY2<-as.data.frame(pc.clr_ITSY2$x[,1:3]) #Get PC1 and PC2 
# pc_out_meta_ITSY2<-as.data.frame(bind_cols(pc_out_ITSY2,metadata_ITSY2)) #Add metadata information
# row.names(pc_out_meta_ITSY2)<-row_ITSY2 #Add rownames to dataframe
# pc_out_meta_ITSY2$Facility<-as.factor(pc_out_meta_ITSY2$Facility)
# pc_out_meta_ITSY2$L..monocytogenes<-as.factor(pc_out_meta_ITSY2$L..monocytogenes)
# 
# # Make PCA plot (Adapted from reference to use ggplot2 features)
# PCA_ITSY2 <- ggplot(pc_out_meta_ITSY2, aes(x=PC1,y=PC2, color=Facility, shape=L..monocytogenes))+
#   geom_point(size=3)+
#   theme(legend.text=element_text(size=13,color='black'), legend.title= element_text(size=13, face='italic')) +
#   theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
#   theme(axis.text = element_blank(), axis.ticks = element_blank()) +
#   theme(axis.title = element_text(size=13,color='black')) +
#   theme(legend.position = 'bottom')+
#   theme(panel.background = element_rect(fill='transparent', color = NA),
#         plot.background = element_rect(fill = 'transparent',color = NA), 
#         panel.border = element_rect(color='black',fill = NA,size=1))+
#   scale_x_continuous(name = paste("PC1: ", round(pc.clr_ITSY2$sdev[1]^2/mvar.clr_ITSY2*100, digits=1), "%", sep="")) +
#   scale_y_continuous(name = paste("PC2: ", round(pc.clr_ITSY2$sdev[2]^2/mvar.clr_ITSY2*100, digits=1), "%", sep=""))+
#   ggtitle("Fungi - Y2")+theme(plot.title = element_text(hjust = 0.5, face = 'bold'))+
#   scale_color_viridis_d(begin = 0.2, end = 0.8, option='inferno')
# PCA_ITSY2
# ggsave("PC_ITSY2_Family.png", plot =PCA_ITSY2, device="png", width=6, height=6, units="in",dpi=600)
# 
# 
# #Combine PCA plots
# PCA_CoDa_a = plot_grid(PCA_16sY1+ theme(legend.position = "none"),
#                          PCA_16sY2+ theme(legend.position = "none"),
#                          PCA_ITSY1+ theme(legend.position = "none") ,
#                          PCA_ITSY2+ theme(legend.position = "none") ,
#                          ncol=2, nrow=2, labels = c("A","B","C","D"), label_size = 20, vjust = 3, hjust = -2)
# PCA_CoDa_b =get_legend(PCA_16sY1)
# PCA_CoDA=plot_grid(PCA_CoDa_a,PCA_CoDa_b, nrow = 2, rel_heights =  c(5,1))
# PCA_CoDA
# 
# ggsave("PCA_Family_Y1Y2.png", plot =PCA_CoDA, device="png", width=8, height=10, units="in",dpi=600)


##First 3 axis (3D plot)
#16sY1
#Get PCA axis labels
# paste("PC1: ", round(pc.clr_16sY1$sdev[1]^2/mvar.clr_16sY1*100, digits=1), "%", sep="")
# paste("PC2: ", round(pc.clr_16sY1$sdev[2]^2/mvar.clr_16sY1*100, digits=1), "%", sep="")
# paste("PC3: ", round(pc.clr_16sY1$sdev[3]^2/mvar.clr_16sY1*100, digits=1), "%", sep="")
# 
# PCA_16sY1_3D<-plot_ly(pc_out_meta_16sY1, x = ~PC1, y = ~PC2, z = ~PC3, color = ~Facility, 
#                       colors = c("#420A68","#BB3754","#FCA50A"), hoverinfo = 'text',
#                       text = ~paste('</br> Facility: ', Facility,
#                                     '</br> Section: ', Section,
#                                     '</br> Week: ', Week,
#                                     '</br> L.monocytogenes: ', L..monocytogenes)) %>% 
#   add_markers(marker = list(size= 5)) %>% 
#   layout(title = 'Bacteria - Y1', scene = list(xaxis = list(title = "PC1: 23.5%"),
#                                                yaxis = list(title = "PC2: 12.1%"),
#                                                zaxis = list(title = "PC3: 10.0%")))
# 
# # save as html to allow rotation of axis
# saveWidget(PCA_16sY1_3D, file="PCA-16sY1 3D.html")
# 
# #16sY2
# #Get PCA axis labels
# paste("PC1: ", round(pc.clr_16sY2$sdev[1]^2/mvar.clr_16sY2*100, digits=1), "%", sep="")
# paste("PC2: ", round(pc.clr_16sY2$sdev[2]^2/mvar.clr_16sY2*100, digits=1), "%", sep="")
# paste("PC3: ", round(pc.clr_16sY2$sdev[3]^2/mvar.clr_16sY2*100, digits=1), "%", sep="")
# 
# PCA_16sY2_3D<-plot_ly(pc_out_meta_16sY2, x = ~PC1, y = ~PC2, z = ~PC3, color = ~Facility, 
#                       colors = c("#420A68","#BB3754","#FCA50A"), hoverinfo = 'text',
#                       text = ~paste('</br> Facility: ', Facility,
#                                     '</br> Section: ', Section,
#                                     '</br> Week: ', Week,
#                                     '</br> L.monocytogenes: ', L..monocytogenes)) %>% 
#   add_markers(marker = list(size= 5)) %>% 
#   layout(title = 'Bacteria - Y2', scene = list(xaxis = list(title = "PC1: 14.9%"),
#                                                yaxis = list(title = "PC2: 11.2%"),
#                                                zaxis = list(title = "PC3: 7.4%")))
# 
# # save as html to allow rotation of axis
# saveWidget(PCA_16sY2_3D, file="PCA-16sY2 3D.html")
# 
# #ITSY1
# #Get PCA axis labels
# paste("PC1: ", round(pc.clr_ITSY1$sdev[1]^2/mvar.clr_ITSY1*100, digits=1), "%", sep="")
# paste("PC2: ", round(pc.clr_ITSY1$sdev[2]^2/mvar.clr_ITSY1*100, digits=1), "%", sep="")
# paste("PC3: ", round(pc.clr_ITSY1$sdev[3]^2/mvar.clr_ITSY1*100, digits=1), "%", sep="")
# 
# PCA_ITSY1_3D<-plot_ly(pc_out_meta_ITSY1, x = ~PC1, y = ~PC2, z = ~PC3, color = ~Facility, 
#                       colors = c("#420A68","#BB3754","#FCA50A"), hoverinfo = 'text',
#                       text = ~paste('</br> Facility: ', Facility,
#                                     '</br> Section: ', Section,
#                                     '</br> Week: ', Week,
#                                     '</br> L.monocytogenes: ', L..monocytogenes)) %>% 
#   add_markers(marker = list(size= 5)) %>% 
#   layout(title = 'Fungi - Y1', scene = list(xaxis = list(title = "PC1: 18.6%"),
#                                                yaxis = list(title = "PC2: 11.6%"),
#                                                zaxis = list(title = "PC3: 6.9%")))
# 
# # save as html to allow rotation of axis
# saveWidget(PCA_ITSY1_3D, file="PCA-ITSY1 3D.html")
# 
# #ITSY2
# #Get PCA axis labels
# paste("PC1: ", round(pc.clr_ITSY2$sdev[1]^2/mvar.clr_ITSY2*100, digits=1), "%", sep="")
# paste("PC2: ", round(pc.clr_ITSY2$sdev[2]^2/mvar.clr_ITSY2*100, digits=1), "%", sep="")
# paste("PC3: ", round(pc.clr_ITSY2$sdev[3]^2/mvar.clr_ITSY2*100, digits=1), "%", sep="")
# 
# PCA_ITSY2_3D<-plot_ly(pc_out_meta_ITSY2, x = ~PC1, y = ~PC2, z = ~PC3, color = ~Facility, 
#                       colors = c("#420A68","#BB3754","#FCA50A"), hoverinfo = 'text',
#                       text = ~paste('</br> Facility: ', Facility,
#                                     '</br> Section: ', Section,
#                                     '</br> Week: ', Week,
#                                     '</br> L.monocytogenes: ', L..monocytogenes)) %>% 
#   add_markers(marker = list(size= 5)) %>% 
#   layout(title = 'Bacteria - Y2', scene = list(xaxis = list(title = "PC1: 15.8%"),
#                                                yaxis = list(title = "PC2: 9.5%"),
#                                                zaxis = list(title = "PC3: 7.3%")))
# 
# # save as html to allow rotation of axis
# saveWidget(PCA_ITSY2_3D, file="PCA-ITSY2 3D.html")

# 
# 
# #### HIERARCHICAL CLUSTERING ####
# #Generate the distance matrix
# dist_16s<-dist(otu.n0.clr_16s, method='euclidean') #Distance matrix based on Aitchinson simplex
# dist_16sv138<-dist(otu.n0.clr_16sv138, method='euclidean')
# # dist_ITSY1<-dist(otu.n0.clr_ITSY1, method='euclidean')
# # dist_ITSY2<-dist(otu.n0.clr_ITSY2, method='euclidean')
# 
# #Cluster the data
# hc_16s<-hclust(dist_16s, method = 'ward.D2')
# hc_16sv138<-hclust(dist_16sv138, method = 'ward.D2')
# # hc_ITSY1<-hclust(dist_ITSY1, method = 'ward.D2')
# # hc_ITSY2<-hclust(dist_ITSY2, method = 'ward.D2')
# 
# #Plot the dendrogram
# #Set group color by Facility
# metadata_16s$Facility<-as.factor(metadata_16s$Facility) #Make Facility a factor
# metadata_16s$Facility<-as.factor(metadata_16s$Facility)
# # metadata_ITSY1$Facility<-as.factor(metadata_ITSY1$Facility)
# # metadata_ITSY2$Facility<-as.factor(metadata_ITSY2$Facility)
# 
# metadata_16s$L..monocytogenes<-as.factor(metadata_16s$L..monocytogenes) #make L.mono a factor
# metadata_16s$L..monocytogenes<-as.factor(metadata_16s$L..monocytogenes)
# # metadata_ITSY1$L..monocytogenes<-as.factor(metadata_ITSY1$L..monocytogenes) 
# # metadata_ITSY2$L..monocytogenes<-as.factor(metadata_ITSY2$L..monocytogenes)
# 
# metadata_16s$Year<-as.factor(metadata_16s$Year) #make Year a factor
# metadata_16s$Year<-as.factor(metadata_16s$Year)
# # metadata_ITSY1$L..monocytogenes<-as.factor(metadata_ITSY1$L..monocytogenes) 
# # metadata_ITSY2$L..monocytogenes<-as.factor(metadata_ITSY2$L..monocytogenes)
# 
# #Make vector with colors y facility and order them as in dendrogram
# cols_fac <- c('#420A68','#BB3754','#FCA50A') #colors for facilities
# col_16s_fac <- (cols_fac[metadata_16s$Facility])[order.dendrogram(as.dendrogram(hc_16s))]
# col_16sv138_fac <- (cols_fac[metadata_16s$Facility])[order.dendrogram(as.dendrogram(hc_16sv138))]
# # col_ITSY1_fac <- (cols_fac[metadata_ITSY1$Facility])[order.dendrogram(as.dendrogram(hc_ITSY1))]
# # col_ITSY2_fac <- (cols_fac[metadata_ITSY2$Facility])[order.dendrogram(as.dendrogram(hc_ITSY2))]
# 
# 
# #Make vector with colors y facility and order them as in dendrogram
# cols_Year <- c('#F8CD9C','#EA7580') #colors for year
# col_16s_Year <- (cols_Year[metadata_16s$Year])[order.dendrogram(as.dendrogram(hc_16s))]
# col_16sv138_Year <- (cols_Year[metadata_16s$Year])[order.dendrogram(as.dendrogram(hc_16sv138))]
# # col_ITSY1_Lm <- (cols_Lm[metadata_ITSY1$L..monocytogenes])[order.dendrogram(as.dendrogram(hc_ITSY1))]
# # col_ITSY2_Lm <- (cols_Lm[metadata_ITSY2$L..monocytogenes])[order.dendrogram(as.dendrogram(hc_ITSY2))]
# 
# #Make dendrogram object
# dendro_16s<- as.dendrogram(hc_16s) %>%
#   set("labels", "")
# dendro_16sv138<- as.dendrogram(hc_16sv138) %>%
#   set("labels", "")
# # dendro_ITSY1<- as.dendrogram(hc_ITSY1) %>%
# #   set("labels", "")
# # dendro_ITSY2<- as.dendrogram(hc_ITSY2) %>%
# #   set("labels", "")
# 
# #plot and save
# png("Dendrogram - 16sY1.png", width = 6, height = 4, units = 'in', res=600)
# par(mar=c(3,6,3,4), xpd=T)
# plot(dendro_16s, horiz=FALSE, main= "Bacteria - Y1", axes=FALSE)
# colored_bars(colors = cbind(col_16s_fac,col_16s_Year), dend=dendro_16s, rowLabels = c('Facility','Year'), 
#                                                                                           sort_by_labels_order = FALSE)
# legend("topright", inset = c(-0.1,-0.1), legend = c('F1','F2','F3','Y1','Y2'), fill = c('#420A68','#BB3754','#FCA50A','#F8CD9C','#EA7580'), bty='n')
# dev.off() 
# 
# png("Dendrogram - 16sY2.png", width = 6, height = 4, units = 'in', res=600)
# par(mar=c(3,6,3,4), xpd=T)
# plot(dendro_16sY2, horiz=FALSE, main= "Bacteria - Y2", axes=FALSE)
# colored_bars(colors = cbind(col_16sY2_fac,col_16sY2_Lm), dend=dendro_16sY2, rowLabels = c('Facility','L.monocytogenes'), 
#              sort_by_labels_order = FALSE)
# legend("topright", inset = c(-0.1,-0.1), legend = c('F1','F2','F3','-','+'), fill = c('#420A68','#BB3754','#FCA50A','#F8CD9C','#EA7580'), bty='n')
# dev.off()
# 
# png("Dendrogram - ITSY1.png", width = 6, height = 4, units = 'in', res=600)
# par(mar=c(3,6,3,4), xpd=T)
# plot(dendro_ITSY1, horiz=FALSE, main= "Fungi - Y1", axes=FALSE)
# colored_bars(colors = cbind(col_ITSY1_fac,col_ITSY1_Lm), dend=dendro_ITSY1, rowLabels = c('Facility','L.monocytogenes'), 
#              sort_by_labels_order = FALSE)
# legend("topright", inset = c(-0.1,-0.1), legend = c('F1','F2','F3','-','+'), fill = c('#420A68','#BB3754','#FCA50A','#F8CD9C','#EA7580'), bty='n')
# dev.off() 
# 
# png("Dendrogram - ITSY2.png", width = 6, height = 4, units = 'in', res=600)
# par(mar=c(3,6,3,4), xpd=T)
# plot(dendro_ITSY2, horiz=FALSE, main= "Fungi - Y2", axes=FALSE)
# colored_bars(colors = cbind(col_ITSY2_fac,col_ITSY2_Lm), dend=dendro_ITSY2, rowLabels = c('Facility','L.monocytogenes'), 
#              sort_by_labels_order = FALSE)
# legend("topright", inset = c(-0.1,-0.1), legend = c('F1','F2','F3','-','+'), fill = c('#420A68','#BB3754','#FCA50A','#F8CD9C','#EA7580'), bty='n')
# dev.off()
# 
# 
# 
# #### PERMANOVA ####
# #Using Aitchinson distances calculated previously for the dendrogram
# 
# #16sY1
# permanova_16sY1_Fac<-pairwise.adonis(dist_16sY1, factors=metadata_16sY1$Facility, perm = 999, p.adjust.m = 'bonferroni')
# permanova_16sY1_Sec<-pairwise.adonis(dist_16sY1, factors=metadata_16sY1$Section, perm = 999, p.adjust.m = 'bonferroni')
# permanova_16sY1_Month<-pairwise.adonis(dist_16sY1, factors=metadata_16sY1$Month, perm = 999, p.adjust.m = 'bonferroni')
# permanova_16sY1_Lm<-pairwise.adonis2(dist_16sY1~L..monocytogenes:Facility+Facility+L..monocytogenes, data=metadata_16sY1, perm = 999, p.adjust.m = 'bonferroni')
# p.adjust(permanova_16sY1_Lm$`+_vs_-`$`Pr(>F)`, method = 'bonferroni')
# 
# permanova_16sY1<-bind_rows(permanova_16sY1_Fac, permanova_16sY1_Sec, permanova_16sY1_Month)
# write.csv(permanova_16sY1, file='permanova_16sY1_family.csv')
# 
# #16sY2
# permanova_16sY2_Fac<-pairwise.adonis(dist_16sY2, factors=metadata_16sY2$Facility, perm = 999, p.adjust.m = 'bonferroni')
# permanova_16sY2_Sec<-pairwise.adonis(dist_16sY2, factors=metadata_16sY2$Section, perm = 999, p.adjust.m = 'bonferroni')
# permanova_16sY2_Month<-pairwise.adonis(dist_16sY2, factors=metadata_16sY2$Month, perm = 999, p.adjust.m = 'bonferroni')
# permanova_16sY2_Lm<-pairwise.adonis2(dist_16sY2~L..monocytogenes:Facility+Facility+L..monocytogenes, data=metadata_16sY2, perm = 999, p.adjust.m = 'bonferroni')
# p.adjust(permanova_16sY2_Lm$`+_vs_-`$`Pr(>F)`, method = 'bonferroni') # 0.003 1.000 0.753
# 
# permanova_16sY2<-bind_rows(permanova_16sY2_Fac, permanova_16sY2_Sec, permanova_16sY2_Month)
# write.csv(permanova_16sY2, file='permanova_16sY2_family.csv')
# 
# #ITSY1
# permanova_ITSY1_Fac<-pairwise.adonis(dist_ITSY1, factors=metadata_ITSY1$Facility, perm = 999, p.adjust.m = 'bonferroni')
# permanova_ITSY1_Sec<-pairwise.adonis(dist_ITSY1, factors=metadata_ITSY1$Section, perm = 999, p.adjust.m = 'bonferroni')
# permanova_ITSY1_Month<-pairwise.adonis(dist_ITSY1, factors=metadata_ITSY1$Month, perm = 999, p.adjust.m = 'bonferroni')
# permanova_ITSY1_Lm<-pairwise.adonis2(dist_ITSY1~L..monocytogenes:Facility+Facility+L..monocytogenes, data=metadata_ITSY1, perm = 999, p.adjust.m = 'bonferroni')
# permanova_ITSY1_Lm
# p.adjust(permanova_ITSY1_Lm$`+_vs_-`$`Pr(>F)`, method = 'bonferroni') # 0.003 1.000 1.000
# 
# permanova_ITSY1<-bind_rows(permanova_ITSY1_Fac, permanova_ITSY1_Sec, permanova_ITSY1_Month)
# write.csv(permanova_ITSY1, file='permanova_ITSY1_family.csv')
# 
# #ITSY2
# permanova_ITSY2_Fac<-pairwise.adonis(dist_ITSY2, factors=metadata_ITSY2$Facility, perm = 999, p.adjust.m = 'bonferroni')
# permanova_ITSY2_Sec<-pairwise.adonis(dist_ITSY2, factors=metadata_ITSY2$Section, perm = 999, p.adjust.m = 'bonferroni')
# permanova_ITSY2_Month<-pairwise.adonis(dist_ITSY2, factors=metadata_ITSY2$Month, perm = 999, p.adjust.m = 'bonferroni')
# permanova_ITSY2_Lm<-pairwise.adonis2(dist_ITSY2~L..monocytogenes:Facility+Facility+L..monocytogenes, data=metadata_ITSY2, perm = 999, p.adjust.m = 'bonferroni')
# permanova_ITSY2_Lm
# p.adjust(permanova_ITSY2_Lm$`+_vs_-`$`Pr(>F)`, method = 'bonferroni') #0.003 1.000 0.198 
# 
# permanova_ITSY2<-bind_rows(permanova_ITSY2_Fac, permanova_ITSY2_Sec, permanova_ITSY2_Month)
# write.csv(permanova_ITSY2, file='permanova_ITSY2_family.csv')
# 



#### MICROBIOTA AND MYCOBIOTA COMPOSITION PLOTS ####
#Note: used compositional approach to transform the sample counts to compositions. 

#Transform sample counts into compositions
otu.n0.acomp_16sY1<-as.data.frame(acomp(t(otu.n0_16sY1)), total=1)
otu.n0.acomp_16sY2<-as.data.frame(acomp(t(otu.n0_16sY2)), total=1)
otu.n0.acomp_ITSY1<-as.data.frame(acomp(t(otu.n0_ITSY1)), total=1)
otu.n0.acomp_ITSY2<-as.data.frame(acomp(t(otu.n0_ITSY2)), total=1)
rowSums(otu.n0.acomp_16sY1)
rowSums(otu.n0.acomp_16sY2)
rowSums(otu.n0.acomp_ITSY1)
rowSums(otu.n0.acomp_ITSY2)

#Make Phyloseq object
OTU_16sY1_family_prop <- otu_table(otu.n0.acomp_16sY1, taxa_are_rows = FALSE)
phyloseq_16sY1_family_prop = phyloseq(OTU_16sY1_family_prop, TAX_16sY1, META_16sY1)
TREE_16sY1_family_prop = rtree(ntaxa(phyloseq_16sY1_family_prop), rooted=TRUE, tip.label = taxa_names(phyloseq_16sY1_family_prop))  
phyloseq16sY1_family_prop <- phyloseq(OTU_16sY1_family_prop, TAX_16sY1, TREE_16sY1, META_16sY1)

OTU_16sY2_family_prop <- otu_table(otu.n0.acomp_16sY2, taxa_are_rows = FALSE)
phyloseq_16sY2_family_prop = phyloseq(OTU_16sY2_family_prop, TAX_16sY2, META_16sY2)
TREE_16sY2_family_prop = rtree(ntaxa(phyloseq_16sY2_family_prop), rooted=TRUE, tip.label = taxa_names(phyloseq_16sY2_family_prop))  
phyloseq16sY2_family_prop <- phyloseq(OTU_16sY2_family_prop, TAX_16sY2, TREE_16sY2, META_16sY2)

OTU_ITSY1_family_prop <- otu_table(otu.n0.acomp_ITSY1, taxa_are_rows = FALSE)
phyloseq_ITSY1_family_prop = phyloseq(OTU_ITSY1_family_prop, TAX_ITSY1, META_ITSY1)
TREE_ITSY1_family_prop = rtree(ntaxa(phyloseq_ITSY1_family_prop), rooted=TRUE, tip.label = taxa_names(phyloseq_ITSY1_family_prop))  
phyloseqITSY1_family_prop <- phyloseq(OTU_ITSY1_family_prop, TAX_ITSY1, TREE_ITSY1, META_ITSY1)

OTU_ITSY2_family_prop <- otu_table(otu.n0.acomp_ITSY2, taxa_are_rows = FALSE)
phyloseq_ITSY2_family_prop = phyloseq(OTU_ITSY2_family_prop, TAX_ITSY2, META_ITSY2)
TREE_ITSY2_family_prop = rtree(ntaxa(phyloseq_ITSY2_family_prop), rooted=TRUE, tip.label = taxa_names(phyloseq_ITSY2_family_prop))  
phyloseqITSY2_family_prop <- phyloseq(OTU_ITSY2_family_prop, TAX_ITSY2, TREE_ITSY2, META_ITSY2)


#Make long format table from Phyloseq object
family_16sY1 <- phyloseq16sY1_family_prop %>%  
  transform_sample_counts(function(x) {x * 100} ) %>% #Recalculates the relative abundance 
  psmelt() %>%  #Melts to long format
  arrange(desc(Abundance)) #Organizes the table by Abundance in descending order

family_16sY2 <- phyloseq16sY2_family_prop %>%  
  transform_sample_counts(function(x) {x * 100} ) %>% #Recalculates the relative abundance 
  psmelt() %>%  #Melts to long format
  arrange(desc(Abundance)) #Organizes the table by Abundance in descending order

family_ITSY1 <- phyloseqITSY1_family_prop %>%  
  transform_sample_counts(function(x) {x * 100} ) %>% #Recalculates the relative abundance 
  psmelt() %>%  #Melts to long format
  arrange(desc(Abundance)) #Organizes the table by Abundance in descending order

family_ITSY2 <- phyloseqITSY2_family_prop %>%  
  transform_sample_counts(function(x) {x *100} ) %>% #Recalculates the relative abundance 
  psmelt() %>%  #Melts to long format
  arrange(desc(Abundance)) #Organizes the table by Abundance in descending order

write.csv(family_16sY1, "family_16sY1.csv")
write.csv(family_16sY2, "family_16sY2.csv")
write.csv(family_ITSY1, "family_ITSY1.csv")
write.csv(family_ITSY2, "family_ITSY2.csv")

##Open in Excel and filter the 'Abundance' column to 'less than' 5
#All of the outcome rows are representative families with abundance lower than 0.05 
#Change all column labels under 'Family' to "Other (less than 5% relative abundance)" 
#Calculate relative abundance of "Others" per sample and collapse to one row
#Add numbers to unique taxonomic families.
#Save as .csv file
#Open file in R

family_16sY1_CoDa5<-read.csv("family_16sY1_5.csv", sep=",", header=T)
family_16sY2_CoDa5<-read.csv("family_16sY2_5.csv", sep=",", header=T)
family_ITSY1_CoDa5<-read.csv("family_ITSY1_5.csv", sep=",", header=T)
family_ITSY2_CoDa5<-read.csv("family_ITSY2_5.csv", sep=",", header=T)


#Plots with only datapoints above 5% relative abundance
#Stacked barplots
#16s - Y1
family_16sY1_barplot5<- ggplot(family_16sY1_CoDa5, aes(x = SampleOrder, y = Abundance , fill = L.Family))  + 
  facet_grid(Facility~ L..monocytogenes)+
  geom_bar(stat = "identity", color='black') + 
  theme(legend.text=element_text(size=9), legend.title= element_text(size=10, face = 'bold')) +
  theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.title.y=element_text(size=15),axis.text.y = element_text(color='black')) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  ylab("Relative Abundance (%)") + 
  theme(panel.background = element_rect(fill="transparent", color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_blank(), strip.text = element_text(size=15),
        panel.border = element_rect(color="black", fill=NA))+
  geom_text(aes(label=Label), position = position_stack(vjust=0.5), size=2)+
  ggtitle("Bacteria - Y1")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=15))+
  scale_fill_discrete(name="Family")
family_16sY1_barplot5

ggsave("Barplot_CoDa_16sY1.png", plot=family_16sY1_barplot5, device="png", width=13, height=8, units="in", dpi=600)


#For an interactive plot where you can hover over the axis and get the family
barplot_16sY1_animated<-ggplotly(family_16sY1_barplot5, tooltip=c('Family','Abundance'))
saveWidget(barplot_16sY1_animated, file="16sY1_animated.html")

#16s - Y2
family_16sY2_barplot5<-  ggplot(family_16sY2_CoDa5, aes(x = SampleOrder, y = Abundance , fill = L.Family))  + 
  facet_grid(Facility~ L..monocytogenes)+
  geom_bar(stat = "identity", color='black') + 
  theme(legend.text=element_text(size=9), legend.title= element_text(size=10, face = 'bold')) +
  theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.title.y=element_text(size=15),axis.text.y = element_text(color='black')) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  ylab("Relative Abundance (%)") + 
  theme(panel.background = element_rect(fill="transparent", color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_blank(), strip.text = element_text(size=15),
        panel.border = element_rect(color="black", fill=NA))+
  geom_text(aes(label=Label), position = position_stack(vjust=0.5), size=2)+
  ggtitle("Bacteria - Y2")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=15))+
  scale_fill_discrete(name="Family")
family_16sY2_barplot5

ggsave("Barplot_CoDa_16sY2.png", plot=family_16sY2_barplot5, device="png", width=13, height=8, units="in", dpi=600)

#Interactive plot
barplot_16sY2_animated<-ggplotly(family_16sY2_barplot5, tooltip=c('Family','Abundance'))
saveWidget(barplot_16sY2_animated, file="16sY2_animated.html")

#ITS - Y1
family_ITSY1_barplot5<- ggplot(family_ITSY1_CoDa5, aes(x = SampleOrder, y = Abundance , fill = L.Family))  + 
  facet_grid(Facility~ L..monocytogenes)+
  geom_bar(stat = "identity", color='black') + 
  theme(legend.text=element_text(size=9), legend.title= element_text(size=10, face = 'bold')) +
  theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.title.y=element_text(size=15),axis.text.y = element_text(color='black')) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  ylab("Relative Abundance (%)") + 
  theme(panel.background = element_rect(fill="transparent", color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_blank(), strip.text = element_text(size=15),
        panel.border = element_rect(color="black", fill=NA))+
  geom_text(aes(label=Label), position = position_stack(vjust=0.5), size=2)+
  ggtitle("Fungi - Y1")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=15))+
  scale_fill_discrete(name="Family")
family_ITSY1_barplot5

ggsave("Barplot_CoDa_ITSY1.png", plot=family_ITSY1_barplot5, device="png", width=13, height=8, units="in", dpi=600)

#Interactive plot
barplot_ITSY1_animated<-ggplotly(family_ITSY1_barplot5, tooltip=c('Family','Abundance'))
saveWidget(barplot_ITSY1_animated, file="ITSY1_animated.html")

#ITS - Y2
family_ITSY2_barplot5<- ggplot(family_ITSY2_CoDa5, aes(x = SampleOrder, y = Abundance , fill = L.Family))  + 
  facet_grid(Facility~ L..monocytogenes)+
  geom_bar(stat = "identity", color='black') + 
  theme(legend.text=element_text(size=9), legend.title= element_text(size=10, face = 'bold')) +
  theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.title.y=element_text(size=15),axis.text.y = element_text(color='black')) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  ylab("Relative Abundance (%)") + 
  theme(panel.background = element_rect(fill="transparent", color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_blank(), strip.text = element_text(size=15),
        panel.border = element_rect(color="black", fill=NA))+
  geom_text(aes(label=Label), position = position_stack(vjust=0.5), size=2)+
  ggtitle("Fungi - Y2")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=15))+
  scale_fill_discrete(name="Family")
family_ITSY2_barplot5

ggsave("Barplot_CoDa_ITSY2.png", plot=family_ITSY2_barplot5, device="png", width=13, height=8, units="in", dpi=600)

#Interactive plot
barplot_ITSY2_animated<-ggplotly(family_ITSY2_barplot5, tooltip=c('Family','Abundance'))
saveWidget(barplot_ITSY2_animated, file="ITSY2_animated.html")

### Boxplots - Plot with datapoints only above 5% relative abundance
#16sY1
boxplot_16sY1<-ggplot(family_16sY1_CoDa5, aes(x=Facility, y=Abundance))+
  geom_boxplot(color='grey')+ facet_wrap(~Family, scales = 'free_y', labeller = label_wrap_gen())+
  geom_jitter(aes(color=Facility))+
  scale_size(range=c(.1,10), name = "Relative abundance (%)")+
  theme(axis.title.x = element_blank(), axis.text.x=element_text(color='black', size=10), axis.ticks=element_line(color='black'),
        axis.title.y=element_text(size=15),axis.text.y = element_text(color='black')) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  ylab("Relative Abundance") + 
  theme(panel.background = element_rect(fill="grey99", color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_rect(fill="white", color = "black"), strip.text = element_text(size=7, face='bold'), 
        panel.border = element_rect(color="black", fill=NA))+
  ggtitle("Bacteria - Y1")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
  scale_color_viridis_d(begin = 0.2, end = 0.8, option='inferno')+
  theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))
boxplot_16sY1
ggsave("Boxplot_CoDa_16sY1.png", plot=boxplot_16sY1, device="png", width=16, height=8, units="in", dpi=600)


#16sY2
boxplot_16sY2<-ggplot(family_16sY2_CoDa5, aes(x=Facility, y=Abundance))+
  geom_boxplot(color='grey')+ facet_wrap(~Family, scales = 'free_y', labeller = label_wrap_gen())+
  geom_jitter(aes(color=Facility))+
  scale_size(range=c(.1,10), name = "Relative abundance (%)")+
  theme(axis.title.x = element_blank(), axis.text.x=element_text(color='black', size=10), axis.ticks=element_line(color='black'),
        axis.title.y=element_text(size=15),axis.text.y = element_text(color='black')) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  ylab("Relative Abundance (%)") + 
  theme(panel.background = element_rect(fill="grey99", color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_rect(fill="white", color = "black"), strip.text = element_text(size=7, face='bold'), 
        panel.border = element_rect(color="black", fill=NA))+
  ggtitle("Bacteria - Y2")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
  scale_color_viridis_d(begin = 0.2, end = 0.8, option='inferno')+
  theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))
boxplot_16sY2
ggsave("Boxplot_CoDa_16sY2.png", plot=boxplot_16sY2, device="png", width=16, height=8, units="in", dpi=600)

#ITSY1
boxplot_ITSY1<-ggplot(family_ITSY1_CoDa5, aes(x=Facility, y=Abundance))+
  geom_boxplot(color='grey')+ facet_wrap(~Family, scales = 'free_y', labeller = label_wrap_gen())+
  geom_jitter(aes(color=Facility))+
  scale_size(range=c(.1,10), name = "Relative abundance (%)")+
  theme(axis.title.x = element_blank(), axis.text.x=element_text(color='black', size=10), axis.ticks=element_line(color='black'),
        axis.title.y=element_text(size=15),axis.text.y = element_text(color='black')) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  ylab("Relative Abundance") + 
  theme(panel.background = element_rect(fill="grey99", color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_rect(fill="white", color = "black"), strip.text = element_text(size=7, face='bold'), 
        panel.border = element_rect(color="black", fill=NA))+
  ggtitle("Fungi - Y1")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
  scale_color_viridis_d(begin = 0.2, end = 0.8, option='inferno')+
  theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))
boxplot_ITSY1
ggsave("Boxplot_CoDa_ITSY1.png", plot=boxplot_ITSY1, device="png", width=16, height=8, units="in", dpi=600)



#ITSY2
boxplot_ITSY2<-ggplot(family_ITSY2_CoDa5, aes(x=Facility, y=Abundance))+
  geom_boxplot(color='grey')+ facet_wrap(~Family, scales = 'free_y', labeller = label_wrap_gen())+
  geom_jitter(aes(color=Facility))+
  scale_size(range=c(.1,10), name = "Relative abundance (%)")+
  theme(axis.title.x = element_blank(), axis.text.x=element_text(color='black', size=10), axis.ticks=element_line(color='black'),
        axis.title.y=element_text(size=15),axis.text.y = element_text(color='black')) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  ylab("Relative Abundance") + 
  theme(panel.background = element_rect(fill="grey99", color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_rect(fill="white", color = "black"), strip.text = element_text(size=7, face='bold'), 
        panel.border = element_rect(color="black", fill=NA))+
  ggtitle("Fungi - Y2")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
  scale_color_viridis_d(begin = 0.2, end = 0.8, option='inferno')+
  theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))
boxplot_ITSY2
ggsave("Boxplot_CoDa_ITSY2.png", plot=boxplot_ITSY2, device="png", width=16, height=8, units="in", dpi=600)

#Boxplots with all datapoints, for families that appeared over 5% relative abundance
#16sY1
family_16sY1_over5<-c('Aeromonadaceae','Alphaproteobacteria_unclassified','Arcobacteraceae',
                      'Azospirillaceae','Bacteria_unclassified',  'Bdellovibrionaceae','Beijerinckiaceae',
                      'Blastocatellaceae','Burkholderiaceae','Caulobacteraceae','Chitinophagaceae','Enterobacteriaceae',
                      'env.OPS_17','Flavobacteriaceae','Gammaproteobacteria_unclassified','Leuconostocaceae',
                      'Microbacteriaceae','Moraxellaceae','Mycobacteriaceae','Propionibacteriaceae','Pseudomonadaceae',
                      'Rhizobiaceae','Rhodanobacteraceae','Rhodobacteraceae','Rhodocyclaceae','Shewanellaceae',
                      'Sphingobacteriaceae','Sphingomonadaceae','Weeksellaceae','Xanthomonadaceae')



family_16sY1_over5abund <- filter(family_16sY1, Family %in% family_16sY1_over5)

boxplot_16sY1_all<-ggplot(family_16sY1_over5abund, aes(x=Facility, y=Abundance))+
  geom_boxplot(color='grey')+ facet_wrap(~Family, scales = 'free_y', labeller = label_wrap_gen())+
  geom_jitter(aes(color=Facility))+
  scale_size(range=c(.1,10), name = "Relative abundance (%)")+
  theme(axis.title.x = element_blank(), axis.text.x=element_text(color='black', size=10), axis.ticks=element_line(color='black'),
        axis.title.y=element_text(size=15),axis.text.y = element_text(color='black')) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  ylab("Relative Abundance") + 
  theme(panel.background = element_rect(fill="grey99", color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_rect(fill="white", color = "black"), strip.text = element_text(size=7, face='bold'), 
        panel.border = element_rect(color="black", fill=NA))+
  ggtitle("Bacteria - Y1")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
  scale_color_viridis_d(begin = 0.2, end = 0.8, option='inferno')+
  theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))
boxplot_16sY1_all
ggsave("Boxplot_16sY1_all.png", plot=boxplot_16sY1_all, device="png", width=16, height=8, units="in", dpi=600)

#16sY2
family_16sY2_over5<-c('Aeromonadaceae','Arcobacteraceae','Azospirillaceae','Bacteria_unclassified',
'Bacteriovoracaceae','Bacteroidaceae','Bdellovibrionaceae','Beijerinckiaceae','Burkholderiaceae',
'Caulobacteraceae','Chitinophagaceae','Enterobacteriaceae','Flavobacteriaceae','Gammaproteobacteria_unclassified',
'Intrasporangiaceae','JG30-KF-CM45','Leuconostocaceae','Microbacteriaceae','Micrococcaceae','Micrococcales_unclassified',
'Moraxellaceae','Mycobacteriaceae','Nocardiaceae','Nocardioidaceae','Propionibacteriaceae','Pseudomonadaceae',
'Rhizobiaceae','Rhodobacteraceae','Rhodocyclaceae','Rhodopirillaceae','Shewanellaceae','Sphingobacteriaceae',
'Sphingomonadaceae','Thiovulaceae','Weeksellaceae','Xanthomonadaceae')

family_16sY2_over5abund <- filter(family_16sY2, Family %in% family_16sY2_over5)

boxplot_16sY2_all<-ggplot(family_16sY2_over5abund, aes(x=Facility, y=Abundance))+
  geom_boxplot(color='grey')+ facet_wrap(~Family, scales = 'free_y', labeller = label_wrap_gen())+
  geom_jitter(aes(color=Facility))+
  scale_size(range=c(.1,10), name = "Relative abundance (%)")+
  theme(axis.title.x = element_blank(), axis.text.x=element_text(color='black', size=10), axis.ticks=element_line(color='black'),
        axis.title.y=element_text(size=15),axis.text.y = element_text(color='black')) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  ylab("Relative Abundance") + 
  theme(panel.background = element_rect(fill="grey99", color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_rect(fill="white", color = "black"), strip.text = element_text(size=7, face='bold'), 
        panel.border = element_rect(color="black", fill=NA))+
  ggtitle("Bacteria - Y2")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
  scale_color_viridis_d(begin = 0.2, end = 0.8, option='inferno')+
  theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))
boxplot_16sY2_all
ggsave("Boxplot_16sY2_all.png", plot=boxplot_16sY2_all, device="png", width=16, height=8, units="in", dpi=600)


#ITSY1
family_ITSY1_over5<-c('Ascomycota_unclassified','Aspergillaceae','Aureobasidiaceae','Botryosphaeriaceae',
'Bulleribasidiaceae','Cucurbitariaceae','Didymellaceae','Didymosphaeriaceae','Dipodascaceae','Helotiales_fam_Incertae_sedis',
'Mycosphaerellaceae','Nectriaceae','Phacidiaceae','Phaeosphaeriaceae','Pichiaceae','Pleosporaceae',
'Pleosporales_unclassified','Sclerotiniaceae','Sordariomycetes_unclassified','Trichomeriaceae',
'Trichosporonaceae','unclassified_Saccharomycetales')


family_ITSY1_over5abund <- filter(family_ITSY1, Family %in% family_ITSY1_over5)

boxplot_ITSY1_all<-ggplot(family_ITSY1_over5abund, aes(x=Facility, y=Abundance))+
  geom_boxplot(color='grey')+ facet_wrap(~Family, scales = 'free_y', labeller = label_wrap_gen())+
  geom_jitter(aes(color=Facility))+
  scale_size(range=c(.1,10), name = "Relative abundance (%)")+
  theme(axis.title.x = element_blank(), axis.text.x=element_text(color='black', size=10), axis.ticks=element_line(color='black'),
        axis.title.y=element_text(size=15),axis.text.y = element_text(color='black')) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  ylab("Relative Abundance") + 
  theme(panel.background = element_rect(fill="grey99", color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_rect(fill="white", color = "black"), strip.text = element_text(size=7, face='bold'), 
        panel.border = element_rect(color="black", fill=NA))+
  ggtitle("Fungi - Y1")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
  scale_color_viridis_d(begin = 0.2, end = 0.8, option='inferno')+
  theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))
boxplot_ITSY1_all
ggsave("Boxplot_ITSY1_all.png", plot=boxplot_ITSY1_all, device="png", width=16, height=8, units="in", dpi=600)

#ITSY2
family_ITSY2_over5<-c('Ascomycota_unclassified','Aspergillaceae','Aureobasidiaceae','Bulleribasidiaceae',
'Cucurbitariaceae','Cyphellophoraceae','Cystobasidiaceae','Didymellaceae','Didymosphaeriaceae',
'Dipodascaceae','Filobasidiaceae','Herpotrichiellaceae','Mrakiaceae','Mucoraceae','Mycosphaerellaceae',
'Nectriaceae','Peniophoraceae','Phacidiaceae','Phaeosphaeriaceae','Phaffomycetaceae','Physalacriaceae','Pichiaceae',
'Pleosporaceae','Pleosporales_fam_Incertae_sedis','Pleosporales_unclassified','Saccharomycetales_fam_Incertae_sedis',
'Saccharomycetales_unclassified','Sclerotiniaceae','Sporidiobolaceae','Tremellaceae','Trichomeriaceae',
'Trichosporonaceae','unclassified_Saccharomycetales','Wallemiaceae')
      
family_ITSY2_over5abund <- filter(family_ITSY2, Family %in% family_ITSY2_over5)

boxplot_ITSY2_all<-ggplot(family_ITSY2_over5abund, aes(x=Facility, y=Abundance))+
  geom_boxplot(color='grey')+ facet_wrap(~Family, scales = 'free_y', labeller = label_wrap_gen())+
  geom_jitter(aes(color=Facility))+
  scale_size(range=c(.1,10), name = "Relative abundance (%)")+
  theme(axis.title.x = element_blank(), axis.text.x=element_text(color='black', size=10), axis.ticks=element_line(color='black'),
        axis.title.y=element_text(size=15),axis.text.y = element_text(color='black')) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  ylab("Relative Abundance") + 
  theme(panel.background = element_rect(fill="grey99", color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_rect(fill="white", color = "black"), strip.text = element_text(size=7, face='bold'), 
        panel.border = element_rect(color="black", fill=NA))+
  ggtitle("Fungi - Y2")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
  scale_color_viridis_d(begin = 0.2, end = 0.8, option='inferno')+
  theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))
boxplot_ITSY2_all
ggsave("Boxplot_ITSY2_all.png", plot=boxplot_ITSY2_all, device="png", width=16, height=8, units="in", dpi=600)




### Bubbleplots by facility
#16sY1
bubbleplot_16sY1<-ggplot(family_16sY1_CoDa5, aes(x=Facility, y=reorder(Family,desc(Family)), size=Abundance, color=Facility))+
  geom_point()+
  scale_size(range=c(.1,6), name = "Relative abundance (%)")+
  theme(axis.title = element_blank(), axis.text.x=element_text(color='black', size=8), axis.ticks=element_line(color='black'),
        axis.text.y = element_text(color='black', size=6)) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  theme(panel.background = element_rect(fill="grey99", color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_rect(fill="white", color = "black"), strip.text = element_text(size=10, face='bold'), 
        panel.border = element_rect(color="black", fill=NA))+
  ggtitle("Bacteria - Y1")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
  scale_color_viridis_d(begin = 0.2, end = 0.8, option='inferno')+
  theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))
bubbleplot_16sY1


#16sY2
bubbleplot_16sY2<-ggplot(family_16sY2_CoDa5, aes(x=Facility, y=reorder(Family,desc(Family)), size=Abundance, color=Facility))+
  geom_point()+
  scale_size(range=c(.1,8), name = "Relative abundance (%)")+
  theme(axis.title = element_blank(), axis.text.x=element_text(color='black', size=8), axis.ticks=element_line(color='black'),
        axis.text.y = element_text(color='black', size=6)) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  theme(panel.background = element_rect(fill="grey99", color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_rect(fill="white", color = "black"), strip.text = element_text(size=10, face='bold'), 
        panel.border = element_rect(color="black", fill=NA))+
  ggtitle("Bacteria - Y2")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
  scale_color_viridis_d(begin = 0.2, end = 0.8, option='inferno')+
  theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))
bubbleplot_16sY2



#ITSY1
bubbleplot_ITSY1<-ggplot(family_ITSY1_CoDa5, aes(x=Facility, y=reorder(Family,desc(Family)), size=Abundance, color=Facility))+
  geom_point()+
  scale_size(range=c(.1,8), name = "Relative abundance (%)")+
  theme(axis.title= element_blank(), axis.text.x=element_text(color='black', size=8), axis.ticks=element_line(color='black'),
        axis.text.y = element_text(color='black', size=6)) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  theme(panel.background = element_rect(fill="grey99", color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_rect(fill="white", color = "black"), strip.text = element_text(size=10, face='bold'), 
        panel.border = element_rect(color="black", fill=NA))+
  ggtitle("Fungi - Y1")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
  scale_color_viridis_d(begin = 0.2, end = 0.8, option='inferno')+
  theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))
bubbleplot_ITSY1


#ITSY2
bubbleplot_ITSY2<-ggplot(family_ITSY2_CoDa5, aes(x=Facility, y=reorder(Family,desc(Family)), size=Abundance, color=Facility))+
  geom_point()+
  scale_size(range=c(.1,8), name = "Relative abundance (%)")+
  theme(axis.title = element_blank(), axis.text.x=element_text(color='black', size=8), axis.ticks=element_line(color='black'),
        axis.text.y = element_text(color='black', size=6)) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  theme(panel.background = element_rect(fill="grey99", color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_rect(fill="white", color = "black"), strip.text = element_text(size=10, face='bold'), 
        panel.border = element_rect(color="black", fill=NA))+
  ggtitle("Fungi - Y2")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
  scale_color_viridis_d(begin = 0.2, end = 0.8, option='inferno')+
  theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))+
  theme(legend.position = "bottom")
bubbleplot_ITSY2


#Combine bubble plots
Bubbleplots_all = plot_grid(bubbleplot_16sY1 +theme(legend.position = "none"),
                            bubbleplot_16sY2+theme(legend.position = "none"), 
                            bubbleplot_ITSY1+theme(legend.position = "none"), 
                            bubbleplot_ITSY2 +theme(legend.position = "none"),
                       ncol=2, nrow=2, labels = c("A","B","C","D"), label_size = 20, vjust = 2, hjust = -1.5)
Bubbleplots_all

Bubbleplots_all_leg<-get_legend(bubbleplot_ITSY2)

Bubbleplot_final<-plot_grid(Bubbleplots_all, Bubbleplots_all_leg,
                            ncol=1, nrow=2, rel_heights = c(7,1))
Bubbleplot_final


ggsave("Bubbleplot_family_facility.png", plot=Bubbleplot_final, device="png", width=8, height=11, units="in", dpi=600)

### Bubbleplots by facility and Lm
#16sY1
bubbleplot_16sY1_Lm<-ggplot(family_16sY1_CoDa5, aes(x=Facility, y=reorder(Family,desc(Family)), size=Abundance, color=Facility))+
  geom_point()+facet_grid(~L..monocytogenes)+
  scale_size(range=c(.1,6), name = "Relative abundance (%)")+
  theme(axis.title = element_blank(), axis.text.x=element_text(color='black', size=8), axis.ticks=element_line(color='black'),
        axis.text.y = element_text(color='black', size=6)) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  theme(panel.background = element_rect(fill="grey99", color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_rect(fill="white", color = "black"), strip.text = element_text(size=10, face='bold'), 
        panel.border = element_rect(color="black", fill=NA))+
  ggtitle("Bacteria - Y1")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
  scale_color_viridis_d(begin = 0.2, end = 0.8, option='inferno')+
  theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))+
  theme(legend.position = "bottom")
bubbleplot_16sY1_Lm


#16sY2
bubbleplot_16sY2_Lm<-ggplot(family_16sY2_CoDa5, aes(x=Facility, y=reorder(Family,desc(Family)), size=Abundance, color=Facility))+
  geom_point()+facet_grid(~L..monocytogenes)+
  scale_size(range=c(.1,6), name = "Relative abundance (%)")+
  theme(axis.title = element_blank(), axis.text.x=element_text(color='black', size=8), axis.ticks=element_line(color='black'),
        axis.text.y = element_text(color='black', size=6)) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  theme(panel.background = element_rect(fill="grey99", color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_rect(fill="white", color = "black"), strip.text = element_text(size=10, face='bold'), 
        panel.border = element_rect(color="black", fill=NA))+
  ggtitle("Bacteria - Y2")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
  scale_color_viridis_d(begin = 0.2, end = 0.8, option='inferno')+
  theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))+
  theme(legend.position = "bottom")
bubbleplot_16sY2_Lm

#ITSY1
bubbleplot_ITSY1_Lm<-ggplot(family_ITSY1_CoDa5, aes(x=Facility, y=reorder(Family,desc(Family)), size=Abundance, color=Facility))+
  geom_point()+facet_grid(~L..monocytogenes)+
  scale_size(range=c(.1,6), name = "Relative abundance (%)")+
  theme(axis.title = element_blank(), axis.text.x=element_text(color='black', size=8), axis.ticks=element_line(color='black'),
        axis.text.y = element_text(color='black', size=6)) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  theme(panel.background = element_rect(fill="grey99", color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_rect(fill="white", color = "black"), strip.text = element_text(size=10, face='bold'), 
        panel.border = element_rect(color="black", fill=NA))+
  ggtitle("Bacteria - Y1")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
  scale_color_viridis_d(begin = 0.2, end = 0.8, option='inferno')+
  theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))+ 
  theme(legend.position = "bottom")
bubbleplot_ITSY1_Lm

#ITSY2
bubbleplot_ITSY2_Lm<-ggplot(family_ITSY2_CoDa5, aes(x=Facility, y=reorder(Family,desc(Family)), size=Abundance, color=Facility))+
  geom_point()+facet_grid(~L..monocytogenes)+
  scale_size(range=c(.1,6), name = "Relative abundance (%)")+
  theme(axis.title = element_blank(), axis.text.x=element_text(color='black', size=8), axis.ticks=element_line(color='black'),
        axis.text.y = element_text(color='black', size=6)) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  theme(panel.background = element_rect(fill="grey99", color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_rect(fill="white", color = "black"), strip.text = element_text(size=10, face='bold'), 
        panel.border = element_rect(color="black", fill=NA))+
  ggtitle("Bacteria - Y2")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
  scale_color_viridis_d(begin = 0.2, end = 0.8, option='inferno')+
  theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))+
  theme(legend.position = "bottom")
bubbleplot_ITSY2_Lm


#Combine bubble plots
Bubbleplots_all_Lm = plot_grid(bubbleplot_16sY1_Lm +theme(legend.position = "none"),
                            bubbleplot_16sY2_Lm+theme(legend.position = "none"), 
                            bubbleplot_ITSY1_Lm+theme(legend.position = "none"), 
                            bubbleplot_ITSY2_Lm +theme(legend.position = "none"),
                            ncol=2, nrow=2, labels = c("A","B","C","D"), label_size = 20, vjust = 2, hjust = -1.5)
Bubbleplots_all_Lm

Bubbleplots_all_leg_Lm<-get_legend(bubbleplot_ITSY2_Lm)

Bubbleplot_final_Lm<-plot_grid(Bubbleplots_all_Lm, Bubbleplots_all_leg_Lm,
                            ncol=1, nrow=2, rel_heights = c(7,1))
Bubbleplot_final_Lm


ggsave("Bubbleplot_family_facility_Lm.png", plot=Bubbleplot_final_Lm, device="png", width=8, height=11, units="in", dpi=600)


#### CORE MICROBIOME ####
library(microbiome)
library(knitr)


# Subset data for each facility
F1_16sY1 <- subset_samples(phyloseq16sY1_family_prop, Facility == "F1") 
F2_16sY1 <- subset_samples(phyloseq16sY1_family_prop, Facility == "F2") 
F3_16sY1 <- subset_samples(phyloseq16sY1_family_prop, Facility == "F3") 

F1_16sY2 <- subset_samples(phyloseq16sY2_family_prop, Facility == "F1") 
F2_16sY2 <- subset_samples(phyloseq16sY2_family_prop, Facility == "F2") 
F3_16sY2 <- subset_samples(phyloseq16sY2_family_prop, Facility == "F3") 

F1_ITSY1 <- subset_samples(phyloseqITSY1_family_prop, Facility == "F1") 
F2_ITSY1 <- subset_samples(phyloseqITSY1_family_prop, Facility == "F2") 
F3_ITSY1 <- subset_samples(phyloseqITSY1_family_prop, Facility == "F3") 

F1_ITSY2 <- subset_samples(phyloseqITSY2_family_prop, Facility == "F1") 
F2_ITSY2 <- subset_samples(phyloseqITSY2_family_prop, Facility == "F2") 
F3_ITSY2 <- subset_samples(phyloseqITSY2_family_prop, Facility == "F3") 


#Core microbiota by facility
#Families detected in at least 75% of the samples with a relative abundance threshold value above 0.01%
F1_16sY1_core <- core(F1_16sY1, detection = 0.001, prevalence = .75) 
F2_16sY1_core <- core(F2_16sY1, detection = 0.001, prevalence = .75)
F3_16sY1_core <- core(F3_16sY1, detection = 0.001, prevalence = .75)

F1_16sY2_core <- core(F1_16sY2, detection = 0.001, prevalence = .75) 
F2_16sY2_core <- core(F2_16sY2, detection = 0.001, prevalence = .75)
F3_16sY2_core <- core(F3_16sY2, detection = 0.001, prevalence = .75)

F1_ITSY1_core <- core(F1_ITSY1, detection = 0.001, prevalence = .75) 
F2_ITSY1_core <- core(F2_ITSY1, detection = 0.001, prevalence = .75)
F3_ITSY1_core <- core(F3_ITSY1, detection = 0.001, prevalence = .75)

F1_ITSY2_core <- core(F1_ITSY2, detection = 0.001, prevalence = .75) 
F2_ITSY2_core <- core(F2_ITSY2, detection = 0.001, prevalence = .75)
F3_ITSY2_core <- core(F3_ITSY2, detection = 0.001, prevalence = .75)


# get the taxonomy data
F1_16sY1_tax <- as.data.frame(tax_table(F1_16sY1_core)) #17 taxa
F2_16sY1_tax <- as.data.frame(tax_table(F2_16sY1_core)) #14 taxa
F3_16sY1_tax <- as.data.frame(tax_table(F3_16sY1_core)) #18 taxa

F1_16sY2_tax <- as.data.frame(tax_table(F1_16sY2_core)) #23 taxa
F2_16sY2_tax <- as.data.frame(tax_table(F2_16sY2_core)) #29 taxa
F3_16sY2_tax <- as.data.frame(tax_table(F3_16sY2_core)) #32 taxa

F1_ITSY1_tax <- as.data.frame(tax_table(F1_ITSY1_core)) #14 taxa
F2_ITSY1_tax <- as.data.frame(tax_table(F2_ITSY1_core)) #12 taxa
F3_ITSY1_tax <- as.data.frame(tax_table(F3_ITSY1_core)) #13 taxa

F1_ITSY2_tax <- as.data.frame(tax_table(F1_ITSY2_core)) #21 taxa
F2_ITSY2_tax <- as.data.frame(tax_table(F2_ITSY2_core)) #20 taxa
F3_ITSY2_tax <- as.data.frame(tax_table(F3_ITSY2_core)) #19 taxa

#Extract family
F1_16sY1_family<-F1_16sY1_tax$Family
F2_16sY1_family<-F2_16sY1_tax$Family
F3_16sY1_family<-F3_16sY1_tax$Family

F1_16sY2_family<-F1_16sY2_tax$Family
F2_16sY2_family<-F2_16sY2_tax$Family
F3_16sY2_family<-F3_16sY2_tax$Family

F1_ITSY1_family<-F1_ITSY1_tax$Family
F2_ITSY1_family<-F2_ITSY1_tax$Family
F3_ITSY1_family<-F3_ITSY1_tax$Family

F1_ITSY2_family<-F1_ITSY2_tax$Family
F2_ITSY2_family<-F2_ITSY2_tax$Family
F3_ITSY2_family<-F3_ITSY2_tax$Family

#Make list with all core taxa for each category
coretaxa_16sY1_list<-list(F1=F1_16sY1_family, F2=F2_16sY1_family, F3=F3_16sY1_family)
coretaxa_16sY2_list<-list(F1=F1_16sY2_family, F2=F2_16sY2_family, F3=F3_16sY2_family)
coretaxa_ITSY1_list<-list(F1=F1_ITSY1_family, F2=F2_ITSY1_family, F3=F3_ITSY1_family)
coretaxa_ITSY2_list<-list(F1=F1_ITSY2_family, F2=F2_ITSY2_family, F3=F3_ITSY2_family)

library(ggvenn)

venn_16sY1<-ggvenn(coretaxa_16sY1_list, show_percentage = FALSE,
  fill_color = c('#420A68','#BB3754','#FCA50A'), set_name_size = 4, stroke_size = 0.5)

venn_16sY2<-ggvenn(coretaxa_16sY2_list, show_percentage = FALSE,
  fill_color = c('#420A68','#BB3754','#FCA50A'), set_name_size = 4, stroke_size = 0.5)

venn_ITSY1<-ggvenn(coretaxa_ITSY1_list, show_percentage = FALSE,
  fill_color = c('#420A68','#BB3754','#FCA50A'), set_name_size = 4, stroke_size = 0.5)

venn_ITSY2<-ggvenn(coretaxa_ITSY2_list, show_percentage = FALSE,
  fill_color = c('#420A68','#BB3754','#FCA50A'), set_name_size = 4, stroke_size = 0.5)


#Combine Venn diagrams
Venn_all_Fac = plot_grid(venn_16sY1, venn_16sY2, venn_ITSY1, venn_ITSY2,
                        ncol=2, nrow=2, labels = c("A","B","C","D"), label_size = 20, vjust = 2, hjust = -1.5)
Venn_all_Fac

ggsave("Venn_core_facility.png", plot=Venn_all_Fac, device="png", width=8, height=8, units="in", dpi=600)


#Compare core taxa by year  for each facility
coretaxa_16sF1_list<-list(Y1=F1_16sY1_family, Y2=F1_16sY2_family)
coretaxa_16sF2_list<-list(Y1=F2_16sY1_family, Y2=F2_16sY2_family)
coretaxa_16sF3_list<-list(Y1=F3_16sY1_family, Y2=F3_16sY2_family)

coretaxa_ITSF1_list<-list(Y1=F1_ITSY1_family, Y2=F1_ITSY2_family)
coretaxa_ITSF2_list<-list(Y1=F2_ITSY1_family, Y2=F2_ITSY2_family)
coretaxa_ITSF3_list<-list(Y1=F3_ITSY1_family, Y2=F3_ITSY2_family)

#Make Venn diagrams

venn_16sF1<-ggvenn(coretaxa_16sF1_list, show_percentage = FALSE,
                   fill_color = c('#420A68','#8D6CA4'), set_name_size = 4, stroke_size = 0.5)

venn_16sF2<-ggvenn(coretaxa_16sF2_list, show_percentage = FALSE,
                   fill_color = c('#BB3754','#D68798'), set_name_size = 4, stroke_size = 0.5)

venn_16sF3<-ggvenn(coretaxa_16sF3_list, show_percentage = FALSE,
                   fill_color = c('#FCA50A','#FDC96C'), set_name_size = 4, stroke_size = 0.5)


venn_ITSF1<-ggvenn(coretaxa_ITSF1_list, show_percentage = FALSE,
                   fill_color = c('#420A68','#8D6CA4'), set_name_size = 4, stroke_size = 0.5)

venn_ITSF2<-ggvenn(coretaxa_ITSF2_list, show_percentage = FALSE,
                   fill_color = c('#BB3754','#D68798'), set_name_size = 4, stroke_size = 0.5)

venn_ITSF3<-ggvenn(coretaxa_ITSF3_list, show_percentage = FALSE,
                   fill_color = c('#FCA50A','#FDC96C'), set_name_size = 4, stroke_size = 0.5)


#Combine Venn diagrams by year
Venn_all_Year = plot_grid(venn_16sF1, venn_16sF2, venn_16sF3, venn_ITSF1, venn_ITSF2, venn_ITSF3,
                         ncol=3, nrow=2, labels = c("A","B","C","D","E","F"), label_size = 20, vjust = 2, hjust = -1.5)
Venn_all_Year

ggsave("Venn_core_year.png", plot=Venn_all_Year, device="png", width=8, height=4, units="in", dpi=600)


#Make barplots with abundance of core microbiome by facility
#Open in Excel the files family_16sY1.csv,  family_16sY2.csv, family_ITSY1.csv, family_ITSY2.csv for the core bactera shared by all facilities in 2 years.
#Copy and paste in a new sheet and save as .csv

family_16sY1_core<-read.csv('core_16sY1.csv', header = TRUE)
family_16sY2_core<-read.csv('core_16sY2.csv', header = TRUE)
family_ITSY1_core<-read.csv('core_ITSY1.csv', header = TRUE)
family_ITSY2_core<-read.csv('core_ITSY2.csv', header = TRUE)

#Add year for all observations
family_16sY1_core$Year<-rep("Y1", 1170)
family_16sY2_core$Year<-rep("Y2", 1070)
family_ITSY1_core$Year<-rep("Y1", 936)
family_ITSY2_core$Year<-rep("Y2", 856)

#Combine the two years data
family_16s_core<-bind_rows(family_16sY1_core, family_16sY2_core)
family_ITS_core<-bind_rows(family_ITSY1_core, family_ITSY2_core)

#Plot
#16s
boxplot_16s_core<-ggplot(family_16s_core, aes(x=Facility, y=Abundance))+
  geom_boxplot(color='grey')+ facet_wrap(~Family+Year, scales = 'free_y', labeller = label_wrap_gen(), nrow=5, ncol=4)+
  geom_jitter(aes(color=Facility))+
  scale_size(range=c(.1,10), name = "Relative abundance (%)")+
  theme(axis.title.x = element_blank(), axis.text.x=element_text(color='black', size=10), axis.ticks=element_line(color='black'),
        axis.title.y=element_text(size=15),axis.text.y = element_text(color='black')) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  ylab("Relative Abundance") + 
  theme(panel.background = element_rect(fill="grey99", color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_rect(fill="white", color = "black"), strip.text = element_text(size=7, face='bold'), 
        panel.border = element_rect(color="black", fill=NA))+
  ggtitle("Core Bacteria")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
  scale_color_viridis_d(begin = 0.2, end = 0.8, option='inferno')+
  theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))
boxplot_16s_core

ggsave("Boxplot_core_16s.png", plot=boxplot_16s_core, device="png", width=10, height=10, units="in", dpi=600)

#ITS
boxplot_ITS_core<-ggplot(family_ITS_core, aes(x=Facility, y=Abundance))+
  geom_boxplot(color='grey')+ facet_wrap(~Family+Year, scales = 'free_y', labeller = label_wrap_gen(), nrow=4, ncol=4)+
  geom_jitter(aes(color=Facility))+
  scale_size(range=c(.1,10), name = "Relative abundance (%)")+
  theme(axis.title.x = element_blank(), axis.text.x=element_text(color='black', size=10), axis.ticks=element_line(color='black'),
        axis.title.y=element_text(size=15),axis.text.y = element_text(color='black')) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  ylab("Relative Abundance") + 
  theme(panel.background = element_rect(fill="grey99", color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_rect(fill="white", color = "black"), strip.text = element_text(size=7, face='bold'), 
        panel.border = element_rect(color="black", fill=NA))+
  ggtitle("Core Fungi")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
  scale_color_viridis_d(begin = 0.2, end = 0.8, option='inferno')+
  theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))
boxplot_ITS_core

ggsave("Boxplot_core_ITS.png", plot=boxplot_ITS_core, device="png", width=10, height=10, units="in", dpi=600)

#### DIFFERENTIAL ABUNDANCE BY FACILITY ####


# Subset Phyloseq with read counts at family level for each facility
F1_16sY1_count <- subset_samples(phyloseq16sY1_family, Facility == "F1") 
F2_16sY1_count <- subset_samples(phyloseq16sY1_family, Facility == "F2") 
F3_16sY1_count <- subset_samples(phyloseq16sY1_family, Facility == "F3") 

F1_16sY2_count <- subset_samples(phyloseq16sY2_family, Facility == "F1") 
F2_16sY2_count <- subset_samples(phyloseq16sY2_family, Facility == "F2") 
F3_16sY2_count <- subset_samples(phyloseq16sY2_family, Facility == "F3") 

F1_ITSY1_count <- subset_samples(phyloseqITSY1_family, Facility == "F1") 
F2_ITSY1_count <- subset_samples(phyloseqITSY1_family, Facility == "F2") 
F3_ITSY1_count <- subset_samples(phyloseqITSY1_family, Facility == "F3") 

F1_ITSY2_count <- subset_samples(phyloseqITSY2_family, Facility == "F1") 
F2_ITSY2_count <- subset_samples(phyloseqITSY2_family, Facility == "F2") 
F3_ITSY2_count <- subset_samples(phyloseqITSY2_family, Facility == "F3") 

#Take otu table for each facility from Phyloseq
family_F116sY1_otu<-as.data.frame(as(otu_table(F1_16sY1_count), "matrix"))
family_F216sY1_otu<-as.data.frame(as(otu_table(F2_16sY1_count), "matrix"))
family_F316sY1_otu<-as.data.frame(as(otu_table(F3_16sY1_count), "matrix"))

family_F116sY2_otu<-as.data.frame(as(otu_table(F1_16sY2_count), "matrix"))
family_F216sY2_otu<-as.data.frame(as(otu_table(F2_16sY2_count), "matrix"))
family_F316sY2_otu<-as.data.frame(as(otu_table(F3_16sY2_count), "matrix"))

family_F1ITSY1_otu<-as.data.frame(as(otu_table(F1_ITSY1_count), "matrix"))
family_F2ITSY1_otu<-as.data.frame(as(otu_table(F2_ITSY1_count), "matrix"))
family_F3ITSY1_otu<-as.data.frame(as(otu_table(F3_ITSY1_count), "matrix"))

family_F1ITSY2_otu<-as.data.frame(as(otu_table(F1_ITSY2_count), "matrix"))
family_F2ITSY2_otu<-as.data.frame(as(otu_table(F2_ITSY2_count), "matrix"))
family_F3ITSY2_otu<-as.data.frame(as(otu_table(F3_ITSY2_count), "matrix"))

#Merge OTU tables for each pair of facilities
family_otu_16sY1_F1F2<-as.data.frame(bind_cols(family_F116sY1_otu, family_F216sY1_otu)) #Otu table for aldex2 needs to have otus in rows
family_otu_16sY1_F1F3<-as.data.frame(bind_cols(family_F116sY1_otu, family_F316sY1_otu)) #Otu table for aldex2 needs to have otus in rows
family_otu_16sY1_F2F3<-as.data.frame(bind_cols(family_F216sY1_otu, family_F316sY1_otu)) #Otu table for aldex2 needs to have otus in rows

family_otu_16sY2_F1F2<-as.data.frame(bind_cols(family_F116sY2_otu, family_F216sY2_otu)) #Otu table for aldex2 needs to have otus in rows
family_otu_16sY2_F1F3<-as.data.frame(bind_cols(family_F116sY2_otu, family_F316sY2_otu)) #Otu table for aldex2 needs to have otus in rows
family_otu_16sY2_F2F3<-as.data.frame(bind_cols(family_F216sY2_otu, family_F316sY2_otu)) #Otu table for aldex2 needs to have otus in rows

family_otu_ITSY1_F1F2<-as.data.frame(bind_cols(family_F1ITSY1_otu, family_F2ITSY1_otu)) #Otu table for aldex2 needs to have otus in rows
family_otu_ITSY1_F1F3<-as.data.frame(bind_cols(family_F1ITSY1_otu, family_F3ITSY1_otu)) #Otu table for aldex2 needs to have otus in rows
family_otu_ITSY1_F2F3<-as.data.frame(bind_cols(family_F2ITSY1_otu, family_F3ITSY1_otu)) #Otu table for aldex2 needs to have otus in rows

family_otu_ITSY2_F1F2<-as.data.frame(bind_cols(family_F1ITSY2_otu, family_F2ITSY2_otu)) #Otu table for aldex2 needs to have otus in rows
family_otu_ITSY2_F1F3<-as.data.frame(bind_cols(family_F1ITSY2_otu, family_F3ITSY2_otu)) #Otu table for aldex2 needs to have otus in rows
family_otu_ITSY2_F2F3<-as.data.frame(bind_cols(family_F2ITSY2_otu, family_F3ITSY2_otu)) #Otu table for aldex2 needs to have otus in rows

#Subset metadata by pairs of facilities. Aldex needs a category with the same amount of samples as input
metadata_F116sY1<-metadata_16sY1[ which(metadata_16sY1$Facility =="F1"),]
metadata_F216sY1<-metadata_16sY1[ which(metadata_16sY1$Facility =="F2"),]
metadata_F316sY1<-metadata_16sY1[ which(metadata_16sY1$Facility =="F3"),]

metadata_F116sY2<-metadata_16sY2[ which(metadata_16sY2$Facility =="F1"),]
metadata_F216sY2<-metadata_16sY2[ which(metadata_16sY2$Facility =="F2"),]
metadata_F316sY2<-metadata_16sY2[ which(metadata_16sY2$Facility =="F3"),]

metadata_F1ITSY1<-metadata_ITSY1[ which(metadata_ITSY1$Facility =="F1"),]
metadata_F2ITSY1<-metadata_ITSY1[ which(metadata_ITSY1$Facility =="F2"),]
metadata_F3ITSY1<-metadata_ITSY1[ which(metadata_ITSY1$Facility =="F3"),]

metadata_F1ITSY2<-metadata_ITSY2[ which(metadata_ITSY2$Facility =="F1"),]
metadata_F2ITSY2<-metadata_ITSY2[ which(metadata_ITSY2$Facility =="F2"),]
metadata_F3ITSY2<-metadata_ITSY2[ which(metadata_ITSY2$Facility =="F3"),]

#Change Facility to character vector - as.factor affects aldex.effect function
metadata_F116sY1$Facility <- as.character(metadata_F116sY1$Facility)
metadata_F216sY1$Facility <- as.character(metadata_F216sY1$Facility)
metadata_F316sY1$Facility <- as.character(metadata_F316sY1$Facility)

metadata_F116sY2$Facility <- as.character(metadata_F116sY2$Facility)
metadata_F216sY2$Facility <- as.character(metadata_F216sY2$Facility)
metadata_F316sY2$Facility <- as.character(metadata_F316sY2$Facility)

metadata_F1ITSY1$Facility <- as.character(metadata_F1ITSY1$Facility)
metadata_F2ITSY1$Facility <- as.character(metadata_F2ITSY1$Facility)
metadata_F3ITSY1$Facility <- as.character(metadata_F3ITSY1$Facility)

metadata_F1ITSY2$Facility <- as.character(metadata_F1ITSY2$Facility)
metadata_F2ITSY2$Facility <- as.character(metadata_F2ITSY2$Facility)
metadata_F3ITSY2$Facility <- as.character(metadata_F3ITSY2$Facility)

#Bind metadata by pair of facilities
family_metadata_16sY1_F1F2<-bind_rows(metadata_F116sY1, metadata_F216sY1) 
family_metadata_16sY1_F1F3<-bind_rows(metadata_F116sY1, metadata_F316sY1)
family_metadata_16sY1_F2F3<-bind_rows(metadata_F216sY1, metadata_F316sY1)

family_metadata_16sY2_F1F2<-bind_rows(metadata_F116sY2, metadata_F216sY2) 
family_metadata_16sY2_F1F3<-bind_rows(metadata_F116sY2, metadata_F316sY2)
family_metadata_16sY2_F2F3<-bind_rows(metadata_F216sY2, metadata_F316sY2)

family_metadata_ITSY1_F1F2<-bind_rows(metadata_F1ITSY1, metadata_F2ITSY1) 
family_metadata_ITSY1_F1F3<-bind_rows(metadata_F1ITSY1, metadata_F3ITSY1)
family_metadata_ITSY1_F2F3<-bind_rows(metadata_F2ITSY1, metadata_F3ITSY1)

family_metadata_ITSY2_F1F2<-bind_rows(metadata_F1ITSY2, metadata_F2ITSY2) 
family_metadata_ITSY2_F1F3<-bind_rows(metadata_F1ITSY2, metadata_F3ITSY2)
family_metadata_ITSY2_F2F3<-bind_rows(metadata_F2ITSY2, metadata_F3ITSY2)

#Differential abundance analysis 
#Based on Gloor et al 2016. Annals of Epidemiology 26:322-329.
#Supplemental material page 19.

#Generate 128 Dirichlet distributed Monte Carlo instances and center-log ratio transform them
#OTU table needs to be as OTU in rows. Only two conditions can be compared at a time.

#Compare between samples that were positive and negative
#At family level
Aldex_family_16sY1_F1F2.clr<-aldex.clr(family_otu_16sY1_F1F2, mc.samples = 128, conds = family_metadata_16sY1_F1F2$Facility)
Aldex_family_16sY1_F1F3.clr<-aldex.clr(family_otu_16sY1_F1F3, mc.samples = 128, conds = family_metadata_16sY1_F1F3$Facility)
Aldex_family_16sY1_F2F3.clr<-aldex.clr(family_otu_16sY1_F2F3, mc.samples = 128, conds = family_metadata_16sY1_F2F3$Facility)

Aldex_family_16sY2_F1F2.clr<-aldex.clr(family_otu_16sY2_F1F2, mc.samples = 128, conds = family_metadata_16sY2_F1F2$Facility)
Aldex_family_16sY2_F1F3.clr<-aldex.clr(family_otu_16sY2_F1F3, mc.samples = 128, conds = family_metadata_16sY2_F1F3$Facility)
Aldex_family_16sY2_F2F3.clr<-aldex.clr(family_otu_16sY2_F2F3, mc.samples = 128, conds = family_metadata_16sY2_F2F3$Facility)

Aldex_family_ITSY1_F1F2.clr<-aldex.clr(family_otu_ITSY1_F1F2, mc.samples = 128, conds = family_metadata_ITSY1_F1F2$Facility)
Aldex_family_ITSY1_F1F3.clr<-aldex.clr(family_otu_ITSY1_F1F3, mc.samples = 128, conds = family_metadata_ITSY1_F1F3$Facility)
Aldex_family_ITSY1_F2F3.clr<-aldex.clr(family_otu_ITSY1_F2F3, mc.samples = 128, conds = family_metadata_ITSY1_F2F3$Facility)

Aldex_family_ITSY2_F1F2.clr<-aldex.clr(family_otu_ITSY2_F1F2, mc.samples = 128, conds = family_metadata_ITSY2_F1F2$Facility)
Aldex_family_ITSY2_F1F3.clr<-aldex.clr(family_otu_ITSY2_F1F3, mc.samples = 128, conds = family_metadata_ITSY2_F1F3$Facility)
Aldex_family_ITSY2_F2F3.clr<-aldex.clr(family_otu_ITSY2_F2F3, mc.samples = 128, conds = family_metadata_ITSY2_F2F3$Facility)


#Calculate the expected effect size
#At family level
Aldex_family_16sY1_F1F2.e<-aldex.effect(Aldex_family_16sY1_F1F2.clr)
Aldex_family_16sY1_F1F3.e<-aldex.effect(Aldex_family_16sY1_F1F3.clr)
Aldex_family_16sY1_F2F3.e<-aldex.effect(Aldex_family_16sY1_F2F3.clr)

Aldex_family_16sY2_F1F2.e<-aldex.effect(Aldex_family_16sY2_F1F2.clr)
Aldex_family_16sY2_F1F3.e<-aldex.effect(Aldex_family_16sY2_F1F3.clr)
Aldex_family_16sY2_F2F3.e<-aldex.effect(Aldex_family_16sY2_F2F3.clr)

Aldex_family_ITSY1_F1F2.e<-aldex.effect(Aldex_family_ITSY1_F1F2.clr)
Aldex_family_ITSY1_F1F3.e<-aldex.effect(Aldex_family_ITSY1_F1F3.clr)
Aldex_family_ITSY1_F2F3.e<-aldex.effect(Aldex_family_ITSY1_F2F3.clr)

Aldex_family_ITSY2_F1F2.e<-aldex.effect(Aldex_family_ITSY2_F1F2.clr)
Aldex_family_ITSY2_F1F3.e<-aldex.effect(Aldex_family_ITSY2_F1F3.clr)
Aldex_family_ITSY2_F2F3.e<-aldex.effect(Aldex_family_ITSY2_F2F3.clr)

#Generate the P and Benjamini-Hochberg corrected P values.
#At family level
Aldex_family_16sY1_F1F2.t<-aldex.ttest(Aldex_family_16sY1_F1F2.clr)
Aldex_family_16sY1_F1F3.t<-aldex.ttest(Aldex_family_16sY1_F1F3.clr)
Aldex_family_16sY1_F2F3.t<-aldex.ttest(Aldex_family_16sY1_F2F3.clr)

Aldex_family_16sY2_F1F2.t<-aldex.ttest(Aldex_family_16sY2_F1F2.clr)
Aldex_family_16sY2_F1F3.t<-aldex.ttest(Aldex_family_16sY2_F1F3.clr)
Aldex_family_16sY2_F2F3.t<-aldex.ttest(Aldex_family_16sY2_F2F3.clr)

Aldex_family_ITSY1_F1F2.t<-aldex.ttest(Aldex_family_ITSY1_F1F2.clr)
Aldex_family_ITSY1_F1F3.t<-aldex.ttest(Aldex_family_ITSY1_F1F3.clr)
Aldex_family_ITSY1_F2F3.t<-aldex.ttest(Aldex_family_ITSY1_F2F3.clr)

Aldex_family_ITSY2_F1F2.t<-aldex.ttest(Aldex_family_ITSY2_F1F2.clr)
Aldex_family_ITSY2_F1F3.t<-aldex.ttest(Aldex_family_ITSY2_F1F3.clr)
Aldex_family_ITSY2_F2F3.t<-aldex.ttest(Aldex_family_ITSY2_F2F3.clr)

#Merge data frames
#At family level
Aldex_family_16sY1_F1F2.all<-data.frame(Aldex_family_16sY1_F1F2.e,Aldex_family_16sY1_F1F2.t)
Aldex_family_16sY1_F1F3.all<-data.frame(Aldex_family_16sY1_F1F3.e,Aldex_family_16sY1_F1F3.t)
Aldex_family_16sY1_F2F3.all<-data.frame(Aldex_family_16sY1_F2F3.e,Aldex_family_16sY1_F2F3.t)

Aldex_family_16sY2_F1F2.all<-data.frame(Aldex_family_16sY2_F1F2.e,Aldex_family_16sY2_F1F2.t)
Aldex_family_16sY2_F1F3.all<-data.frame(Aldex_family_16sY2_F1F3.e,Aldex_family_16sY2_F1F3.t)
Aldex_family_16sY2_F2F3.all<-data.frame(Aldex_family_16sY2_F2F3.e,Aldex_family_16sY2_F2F3.t)

Aldex_family_ITSY1_F1F2.all<-data.frame(Aldex_family_ITSY1_F1F2.e,Aldex_family_ITSY1_F1F2.t)
Aldex_family_ITSY1_F1F3.all<-data.frame(Aldex_family_ITSY1_F1F3.e,Aldex_family_ITSY1_F1F3.t)
Aldex_family_ITSY1_F2F3.all<-data.frame(Aldex_family_ITSY1_F2F3.e,Aldex_family_ITSY1_F2F3.t)

Aldex_family_ITSY2_F1F2.all<-data.frame(Aldex_family_ITSY2_F1F2.e,Aldex_family_ITSY2_F1F2.t)
Aldex_family_ITSY2_F1F3.all<-data.frame(Aldex_family_ITSY2_F1F3.e,Aldex_family_ITSY2_F1F3.t)
Aldex_family_ITSY2_F2F3.all<-data.frame(Aldex_family_ITSY2_F2F3.e,Aldex_family_ITSY2_F2F3.t)

#Determine which corrected values fall below a threshold
#At family level
Aldex_family_16sY1_F1F2.sig<-which(Aldex_family_16sY1_F1F2.all$wi.eBH <=0.05)
Aldex_family_16sY1_F1F3.sig<-which(Aldex_family_16sY1_F1F3.all$wi.eBH <=0.05)
Aldex_family_16sY1_F2F3.sig<-which(Aldex_family_16sY1_F2F3.all$wi.eBH <=0.05)

Aldex_family_16sY2_F1F2.sig<-which(Aldex_family_16sY2_F1F2.all$wi.eBH <=0.05)
Aldex_family_16sY2_F1F3.sig<-which(Aldex_family_16sY2_F1F3.all$wi.eBH <=0.05)
Aldex_family_16sY2_F2F3.sig<-which(Aldex_family_16sY2_F2F3.all$wi.eBH <=0.05)

Aldex_family_ITSY1_F1F2.sig<-which(Aldex_family_ITSY1_F1F2.all$wi.eBH <=0.05)
Aldex_family_ITSY1_F1F3.sig<-which(Aldex_family_ITSY1_F1F3.all$wi.eBH <=0.05)
Aldex_family_ITSY1_F2F3.sig<-which(Aldex_family_ITSY1_F2F3.all$wi.eBH <=0.05)

Aldex_family_ITSY2_F1F2.sig<-which(Aldex_family_ITSY2_F1F2.all$wi.eBH <=0.05)
Aldex_family_ITSY2_F1F3.sig<-which(Aldex_family_ITSY2_F1F3.all$wi.eBH <=0.05)
Aldex_family_ITSY2_F2F3.sig<-which(Aldex_family_ITSY2_F2F3.all$wi.eBH <=0.05)

#Plot the results - Effect plots described by the documentation

png("Aldex - by Facility - effect plot.png", width = 10, height = 10, units = 'in', res=600)

#16sY1
par(mar=c(4,6,3,4))
par(mfrow=c(4,3))
plot(Aldex_family_16sY1_F1F2.all$diff.win, Aldex_family_16sY1_F1F2.all$diff.btw, pch=19, cex=0.2,
     col='grey', xlab='within group dispersion', ylab='between group difference', main="Bacteria Y1 - F1 v F2")
points(Aldex_family_16sY1_F1F2.all$diff.win[Aldex_family_16sY1_F1F2.sig], 
       Aldex_family_16sY1_F1F2.all$diff.btw[Aldex_family_16sY1_F1F2.sig], pch=19, cex=0.2,col='red')#color significant ones in reg
abline(0,1, col='grey', lty=1)#Add approximate effect size lines
abline(0,-1, col='grey', lty=1)
abline(0,0.5, col='grey', lty=2)
abline(0,-0.5, col='grey', lty=2)

plot(Aldex_family_16sY1_F1F3.all$diff.win, Aldex_family_16sY1_F1F3.all$diff.btw, pch=19, cex=0.2,
     col='grey', xlab="within group dispersion", ylab='between group difference', main="Bacteria Y1 - F1 v F3")
points(Aldex_family_16sY1_F1F3.all$diff.win[Aldex_family_16sY1_F1F3.sig], 
       Aldex_family_16sY1_F1F3.all$diff.btw[Aldex_family_16sY1_F1F3.sig], pch=19, cex=0.2,col='red')#color significant ones in reg
abline(0,1, col='grey', lty=1)#Add approximate effect size lines
abline(0,-1, col='grey', lty=1)
abline(0,0.5, col='grey', lty=2)
abline(0,-0.5, col='grey', lty=2)

plot(Aldex_family_16sY1_F2F3.all$diff.win, Aldex_family_16sY1_F2F3.all$diff.btw, pch=19, cex=0.2,
     col='grey', xlab="within group dispersion", ylab='between group difference', main="Bacteria Y1 - F2 v F3")
points(Aldex_family_16sY1_F2F3.all$diff.win[Aldex_family_16sY1_F2F3.sig], 
       Aldex_family_16sY1_F2F3.all$diff.btw[Aldex_family_16sY1_F2F3.sig], pch=19, cex=0.2,col='red')#color significant ones in reg
abline(0,1, col='grey', lty=1)#Add approximate effect size lines
abline(0,-1, col='grey', lty=1)
abline(0,0.5, col='grey', lty=2)
abline(0,-0.5, col='grey', lty=2)

#16sY2
plot(Aldex_family_16sY2_F1F2.all$diff.win, Aldex_family_16sY2_F1F2.all$diff.btw, pch=19, cex=0.2,
     col='grey', xlab="within group dispersion", ylab='between group difference', main="Bacteria Y2 - F1 v F2")
points(Aldex_family_16sY2_F1F2.all$diff.win[Aldex_family_16sY2_F1F2.sig], 
       Aldex_family_16sY2_F1F2.all$diff.btw[Aldex_family_16sY2_F1F2.sig], pch=19, cex=0.2,col='red')#color significant ones in reg
abline(0,1, col='grey', lty=1)#Add approximate effect size lines
abline(0,-1, col='grey', lty=1)
abline(0,0.5, col='grey', lty=2)
abline(0,-0.5, col='grey', lty=2)

plot(Aldex_family_16sY2_F1F3.all$diff.win, Aldex_family_16sY2_F1F3.all$diff.btw, pch=19, cex=0.2,
     col='grey', xlab="within group dispersion", ylab='between group difference', main="Bacteria Y2 - F1 v F3")
points(Aldex_family_16sY2_F1F3.all$diff.win[Aldex_family_16sY2_F1F3.sig], 
       Aldex_family_16sY2_F1F3.all$diff.btw[Aldex_family_16sY2_F1F3.sig], pch=19, cex=0.2,col='red')#color significant ones in reg
abline(0,1, col='grey', lty=1)#Add approximate effect size lines
abline(0,-1, col='grey', lty=1)
abline(0,0.5, col='grey', lty=2)
abline(0,-0.5, col='grey', lty=2)

plot(Aldex_family_16sY2_F2F3.all$diff.win, Aldex_family_16sY2_F2F3.all$diff.btw, pch=19, cex=0.2,
     col='grey', xlab="within group dispersion", ylab='between group difference', main="Bacteria Y2 - F2 v F3")
points(Aldex_family_16sY2_F2F3.all$diff.win[Aldex_family_16sY2_F2F3.sig], 
       Aldex_family_16sY2_F2F3.all$diff.btw[Aldex_family_16sY2_F2F3.sig], pch=19, cex=0.2,col='red')#color significant ones in reg
abline(0,1, col='grey', lty=1)#Add approximate effect size lines
abline(0,-1, col='grey', lty=1)
abline(0,0.5, col='grey', lty=2)
abline(0,-0.5, col='grey', lty=2)

#ITSY1
plot(Aldex_family_ITSY1_F1F2.all$diff.win, Aldex_family_ITSY1_F1F2.all$diff.btw, pch=19, cex=0.2,
     col='grey', xlab="within group dispersion", ylab='between group difference', main="Fungi Y1 - F1 v F2")
points(Aldex_family_ITSY1_F1F2.all$diff.win[Aldex_family_ITSY1_F1F2.sig], 
       Aldex_family_ITSY1_F1F2.all$diff.btw[Aldex_family_ITSY1_F1F2.sig], pch=19, cex=0.2,col='red')#color significant ones in reg
abline(0,1, col='grey', lty=1)#Add approximate effect size lines
abline(0,-1, col='grey', lty=1)
abline(0,0.5, col='grey', lty=2)
abline(0,-0.5, col='grey', lty=2)

plot(Aldex_family_ITSY1_F1F3.all$diff.win, Aldex_family_ITSY1_F1F3.all$diff.btw, pch=19, cex=0.2,
     col='grey', xlab="within group dispersion", ylab='between group difference', main="Fungi Y1 - F1 v F3")
points(Aldex_family_ITSY1_F1F3.all$diff.win[Aldex_family_ITSY1_F1F3.sig], 
       Aldex_family_ITSY1_F1F3.all$diff.btw[Aldex_family_ITSY1_F1F3.sig], pch=19, cex=0.2,col='red')#color significant ones in reg
abline(0,1, col='grey', lty=1)#Add approximate effect size lines
abline(0,-1, col='grey', lty=1)
abline(0,0.5, col='grey', lty=2)
abline(0,-0.5, col='grey', lty=2)

plot(Aldex_family_ITSY1_F2F3.all$diff.win, Aldex_family_ITSY1_F2F3.all$diff.btw, pch=19, cex=0.2,
     col='grey', xlab="within group dispersion", ylab='between group difference', main="Fungi Y1 - F2 v F3")
points(Aldex_family_ITSY1_F2F3.all$diff.win[Aldex_family_ITSY1_F2F3.sig], 
       Aldex_family_ITSY1_F2F3.all$diff.btw[Aldex_family_ITSY1_F2F3.sig], pch=19, cex=0.2,col='red')#color significant ones in reg
abline(0,1, col='grey', lty=1)#Add approximate effect size lines
abline(0,-1, col='grey', lty=1)
abline(0,0.5, col='grey', lty=2)
abline(0,-0.5, col='grey', lty=2)

#ITSY2
plot(Aldex_family_ITSY2_F1F2.all$diff.win, Aldex_family_ITSY2_F1F2.all$diff.btw, pch=19, cex=0.2,
     col='grey', xlab="within group dispersion", ylab='between group difference', main="Fungi Y2 - F1 v F2")
points(Aldex_family_ITSY2_F1F2.all$diff.win[Aldex_family_ITSY2_F1F2.sig], 
       Aldex_family_ITSY2_F1F2.all$diff.btw[Aldex_family_ITSY2_F1F2.sig], pch=19, cex=0.2,col='red')#color significant ones in reg
abline(0,1, col='grey', lty=1)#Add approximate effect size lines
abline(0,-1, col='grey', lty=1)
abline(0,0.5, col='grey', lty=2)
abline(0,-0.5, col='grey', lty=2)

plot(Aldex_family_ITSY2_F1F3.all$diff.win, Aldex_family_ITSY2_F1F3.all$diff.btw, pch=19, cex=0.2,
     col='grey', xlab="within group dispersion", ylab='between group difference', main="Fungi Y2 - F1 v F3")
points(Aldex_family_ITSY2_F1F3.all$diff.win[Aldex_family_ITSY2_F1F3.sig], 
       Aldex_family_ITSY2_F1F3.all$diff.btw[Aldex_family_ITSY2_F1F3.sig], pch=19, cex=0.2,col='red')#color significant ones in reg
abline(0,1, col='grey', lty=1)#Add approximate effect size lines
abline(0,-1, col='grey', lty=1)
abline(0,0.5, col='grey', lty=2)
abline(0,-0.5, col='grey', lty=2)

plot(Aldex_family_ITSY2_F2F3.all$diff.win, Aldex_family_ITSY2_F2F3.all$diff.btw, pch=19, cex=0.2,
     col='grey', xlab="within group dispersion", ylab='between group difference', main="Fungi Y2 - F2 v F3")
points(Aldex_family_ITSY2_F2F3.all$diff.win[Aldex_family_ITSY2_F2F3.sig], 
       Aldex_family_ITSY2_F2F3.all$diff.btw[Aldex_family_ITSY2_F2F3.sig], pch=19, cex=0.2,col='red')#color significant ones in reg
abline(0,1, col='grey', lty=1)#Add approximate effect size lines
abline(0,-1, col='grey', lty=1)
abline(0,0.5, col='grey', lty=2)
abline(0,-0.5, col='grey', lty=2)

dev.off()


#Plots of significant families relative abundances
#Extract significant OTU data from ALDEx2 output
# 16sY1 
Aldex_family_16sY1_F1F2.sig.row<-rownames(Aldex_family_16sY1_F1F2.all)[which(Aldex_family_16sY1_F1F2.all$wi.eBH <=0.05)] #select the rownames of significant families
# Aldex_family_16sY1_F1F2.sig.table<-subset(Aldex_family_16sY1_F1F2.all, rownames(Aldex_family_16sY1_F1F2.all) %in% Aldex_family_16sY1_F1F2.sig.row) #Subset significant families
# Aldex_family_16sY1_F1F2.sig.otu<-as.data.frame(subset((t(otu.n0.acomp_16sY1)), rownames(t(otu.n0.acomp_16sY1)) %in% Aldex_family_16sY1_F1F2.sig.row)) #Subset table with relative abundance per sample, only for significant fammilies
# Aldex_family_16sY1_F1F2.sig.taxon<-subset(taxon_16sY1, rownames(taxon_16sY1) %in% Aldex_family_16sY1_F1F2.sig.row) #Subset taxon table to get the taxonomy of significant families
# Aldex_family_16sY1_F1F2.sig.table.all<-bind_cols(Aldex_family_16sY1_F1F2.sig.taxon, Aldex_family_16sY1_F1F2.sig.table, Aldex_family_16sY1_F1F2.sig.otu) #combine tables
# Aldex_family_16sY1_F1F2.sig.table.all$Comparison<-rep("F1 v F2", 37)

Aldex_family_16sY1_F1F3.sig.row<-rownames(Aldex_family_16sY1_F1F3.all)[which(Aldex_family_16sY1_F1F3.all$wi.eBH <=0.05)] #select the rownames of significant families
# Aldex_family_16sY1_F1F3.sig.table<-subset(Aldex_family_16sY1_F1F3.all, rownames(Aldex_family_16sY1_F1F3.all) %in% Aldex_family_16sY1_F1F3.sig.row) #Subset significant families
# Aldex_family_16sY1_F1F3.sig.otu<-as.data.frame(subset((t(otu.n0.acomp_16sY1)), rownames(t(otu.n0.acomp_16sY1)) %in% Aldex_family_16sY1_F1F3.sig.row)) #Subset table with relative abundance per sample, only for significant fammilies
# Aldex_family_16sY1_F1F3.sig.taxon<-subset(taxon_16sY1, rownames(taxon_16sY1) %in% Aldex_family_16sY1_F1F3.sig.row) #Subset taxon table to get the taxonomy of significant families
# Aldex_family_16sY1_F1F3.sig.table.all<-bind_cols(Aldex_family_16sY1_F1F3.sig.taxon, Aldex_family_16sY1_F1F3.sig.table, Aldex_family_16sY1_F1F3.sig.otu) #combine tables
# Aldex_family_16sY1_F1F3.sig.table.all$Comparison<-rep("F1 v F3", 39)

Aldex_family_16sY1_F2F3.sig.row<-rownames(Aldex_family_16sY1_F2F3.all)[which(Aldex_family_16sY1_F2F3.all$wi.eBH <=0.05)] #select the rownames of significant families
# Aldex_family_16sY1_F2F3.sig.table<-subset(Aldex_family_16sY1_F2F3.all, rownames(Aldex_family_16sY1_F2F3.all) %in% Aldex_family_16sY1_F2F3.sig.row) #Subset significant families
# Aldex_family_16sY1_F2F3.sig.otu<-as.data.frame(subset((t(otu.n0.acomp_16sY1)), rownames(t(otu.n0.acomp_16sY1)) %in% Aldex_family_16sY1_F2F3.sig.row)) #Subset table with relative abundance per sample, only for significant fammilies
# Aldex_family_16sY1_F2F3.sig.taxon<-subset(taxon_16sY1, rownames(taxon_16sY1) %in% Aldex_family_16sY1_F2F3.sig.row) #Subset taxon table to get the taxonomy of significant families
# Aldex_family_16sY1_F2F3.sig.table.all<-bind_cols(Aldex_family_16sY1_F2F3.sig.taxon, Aldex_family_16sY1_F2F3.sig.table, Aldex_family_16sY1_F2F3.sig.otu) #combine tables
# Aldex_family_16sY1_F2F3.sig.table.all$Comparison<-rep("F3 v F3", 71)

Aldex_family_16sY2_F1F2.sig.row<-rownames(Aldex_family_16sY2_F1F2.all)[which(Aldex_family_16sY2_F1F2.all$wi.eBH <=0.05)] #select the rownames of significant families
Aldex_family_16sY2_F1F3.sig.row<-rownames(Aldex_family_16sY2_F1F3.all)[which(Aldex_family_16sY2_F1F3.all$wi.eBH <=0.05)] #select the rownames of significant families
Aldex_family_16sY2_F2F3.sig.row<-rownames(Aldex_family_16sY2_F2F3.all)[which(Aldex_family_16sY2_F2F3.all$wi.eBH <=0.05)] #select the rownames of significant families

Aldex_family_ITSY1_F1F2.sig.row<-rownames(Aldex_family_ITSY1_F1F2.all)[which(Aldex_family_ITSY1_F1F2.all$wi.eBH <=0.05)] #select the rownames of significant families
Aldex_family_ITSY1_F1F3.sig.row<-rownames(Aldex_family_ITSY1_F1F3.all)[which(Aldex_family_ITSY1_F1F3.all$wi.eBH <=0.05)] #select the rownames of significant families
Aldex_family_ITSY1_F2F3.sig.row<-rownames(Aldex_family_ITSY1_F2F3.all)[which(Aldex_family_ITSY1_F2F3.all$wi.eBH <=0.05)] #select the rownames of significant families

Aldex_family_ITSY2_F1F2.sig.row<-rownames(Aldex_family_ITSY2_F1F2.all)[which(Aldex_family_ITSY2_F1F2.all$wi.eBH <=0.05)] #select the rownames of significant families
Aldex_family_ITSY2_F1F3.sig.row<-rownames(Aldex_family_ITSY2_F1F3.all)[which(Aldex_family_ITSY2_F1F3.all$wi.eBH <=0.05)] #select the rownames of significant families
Aldex_family_ITSY2_F2F3.sig.row<-rownames(Aldex_family_ITSY2_F2F3.all)[which(Aldex_family_ITSY2_F2F3.all$wi.eBH <=0.05)] #select the rownames of significant families

#Combine list of significant OTU names and remove duplicates
Aldex_family_16sY1.sig.row<-unique(c(Aldex_family_16sY1_F1F2.sig.row, Aldex_family_16sY1_F1F3.sig.row, Aldex_family_16sY1_F2F3.sig.row))
Aldex_family_16sY2.sig.row<-unique(c(Aldex_family_16sY2_F1F2.sig.row, Aldex_family_16sY2_F1F3.sig.row, Aldex_family_16sY2_F2F3.sig.row))
Aldex_family_ITSY1.sig.row<-unique(c(Aldex_family_ITSY1_F1F2.sig.row, Aldex_family_ITSY1_F1F3.sig.row, Aldex_family_ITSY1_F2F3.sig.row))
Aldex_family_ITSY2.sig.row<-unique(c(Aldex_family_ITSY2_F1F2.sig.row, Aldex_family_ITSY2_F1F3.sig.row, Aldex_family_ITSY2_F2F3.sig.row))

#Subset table with relative abundance to keep only differential abundant families
SigFamilies_16sY1<-as.data.frame(subset(family_16sY1, family_16sY1$OTU %in% Aldex_family_16sY1.sig.row)) 
SigFamilies_16sY2<-as.data.frame(subset(family_16sY2, family_16sY2$OTU %in% Aldex_family_16sY2.sig.row)) 
SigFamilies_ITSY1<-as.data.frame(subset(family_ITSY1, family_ITSY1$OTU %in% Aldex_family_ITSY1.sig.row)) 
SigFamilies_ITSY2<-as.data.frame(subset(family_ITSY2, family_ITSY2$OTU %in% Aldex_family_ITSY2.sig.row)) 

#Calculate mean, sd and se for each Family by Facility
SigFamilies_16sY1_summarystat <- describeBy(SigFamilies_16sY1$Abundance, list(SigFamilies_16sY1$Family,SigFamilies_16sY1$Facility), mat = TRUE)
SigFamilies_16sY2_summarystat <- describeBy(SigFamilies_16sY2$Abundance, list(SigFamilies_16sY2$Family,SigFamilies_16sY2$Facility), mat = TRUE)
SigFamilies_ITSY1_summarystat <- describeBy(SigFamilies_ITSY1$Abundance, list(SigFamilies_ITSY1$Family,SigFamilies_ITSY1$Facility), mat = TRUE)
SigFamilies_ITSY2_summarystat <- describeBy(SigFamilies_ITSY2$Abundance, list(SigFamilies_ITSY2$Family,SigFamilies_ITSY2$Facility), mat = TRUE)


#Plot Significant families barplots of mean relative abundance by facility

SigFamilies_16sY1_barplot<-ggplot(SigFamilies_16sY1_summarystat, aes(x=group2, y=mean, fill=group2)) + 
  geom_bar(stat="identity", color="black", position=position_dodge()) + 
  facet_wrap(~group1, scales='free_y')+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(.9))+
  theme(axis.title.x = element_blank(), axis.ticks=element_line(color='black'),
  axis.text=element_text(size=12, color='black')) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  ylab("Relative Abundance (%)") + 
  theme(panel.background = element_rect(fill=NA, color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background = element_rect(color='black', fill='white', size=1, linetype = 'solid'))+
  theme(strip.text = element_text(size=10, face='italic'),
        panel.border = element_rect(color="black", fill=NA))+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
    scale_fill_viridis_d(begin = 0.2, end = 0.8, option='inferno')+
  theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))+
  ggtitle("Bacteria - Y1", subtitle = "Differential abundant families")

ggsave("Barplot_SigFam_16sY1.png", plot=SigFamilies_16sY1_barplot, device="png", width=24, height=20, units="in", dpi=600)

SigFamilies_16sY2_barplot<-ggplot(SigFamilies_16sY2_summarystat, aes(x=group2, y=mean, fill=group2)) + 
  geom_bar(stat="identity", color="black", position=position_dodge()) + 
  facet_wrap(~group1, scales='free_y')+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(.9))+
  theme(axis.title.x = element_blank(), axis.ticks=element_line(color='black'),
        axis.text=element_text(size=12, color='black')) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  ylab("Relative Abundance (%)") + 
  theme(panel.background = element_rect(fill=NA, color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background = element_rect(color='black', fill='white', size=1, linetype = 'solid'))+
  theme(strip.text = element_text(size=10, face='italic'),
        panel.border = element_rect(color="black", fill=NA))+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
  scale_fill_viridis_d(begin = 0.2, end = 0.8, option='inferno')+
  theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))+
  ggtitle("Bacteria - Y2", subtitle = "Differential abundant families")

ggsave("Barplot_SigFam_16sY2.png", plot=SigFamilies_16sY2_barplot, device="png", width=24, height=20, units="in", dpi=600)

SigFamilies_ITSY1_barplot<-ggplot(SigFamilies_ITSY1_summarystat, aes(x=group2, y=mean, fill=group2)) + 
  geom_bar(stat="identity", color="black", position=position_dodge()) + 
  facet_wrap(~group1, scales='free_y')+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(.9))+
  theme(axis.title.x = element_blank(), axis.ticks=element_line(color='black'),
        axis.text=element_text(size=12, color='black')) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  ylab("Relative Abundance (%)") + 
  theme(panel.background = element_rect(fill=NA, color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background = element_rect(color='black', fill='white', size=1, linetype = 'solid'))+
  theme(strip.text = element_text(size=10, face='italic'),
        panel.border = element_rect(color="black", fill=NA))+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
  scale_fill_viridis_d(begin = 0.2, end = 0.8, option='inferno')+
  theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))+
  ggtitle("Fungi - Y1", subtitle = "Differential abundant families")

ggsave("Barplot_SigFam_ITSY1.png", plot=SigFamilies_ITSY1_barplot, device="png", width=24, height=20, units="in", dpi=600)

SigFamilies_ITSY2_barplot<-ggplot(SigFamilies_ITSY2_summarystat, aes(x=group2, y=mean, fill=group2)) + 
  geom_bar(stat="identity", color="black", position=position_dodge()) + 
  facet_wrap(~group1, scales='free_y')+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(.9))+
  theme(axis.title.x = element_blank(), axis.ticks=element_line(color='black'),
        axis.text=element_text(size=12, color='black')) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  ylab("Relative Abundance (%)") + 
  theme(panel.background = element_rect(fill=NA, color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background = element_rect(color='black', fill='white', size=1, linetype = 'solid'))+
  theme(strip.text = element_text(size=10, face='italic'),
        panel.border = element_rect(color="black", fill=NA))+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
  scale_fill_viridis_d(begin = 0.2, end = 0.8, option='inferno')+
  theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))+
  ggtitle("Fungi - Y2", subtitle = "Differential abundant families")

ggsave("Barplot_SigFam_ITSY2.png", plot=SigFamilies_ITSY2_barplot, device="png", width=24, height=20, units="in", dpi=600)


#Plots for differential abundant families above 5% relative abundance
Sigfamily_16sY1_over5<-c('Arcobacteraceae','Azospirillaceae','Caulobacteraceae','Chitinophagaceae','Enterobacteriaceae',
                      'Moraxellaceae','Pseudomonadaceae','Rhizobiaceae','Rhodobacteraceae','Sphingobacteriaceae',
                      'Sphingomonadaceae','Xanthomonadaceae')

Sigfamily_16sY2_over5<-c('Arcobacteraceae','Azospirillaceae','Beijerinckiaceae','Microbacteriaceae',
                         'Micrococcaceae','Mycobacteriaceae','Propionibacteriaceae',
                         'Pseudomonadaceae','Rhizobiaceae','Rhodobacteraceae','Sphingomonadaceae')

Sigfamily_ITSY1_over5<-c('Aspergillaceae','Aureobasidiaceae','Bulleribasidiaceae','Cucurbitariaceae','Didymellaceae','Didymosphaeriaceae','Dipodascaceae','Helotiales_fam_Incertae_sedis',
                      'Mycosphaerellaceae','Pleosporaceae','Pleosporales_unclassified','Trichosporonaceae','Ascomycota_unclassified')


Sigfamily_ITSY2_over5<-c('Ascomycota_unclassified','Aspergillaceae','Aureobasidiaceae','Bulleribasidiaceae',
  'Cucurbitariaceae','Didymosphaeriaceae','Dipodascaceae','Mucoraceae','Pichiaceae','Saccharomycetales_fam_Incertae_sedis',
  'Sporidiobolaceae','Trichomeriaceae','Trichosporonaceae','Pleosporales_unclassified')

Sigfamily_16sY1_over5abund <- filter(SigFamilies_16sY1, Family %in% Sigfamily_16sY1_over5)
Sigfamily_16sY2_over5abund <- filter(SigFamilies_16sY2, Family %in% Sigfamily_16sY2_over5)
Sigfamily_ITSY1_over5abund <- filter(SigFamilies_ITSY1, Family %in% Sigfamily_ITSY1_over5)
Sigfamily_ITSY2_over5abund <- filter(SigFamilies_ITSY2, Family %in% Sigfamily_ITSY2_over5)



#Calculate mean, sd and se for each Family by Facility
SigFamilies_16sY1_over5_summarystat <- describeBy(Sigfamily_16sY1_over5abund$Abundance, list(Sigfamily_16sY1_over5abund$Family,Sigfamily_16sY1_over5abund$Facility), mat = TRUE)
SigFamilies_16sY2_over5_summarystat <- describeBy(Sigfamily_16sY2_over5abund$Abundance, list(Sigfamily_16sY2_over5abund$Family,Sigfamily_16sY2_over5abund$Facility), mat = TRUE)
SigFamilies_ITSY1_over5_summarystat <- describeBy(Sigfamily_ITSY1_over5abund$Abundance, list(Sigfamily_ITSY1_over5abund$Family,Sigfamily_ITSY1_over5abund$Facility), mat = TRUE)
SigFamilies_ITSY2_over5_summarystat <- describeBy(Sigfamily_ITSY2_over5abund$Abundance, list(Sigfamily_ITSY2_over5abund$Family,Sigfamily_ITSY2_over5abund$Facility), mat = TRUE)



#Plot Significant families barplots of mean relative abundance by facility

SigFamilies_16sY1_over5_barplot<-ggplot(SigFamilies_16sY1_over5_summarystat, aes(x=group2, y=mean, fill=group2)) + 
  geom_bar(stat="identity", color="black", position=position_dodge()) + 
  facet_wrap(~group1, scales='free_y')+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(.9))+
  theme(axis.title.x = element_blank(), axis.ticks=element_line(color='black'),
        axis.text=element_text(size=12, color='black')) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  ylab("Relative Abundance (%)") + 
  theme(panel.background = element_rect(fill=NA, color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background = element_rect(color='black', fill='white', size=1, linetype = 'solid'))+
  theme(strip.text = element_text(size=10, face='italic'),
        panel.border = element_rect(color="black", fill=NA))+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
  scale_fill_viridis_d(begin = 0.2, end = 0.8, option='inferno')+
  theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))+
  guides(fill=guide_legend(title="Facility"))+
  ggtitle("Bacteria - Y1", subtitle = "Differential abundant families above 5% relative abundace")
ggsave("Barplot_SigFam_over5_16sY1.png", plot=SigFamilies_16sY1_over5_barplot, device="png", width=10, height=8, units="in", dpi=600)

SigFamilies_16sY2_over5_barplot<-ggplot(SigFamilies_16sY2_over5_summarystat, aes(x=group2, y=mean, fill=group2)) + 
  geom_bar(stat="identity", color="black", position=position_dodge()) + 
  facet_wrap(~group1, scales='free_y')+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(.9))+
  theme(axis.title.x = element_blank(), axis.ticks=element_line(color='black'),
        axis.text=element_text(size=12, color='black')) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  ylab("Relative Abundance (%)") + 
  theme(panel.background = element_rect(fill=NA, color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background = element_rect(color='black', fill='white', size=1, linetype = 'solid'))+
  theme(strip.text = element_text(size=10, face='italic'),
        panel.border = element_rect(color="black", fill=NA))+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
  scale_fill_viridis_d(begin = 0.2, end = 0.8, option='inferno')+
  theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))+
  guides(fill=guide_legend(title="Facility"))+
  ggtitle("Bacteria - Y2", subtitle = "Differential abundant families above 5% relative abundace")
ggsave("Barplot_SigFam_over5_16sY2.png", plot=SigFamilies_16sY2_over5_barplot, device="png", width=10, height=8, units="in", dpi=600)

SigFamilies_ITSY1_over5_barplot<-ggplot(SigFamilies_ITSY1_over5_summarystat, aes(x=group2, y=mean, fill=group2)) + 
  geom_bar(stat="identity", color="black", position=position_dodge()) + 
  facet_wrap(~group1, scales='free_y')+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(.9))+
  theme(axis.title.x = element_blank(), axis.ticks=element_line(color='black'),
        axis.text=element_text(size=12, color='black')) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  ylab("Relative Abundance (%)") + 
  theme(panel.background = element_rect(fill=NA, color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background = element_rect(color='black', fill='white', size=1, linetype = 'solid'))+
  theme(strip.text = element_text(size=10, face='italic'),
        panel.border = element_rect(color="black", fill=NA))+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
  scale_fill_viridis_d(begin = 0.2, end = 0.8, option='inferno')+
  theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))+
  guides(fill=guide_legend(title="Facility"))+
  ggtitle("Fungi - Y1", subtitle = "Differential abundant families above 5% relative abundace")
ggsave("Barplot_SigFam_over5_ITSY1.png", plot=SigFamilies_ITSY1_over5_barplot, device="png", width=10, height=8, units="in", dpi=600)

SigFamilies_ITSY2_over5_barplot<-ggplot(SigFamilies_ITSY2_over5_summarystat, aes(x=group2, y=mean, fill=group2)) + 
  geom_bar(stat="identity", color="black", position=position_dodge()) + 
  facet_wrap(~group1, scales='free_y')+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(.9))+
  theme(axis.title.x = element_blank(), axis.ticks=element_line(color='black'),
        axis.text=element_text(size=12, color='black')) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  ylab("Relative Abundance (%)") + 
  theme(panel.background = element_rect(fill=NA, color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background = element_rect(color='black', fill='white', size=1, linetype = 'solid'))+
  theme(strip.text = element_text(size=10, face='italic'),
        panel.border = element_rect(color="black", fill=NA))+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
  scale_fill_viridis_d(begin = 0.2, end = 0.8, option='inferno')+
  theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))+
  guides(fill=guide_legend(title="Facility"))+
  ggtitle("Fungi - Y2", subtitle = "Differential abundant families above 5% relative abundace")
ggsave("Barplot_SigFam_over5_ITSY2.png", plot=SigFamilies_ITSY2_over5_barplot, device="png", width=10, height=8, units="in", dpi=600)



#### DIFFERENTIAL ABUNDANCE BY Lm ####
#All facilities together, Lm +/- is the condition set

#Run for F1 and F3 separatedly
#PERMANOVA showed the facilities were different before, it doesn't make sense to pool the data together.

#FINISHED CODING HERE#


#### DIFFERENTIAL ABUNDANCE BY SEASON ####
#Modify OTU tables to replace the row names from OTU identifier to Family name
#Add Year info to Metadata and merge
#Merge OTU tables for both years

#Run Aldex2 by Year


#### NETWORK ANALYSIS ####




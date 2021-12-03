#Two-year monitoring of Lm in apple packing houses
#Differential abundance analysis at OTU level
#Laura Rolon
#Last updated: 08/19/21

#Load packages
library(ggplot2)
library(ALDEx2)
library(viridis)
library(dendextend)
library(readxl)
library(gplots)
library(BiocParallel)
library(phyloseq)
library(ape)
library(compositions)
library(dplyr)
library(tidyr)
library(psych)
library(LaCroixColoR)
library(cowplot)
library(reshape2)
library(svglite)

set.seed(336)

#Set working directory to where files are located
setwd("/storage/work/m/mlr355/Apple/Differential Abundance")

#### IMPORT DATA ####

##16s
#Import OTU table - 16s
otus_16s<-import_mothur(mothur_shared_file = 'apple16s_132.shared')

#Import taxonomy table -16s
taxon_16s <- as.data.frame(import_mothur(mothur_constaxonomy_file = 'apple16s_132.cons.taxonomy'))
colnames(taxon_16s) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus")
TAX_16s =tax_table(as.matrix(taxon_16s))

#Import metadata -16s
metadata_16s <-read.csv("metadata_apple_16s.csv", header=TRUE, row.names=1)
META_16s = sample_data(metadata_16s)

##ITS
#Import OTU table - ITS
otus_ITS<-import_mothur(mothur_shared_file = 'appleITS.shared')

#Import taxonomy table -ITS
taxon_ITS <- as.data.frame(import_mothur(mothur_constaxonomy_file = 'appleITS.cons.taxonomy'))
colnames(taxon_ITS) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus")
TAX_ITS =tax_table(as.matrix(taxon_ITS))

#Import metadata -ITS
metadata_ITS <-read.csv("metadata_apple_ITS.csv", header=TRUE, row.names=1)
META_ITS = sample_data(metadata_ITS)


#### Differential abundance by LM PRESENCE at OTU level - All FAC and both seasons together ####

#Compare between samples that were positive and negative
Aldex_16s_Lm.clr<-aldex.clr(otus_16s, mc.samples = 128, conds = metadata_16s$L..monocytogenes)
Aldex_ITS_Lm.clr<-aldex.clr(otus_ITS, mc.samples = 128, conds = metadata_ITS$L..monocytogenes)

#Calculate the expected effect size
Aldex_16s_Lm.e<-aldex.effect(Aldex_16s_Lm.clr)
Aldex_ITS_Lm.e<-aldex.effect(Aldex_ITS_Lm.clr)

#Generate the P and Benjamini-Hochberg corrected P values.
Aldex_16s_Lm.t<-aldex.ttest(Aldex_16s_Lm.clr)
Aldex_ITS_Lm.t<-aldex.ttest(Aldex_ITS_Lm.clr)

#Merge data frames
Aldex_16s_Lm.all<-data.frame(Aldex_16s_Lm.e,Aldex_16s_Lm.t)
Aldex_ITS_Lm.all<-data.frame(Aldex_ITS_Lm.e,Aldex_ITS_Lm.t)

#Determine which corrected values fall below a threshold
Aldex_16s_Lm.sig<-which(Aldex_16s_Lm.all$wi.eBH <=0.05)
Aldex_ITS_Lm.sig<-which(Aldex_ITS_Lm.all$wi.eBH <=0.05)

#Plot the results - Effect plots described by the documentation

png("Aldex_Lm_16s_ITS_all.png", width = 10, height = 8, units = 'in', res=600)

#16s
par(mar=c(4,6,3,4))
par(mfrow=c(1,2))
plot(Aldex_16s_Lm.all$diff.win, Aldex_16s_Lm.all$diff.btw, pch=19, cex=0.2,
     col='grey', xlab='within group dispersion', ylab='between group difference', main="Bacteria")
points(Aldex_16s_Lm.all$diff.win[Aldex_16s_Lm.sig], 
       Aldex_16s_Lm.all$diff.btw[Aldex_16s_Lm.sig], pch=19, cex=0.2,col='red')#color significant ones in reg
abline(0,1, col='grey', lty=1)#Add approximate effect size lines
abline(0,-1, col='grey', lty=1)
abline(0,0.5, col='grey', lty=2)
abline(0,-0.5, col='grey', lty=2)

#ITS
plot(Aldex_ITS_Lm.all$diff.win, Aldex_ITS_Lm.all$diff.btw, pch=19, cex=0.2,
     col='grey', xlab='within group dispersion', ylab='between group difference', main="Fungi Y1")
points(Aldex_ITS_Lm.all$diff.win[Aldex_ITS_Lm.sig], 
       Aldex_ITS_Lm.all$diff.btw[Aldex_ITS_Lm.sig], pch=19, cex=0.2,col='red')#color significant ones in reg
abline(0,1, col='grey', lty=1)#Add approximate effect size lines
abline(0,-1, col='grey', lty=1)
abline(0,0.5, col='grey', lty=2)
abline(0,-0.5, col='grey', lty=2)

dev.off()

#Prepare data for plotting differentiall abundant OTUs
Aldex_16s_Lm.sig.row<-rownames(Aldex_16s_Lm.all)[which(Aldex_16s_Lm.all$wi.eBH <=0.05)] #select the rownames of significant families
Aldex_16s_Lm.sig.table<-subset(Aldex_16s_Lm.all, rownames(Aldex_16s_Lm.all) %in% Aldex_16s_Lm.sig.row) #Subset significant families
# Aldex_16s_Lm.sig.otu<-as.data.frame(subset((t(otu.n0.acomp_16s)), rownames(t(otu.n0.acomp_16s)) %in% Aldex_16s_Lm.sig.row)) #Subset table with relative abundance per sample, only for significant fammilies
Aldex_16s_Lm.sig.taxon<-subset(taxon_16s, rownames(taxon_16s) %in% Aldex_16s_Lm.sig.row) #Subset taxon table to get the taxonomy of significant families
Aldex_16s_Lm.sig.table.all<-bind_cols(Aldex_16s_Lm.sig.taxon, Aldex_16s_Lm.sig.table) #combine tables
Aldex_16s_Lm.sig.table.all$Lm<-ifelse(Aldex_16s_Lm.sig.table.all$effect<0, "+", "-") #Add direction of significance
Aldex_16s_Lm.sig.table.all$OTU<-rownames(Aldex_16s_Lm.sig.table.all)

Aldex_ITS_Lm.sig.row<-rownames(Aldex_ITS_Lm.all)[which(Aldex_ITS_Lm.all$wi.eBH <=0.05)] #select the rownames of significant families
Aldex_ITS_Lm.sig.table<-subset(Aldex_ITS_Lm.all, rownames(Aldex_ITS_Lm.all) %in% Aldex_ITS_Lm.sig.row) #Subset significant families
# Aldex_ITS_Lm.sig.otu<-as.data.frame(subset((t(otu.n0.acomp_ITS)), rownames(t(otu.n0.acomp_ITS)) %in% Aldex_ITS_Lm.sig.row)) #Subset table with relative abundance per sample, only for significant fammilies
Aldex_ITS_Lm.sig.taxon<-subset(taxon_ITS, rownames(taxon_ITS) %in% Aldex_ITS_Lm.sig.row) #Subset taxon table to get the taxonomy of significant families
Aldex_ITS_Lm.sig.table.all<-bind_cols(Aldex_ITS_Lm.sig.taxon, Aldex_ITS_Lm.sig.table) #combine tables
Aldex_ITS_Lm.sig.table.all$Lm<-ifelse(Aldex_ITS_Lm.sig.table.all$effect<0, "+", "-") #Add direction of significance
Aldex_ITS_Lm.sig.table.all$OTU<-rownames(Aldex_ITS_Lm.sig.table.all)

#Make vector with OTU names for differential abundant OTUs 
otu_16s_sig_Lm<-Aldex_16s_Lm.sig.table.all$OTU
otu_ITS_sig_Lm<-Aldex_ITS_Lm.sig.table.all$OTU

#Replace zero values with a small non- zero value before clr transformation. 
#Uses Count Zero Multiplicative method to replace zeros and outputs pseudo-counts\
library(zCompositions)
otu.n0_16s<-t(cmultRepl(t(otus_16s), label=0, method="CZM", output="p-counts")) #349844  corrected values
otu.n0_ITS<-t(cmultRepl(t(otus_ITS), label=0, method="CZM", output="p-counts")) #184651  corrected values

#Note: used compositional approach to transform the sample counts to compositions. 
#Transform sample counts into compositions
otu.n0.acomp_16s<-as.data.frame(acomp(t(otu.n0_16s)), total=1)
otu.n0.acomp_ITS<-as.data.frame(acomp(t(otu.n0_ITS)), total=1)

#Make Phyloseq object
phyloseq16s_acomp <- phyloseq(otu_table(otu.n0.acomp_16s, taxa_are_rows = FALSE), TAX_16s, sample_data(metadata_16s))
phyloseqITS_acomp <- phyloseq(otu_table(otu.n0.acomp_ITS, taxa_are_rows = FALSE), TAX_ITS, sample_data(metadata_ITS))

#Make long format table from Phyloseq object
otu_16s_long <- phyloseq16s_acomp %>%  
        transform_sample_counts(function(x) {x * 100} ) %>% #Recalculates the relative abundance 
        psmelt() %>%  #Melts to long format
        arrange(desc(Abundance))

otu_ITS_long <- phyloseqITS_acomp %>%  
        transform_sample_counts(function(x) {x * 100} ) %>% #Recalculates the relative abundance 
        psmelt() %>%  #Melts to long format
        arrange(desc(Abundance))


#Filter table to obtain only significant OTUs
otu_16s_filter_sig_Lm <- filter(otu_16s_long, OTU %in% otu_16s_sig_Lm)
otu_ITS_filter_sig_Lm <- filter(otu_ITS_long, OTU %in% otu_ITS_sig_Lm)

#calculate statistics by OTU and presence of Lm
stat_16s_Lm<-otu_16s_filter_sig_Lm %>%
        group_by(OTU, L..monocytogenes, Family)%>%
        summarise(Mean=mean(Abundance, na.rm=TRUE))

stat_ITS_Lm<-otu_ITS_filter_sig_Lm %>%
        group_by(OTU, L..monocytogenes, Family)%>%
        summarise(Mean=mean(Abundance, na.rm=TRUE))

#Reshape table to wide format
stat_16s_Lm_mean<-dcast(stat_16s_Lm, OTU + Family~ L..monocytogenes , value.var='Mean')
stat_ITS_Lm_mean<-dcast(stat_ITS_Lm, OTU + Family~ L..monocytogenes , value.var='Mean')

#Calculate log fold chage - log2 of the ratio of mean RA of first facility by second facility in comparison
stat_16s_Lm_mean$logchange<-log2(stat_16s_Lm_mean$'-'/stat_16s_Lm_mean$'+')
stat_ITS_Lm_mean$logchange<-log2(stat_ITS_Lm_mean$'-'/stat_ITS_Lm_mean$'+')

#Add Lm identifier 
stat_16s_Lm_mean$L..monocytogenes<-ifelse(stat_16s_Lm_mean$logchange <0 , "+", "-")
stat_ITS_Lm_mean$L..monocytogenes<-ifelse(stat_ITS_Lm_mean$logchange <0 , "+", "-")

#Plot 
#Fig 5A
Logfold_16s_Lm<-ggplot(stat_16s_Lm_mean, aes(x=logchange, y=reorder(OTU,logchange), fill=L..monocytogenes))+
        geom_bar(stat='identity', color='black')+
        geom_text(aes(label=Family), position = position_stack(vjust=0.5), size=2)+
        theme(axis.title.y = element_blank(), axis.ticks=element_line(color='black'),
              axis.text=element_text(size=8, color='black')) + 
        guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
        xlab("Log fold change (log2 Negative/Positive)") + 
        scale_x_continuous(limits=c(-9,8), breaks = c(-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8))+
        theme(panel.background = element_rect(fill=NA, color =NA),
              plot.background = element_rect(fill="transparent", color =NA),
              panel.border = element_rect(color="black", fill=NA)) +
        theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
        scale_fill_manual(values=c('#EA7580'))+
        theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))+
        geom_text(aes(label=Family), position = position_stack(vjust=0.5), size=2)+
        ggtitle("Bacteria - Lm", subtitle = "All Facilities - Y1Y2")

#Fig 6A
Logfold_ITS_Lm<-ggplot(stat_ITS_Lm_mean, aes(x=logchange, y=reorder(OTU,logchange), fill=L..monocytogenes))+
        geom_bar(stat='identity', color='black')+
        geom_text(aes(label=Family), position = position_stack(vjust=0.5), size=2)+
        theme(axis.title.y = element_blank(), axis.ticks=element_line(color='black'),
              axis.text=element_text(size=8, color='black')) + 
        guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
        xlab("Log fold change (log2 Negative/Positive)") + 
        scale_x_continuous(limits=c(-9,8), breaks = c(-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8))+
        theme(panel.background = element_rect(fill=NA, color =NA),
              plot.background = element_rect(fill="transparent", color =NA),
              panel.border = element_rect(color="black", fill=NA)) +
        theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
        scale_fill_manual(values=c('#EA7580'))+
        theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))+
        geom_text(aes(label=Family), position = position_stack(vjust=0.5), size=2)+
        ggtitle("Fungi - Lm", subtitle = "All Facilities - Y1Y2")


#Combine plots

Sig_Lm_all = plot_grid(Logfold_16s_Lm, Logfold_ITS_Lm,  
                       ncol=2, nrow=1, labels = c("A","B"), label_size = 20, vjust = 2, hjust = -1.5)
Sig_Lm_all

ggsave("Aldex_Logfold_LmY1Y2.png", plot=Sig_Lm_all, device="png", width=13, height=11, units="in", dpi=600)
ggsave("Aldex_Logfold_LmY1Y2.svg", plot=Sig_Lm_all, device="svg", width=13, height=11, units="in", dpi=600)


#### Differential abundance by LM PRESENCE at OTU level - Both seasons together excluding F2 ####

#Make Phyloseq with count data
phyloseq16s <- phyloseq(otu_table(otus_16s, taxa_are_rows = TRUE), TAX_16s, sample_data(metadata_16s))
phyloseqITS <- phyloseq(otu_table(otus_ITS, taxa_are_rows = TRUE), TAX_ITS, sample_data(metadata_ITS))

#Prune samples to remove F2
phyloseq16s_F1F3<-subset_samples(phyloseq16s, Facility != "F2")
phyloseqITS_F1F3<-subset_samples(phyloseqITS, Facility != "F2")

#Take otu table for each facility from Phyloseq
otus_16s_F1F3<-as.data.frame(as(otu_table(phyloseq16s_F1F3), "matrix"))
otus_ITS_F1F3<-as.data.frame(as(otu_table(phyloseqITS_F1F3), "matrix"))

#Remove OTUs that have count zero in all samples - Necessary step for the zero imputation function
otus_16s_F1F3<-otus_16s_F1F3[ which(rowSums(otus_16s_F1F3)>0),]
otus_ITS_F1F3<-otus_ITS_F1F3[ which(rowSums(otus_ITS_F1F3)>0),]

#Obtain metadata for F1&F3
metadata_16s_F1F3<-subset(metadata_16s, Facility =="F1" | Facility=="F3")
metadata_ITS_F1F3<-subset(metadata_ITS, Facility =="F1" | Facility=="F3")

#Compare between samples that were positive and negative
Aldex_16s_F1F3_Lm.clr<-aldex.clr(otus_16s_F1F3, mc.samples = 128, conds = metadata_16s_F1F3$L..monocytogenes)
Aldex_ITS_F1F3_Lm.clr<-aldex.clr(otus_ITS_F1F3, mc.samples = 128, conds = metadata_ITS_F1F3$L..monocytogenes)

#Calculate the expected effect size
Aldex_16s_F1F3_Lm.e<-aldex.effect(Aldex_16s_F1F3_Lm.clr)
Aldex_ITS_F1F3_Lm.e<-aldex.effect(Aldex_ITS_F1F3_Lm.clr)

#Generate the P and Benjamini-Hochberg corrected P values.
Aldex_16s_F1F3_Lm.t<-aldex.ttest(Aldex_16s_F1F3_Lm.clr)
Aldex_ITS_F1F3_Lm.t<-aldex.ttest(Aldex_ITS_F1F3_Lm.clr)

#Merge data frames
Aldex_16s_F1F3_Lm.all<-data.frame(Aldex_16s_F1F3_Lm.e,Aldex_16s_F1F3_Lm.t)
Aldex_ITS_F1F3_Lm.all<-data.frame(Aldex_ITS_F1F3_Lm.e,Aldex_ITS_F1F3_Lm.t)

#Determine which corrected values fall below a threshold
Aldex_16s_F1F3_Lm.sig<-which(Aldex_16s_F1F3_Lm.all$wi.eBH <=0.05)
Aldex_ITS_F1F3_Lm.sig<-which(Aldex_ITS_F1F3_Lm.all$wi.eBH <=0.05)


#### Differential abundance by LM PRESENCE at OTU level - All FAC and split by season ####

#Split OTU table by year
#Make Phyloseq with count data
phyloseq16s <- phyloseq(otu_table(otus_16s, taxa_are_rows = TRUE), TAX_16s, META_16s)
phyloseqITS <- phyloseq(otu_table(otus_ITS, taxa_are_rows = TRUE), TAX_ITS, META_ITS)

# Subset Phyloseq for each year
physeq_16sY1 <- subset_samples(phyloseq16s, Year == "Y1") 
physeq_16sY2 <- subset_samples(phyloseq16s, Year == "Y2") 

physeq_ITSY1 <- subset_samples(phyloseqITS, Year == "Y1") 
physeq_ITSY2 <- subset_samples(phyloseqITS, Year == "Y2") 

#Take otu table for each facility from Phyloseq
otus_16sY1<-as.data.frame(as(otu_table(physeq_16sY1), "matrix"))
otus_16sY2<-as.data.frame(as(otu_table(physeq_16sY2), "matrix"))

otus_ITSY1<-as.data.frame(as(otu_table(physeq_ITSY1), "matrix"))
otus_ITSY2<-as.data.frame(as(otu_table(physeq_ITSY2), "matrix"))

#Remove OTUs that have count zero in all samples - Necessary step for the zero imputation function
otus_16sY1<-otus_16sY1[ which(rowSums(otus_16sY1)>0),]
otus_16sY2<-otus_16sY2[ which(rowSums(otus_16sY2)>0),]

otus_ITSY1<-otus_ITSY1[ which(rowSums(otus_ITSY1)>0),]
otus_ITSY2<-otus_ITSY2[ which(rowSums(otus_ITSY2)>0),]

#Subset metadata
metadata_16sY1<-subset(metadata_16s, Year=="Y1")
metadata_16sY2<-subset(metadata_16s, Year=="Y2")

metadata_ITSY1<-subset(metadata_ITS, Year=="Y1")
metadata_ITSY2<-subset(metadata_ITS, Year=="Y2")

#Compare between samples that were positive and negative
Aldex_16sY1_Lm.clr<-aldex.clr(otus_16sY1, mc.samples = 128, conds = metadata_16sY1$L..monocytogenes)
Aldex_16sY2_Lm.clr<-aldex.clr(otus_16sY2, mc.samples = 128, conds = metadata_16sY2$L..monocytogenes)
Aldex_ITSY1_Lm.clr<-aldex.clr(otus_ITSY1, mc.samples = 128, conds = metadata_ITSY1$L..monocytogenes)
Aldex_ITSY2_Lm.clr<-aldex.clr(otus_ITSY2, mc.samples = 128, conds = metadata_ITSY2$L..monocytogenes)

#Calculate the expected effect size
Aldex_16sY1_Lm.e<-aldex.effect(Aldex_16sY1_Lm.clr)
Aldex_16sY2_Lm.e<-aldex.effect(Aldex_16sY2_Lm.clr)
Aldex_ITSY1_Lm.e<-aldex.effect(Aldex_ITSY1_Lm.clr)
Aldex_ITSY2_Lm.e<-aldex.effect(Aldex_ITSY2_Lm.clr)

#Generate the P and Benjamini-Hochberg corrected P values.
Aldex_16sY1_Lm.t<-aldex.ttest(Aldex_16sY1_Lm.clr)
Aldex_16sY2_Lm.t<-aldex.ttest(Aldex_16sY2_Lm.clr)
Aldex_ITSY1_Lm.t<-aldex.ttest(Aldex_ITSY1_Lm.clr)
Aldex_ITSY2_Lm.t<-aldex.ttest(Aldex_ITSY2_Lm.clr)

#Merge data frames
Aldex_16sY1_Lm.all<-data.frame(Aldex_16sY1_Lm.e,Aldex_16sY1_Lm.t)
Aldex_16sY2_Lm.all<-data.frame(Aldex_16sY2_Lm.e,Aldex_16sY2_Lm.t)
Aldex_ITSY1_Lm.all<-data.frame(Aldex_ITSY1_Lm.e,Aldex_ITSY1_Lm.t)
Aldex_ITSY2_Lm.all<-data.frame(Aldex_ITSY2_Lm.e,Aldex_ITSY2_Lm.t)

#Determine which corrected values fall below a threshold
Aldex_16sY1_Lm.sig<-which(Aldex_16sY1_Lm.all$wi.eBH <=0.05)
Aldex_16sY2_Lm.sig<-which(Aldex_16sY2_Lm.all$wi.eBH <=0.05) #No significant OTU
Aldex_ITSY1_Lm.sig<-which(Aldex_ITSY1_Lm.all$wi.eBH <=0.05)
Aldex_ITSY2_Lm.sig<-which(Aldex_ITSY2_Lm.all$wi.eBH <=0.05) #No significant OTU

#Plot the results - Effect plots described by the documentation

png("Aldex_Lm_16sY1_16sY2_ITSY1_ITSY2_all facilities.png", width = 10, height = 8, units = 'in', res=600)

#16sY1
par(mar=c(4,6,3,4))
par(mfrow=c(2,2))
plot(Aldex_16sY1_Lm.all$diff.win, Aldex_16sY1_Lm.all$diff.btw, pch=19, cex=0.2,
     col='grey', xlab='within group dispersion', ylab='between group difference', main="Bacteria Y1")
points(Aldex_16sY1_Lm.all$diff.win[Aldex_16sY1_Lm.sig], 
       Aldex_16sY1_Lm.all$diff.btw[Aldex_16sY1_Lm.sig], pch=19, cex=0.2,col='red')#color significant ones in reg
abline(0,1, col='grey', lty=1)#Add approximate effect size lines
abline(0,-1, col='grey', lty=1)
abline(0,0.5, col='grey', lty=2)
abline(0,-0.5, col='grey', lty=2)

#16sY2
plot(Aldex_16sY2_Lm.all$diff.win, Aldex_16sY2_Lm.all$diff.btw, pch=19, cex=0.2,
     col='grey', xlab='within group dispersion', ylab='between group difference', main="Bacteria Y2")
points(Aldex_16sY2_Lm.all$diff.win[Aldex_16sY2_Lm.sig], 
       Aldex_16sY2_Lm.all$diff.btw[Aldex_16sY2_Lm.sig], pch=19, cex=0.2,col='red')#color significant ones in reg
abline(0,1, col='grey', lty=1)#Add approximate effect size lines
abline(0,-1, col='grey', lty=1)
abline(0,0.5, col='grey', lty=2)
abline(0,-0.5, col='grey', lty=2)

#ITSY1
plot(Aldex_ITSY1_Lm.all$diff.win, Aldex_ITSY1_Lm.all$diff.btw, pch=19, cex=0.2,
     col='grey', xlab='within group dispersion', ylab='between group difference', main="Fungi Y1")
points(Aldex_ITSY1_Lm.all$diff.win[Aldex_ITSY1_Lm.sig], 
       Aldex_ITSY1_Lm.all$diff.btw[Aldex_ITSY1_Lm.sig], pch=19, cex=0.2,col='red')#color significant ones in reg
abline(0,1, col='grey', lty=1)#Add approximate effect size lines
abline(0,-1, col='grey', lty=1)
abline(0,0.5, col='grey', lty=2)
abline(0,-0.5, col='grey', lty=2)

#ITSY2
plot(Aldex_ITSY2_Lm.all$diff.win, Aldex_ITSY2_Lm.all$diff.btw, pch=19, cex=0.2,
     col='grey', xlab='within group dispersion', ylab='between group difference', main="Fungi Y2")
points(Aldex_ITSY2_Lm.all$diff.win[Aldex_ITSY2_Lm.sig], 
       Aldex_ITSY2_Lm.all$diff.btw[Aldex_ITSY2_Lm.sig], pch=19, cex=0.2,col='red')#color significant ones in reg
abline(0,1, col='grey', lty=1)#Add approximate effect size lines
abline(0,-1, col='grey', lty=1)
abline(0,0.5, col='grey', lty=2)
abline(0,-0.5, col='grey', lty=2)

dev.off()


Aldex_16sY1_Lm.sig.row<-rownames(Aldex_16sY1_Lm.all)[which(Aldex_16sY1_Lm.all$wi.eBH <=0.05)] #select the rownames of significant families
Aldex_16sY1_Lm.sig.table<-subset(Aldex_16sY1_Lm.all, rownames(Aldex_16sY1_Lm.all) %in% Aldex_16sY1_Lm.sig.row) #Subset significant families
# Aldex_16sY1_Lm.sig.otu<-as.data.frame(subset((t(otu.n0.acomp_16sY1)), rownames(t(otu.n0.acomp_16sY1)) %in% Aldex_16sY1_Lm.sig.row)) #Subset table with relative abundance per sample, only for significant fammilies
Aldex_16sY1_Lm.sig.taxon<-subset(taxon_16s, rownames(taxon_16s) %in% Aldex_16sY1_Lm.sig.row) #Subset taxon table to get the taxonomy of significant families
Aldex_16sY1_Lm.sig.table.all<-bind_cols(Aldex_16sY1_Lm.sig.taxon, Aldex_16sY1_Lm.sig.table) #combine tables
Aldex_16sY1_Lm.sig.table.all$Lm<-ifelse(Aldex_16sY1_Lm.sig.table.all$effect<0, "+", "-") #Add direction of significance
Aldex_16sY1_Lm.sig.table.all$OTU<-rownames(Aldex_16sY1_Lm.sig.table.all)

Aldex_ITSY1_Lm.sig.row<-rownames(Aldex_ITSY1_Lm.all)[which(Aldex_ITSY1_Lm.all$wi.eBH <=0.05)] #select the rownames of significant families
Aldex_ITSY1_Lm.sig.table<-subset(Aldex_ITSY1_Lm.all, rownames(Aldex_ITSY1_Lm.all) %in% Aldex_ITSY1_Lm.sig.row) #Subset significant families
# Aldex_ITSY1_Lm.sig.otu<-as.data.frame(subset((t(otu.n0.acomp_ITSY1)), rownames(t(otu.n0.acomp_ITSY1)) %in% Aldex_ITSY1_Lm.sig.row)) #Subset table with relative abundance per sample, only for significant fammilies
Aldex_ITSY1_Lm.sig.taxon<-subset(taxon_ITS, rownames(taxon_ITS) %in% Aldex_ITSY1_Lm.sig.row) #Subset taxon table to get the taxonomy of significant families
Aldex_ITSY1_Lm.sig.table.all<-bind_cols(Aldex_ITSY1_Lm.sig.taxon, Aldex_ITSY1_Lm.sig.table) #combine tables
Aldex_ITSY1_Lm.sig.table.all$Lm<-ifelse(Aldex_ITSY1_Lm.sig.table.all$effect<0, "+", "-") #Add direction of significance
Aldex_ITSY1_Lm.sig.table.all$OTU<-rownames(Aldex_ITSY1_Lm.sig.table.all)


#Make vector with OTU names for differential abundant OTUs with effect size over 1 for each pair of comparisons
otu_16sY1_sig_Lm<-Aldex_16sY1_Lm.sig.table.all$OTU
otu_ITSY1_sig_Lm<-Aldex_ITSY1_Lm.sig.table.all$OTU


# Replace zero values with a small non- zero value before clr transformation. 
#Uses Count Zero Multiplicative method to replace zeros and outputs pseudo-counts\
otu.n0_16sY1<-t(cmultRepl(t(otus_16sY1), label=0, method="CZM", output="p-counts")) #349844  corrected values
otu.n0_ITSY1<-t(cmultRepl(t(otus_ITSY1), label=0, method="CZM", output="p-counts")) #184651  corrected values

#Note: used compositional approach to transform the sample counts to compositions. 
#Transform sample counts into compositions
otu.n0.acomp_16sY1<-as.data.frame(acomp(t(otu.n0_16sY1)), total=1)
otu.n0.acomp_ITSY1<-as.data.frame(acomp(t(otu.n0_ITSY1)), total=1)

#Make Phyloseq object
phyloseq16sY1_acomp <- phyloseq(otu_table(otu.n0.acomp_16sY1, taxa_are_rows = FALSE), TAX_16s, sample_data(metadata_16sY1))
phyloseqITSY1_acomp <- phyloseq(otu_table(otu.n0.acomp_ITSY1, taxa_are_rows = FALSE), TAX_ITS, sample_data(metadata_ITSY1))

#Make long format table from Phyloseq object
otu_16sY1_long <- phyloseq16sY1_acomp %>%  
        transform_sample_counts(function(x) {x * 100} ) %>% #Recalculates the relative abundance 
        psmelt() %>%  #Melts to long format
        arrange(desc(Abundance))

otu_ITSY1_long <- phyloseqITSY1_acomp %>%  
        transform_sample_counts(function(x) {x * 100} ) %>% #Recalculates the relative abundance 
        psmelt() %>%  #Melts to long format
        arrange(desc(Abundance))

#Filter table to obtain only significant OTUs
otu_16sY1_filter_sig_Lm <- filter(otu_16sY1_long, OTU %in% otu_16sY1_sig_Lm)
otu_ITSY1_filter_sig_Lm <- filter(otu_ITSY1_long, OTU %in% otu_ITSY1_sig_Lm)

#calculate statistics by OTU and presence of Lm
stat_16sY1_Lm<-otu_16sY1_filter_sig_Lm %>%
        group_by(OTU, L..monocytogenes, Family)%>%
        summarise(Mean=mean(Abundance, na.rm=TRUE))

stat_ITSY1_Lm<-otu_ITSY1_filter_sig_Lm %>%
        group_by(OTU, L..monocytogenes, Family)%>%
        summarise(Mean=mean(Abundance, na.rm=TRUE))

#Reshape table to wide format
stat_16sY1_Lm_mean<-dcast(stat_16sY1_Lm, OTU + Family~ L..monocytogenes , value.var='Mean')
stat_ITSY1_Lm_mean<-dcast(stat_ITSY1_Lm, OTU + Family~ L..monocytogenes , value.var='Mean')

#Calculate log fold chage - log2 of the ratio of mean RA of first facility by second facility in comparison
stat_16sY1_Lm_mean$logchange<-log2(stat_16sY1_Lm_mean$'-'/stat_16sY1_Lm_mean$'+')
stat_ITSY1_Lm_mean$logchange<-log2(stat_ITSY1_Lm_mean$'-'/stat_ITSY1_Lm_mean$'+')

#Add Lm identifier 
stat_16sY1_Lm_mean$L..monocytogenes<-ifelse(stat_16sY1_Lm_mean$logchange <0 , "+", "-")
stat_ITSY1_Lm_mean$L..monocytogenes<-ifelse(stat_ITSY1_Lm_mean$logchange <0 , "+", "-")

#Plot 
#Fig 5C
Logfold_16sY1_Lm<-ggplot(stat_16sY1_Lm_mean, aes(x=logchange, y=reorder(OTU,logchange), fill=L..monocytogenes))+
        geom_bar(stat='identity', color='black')+
        geom_text(aes(label=Family), position = position_stack(vjust=0.5), size=2)+
        theme(axis.title.y = element_blank(), axis.ticks=element_line(color='black'),
              axis.text=element_text(size=8, color='black')) + 
        guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
        xlab("Log fold change (log2 Negative/Positive)") + 
        scale_x_continuous(limits=c(-9,8), breaks = c(-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8))+
        theme(panel.background = element_rect(fill=NA, color =NA),
              plot.background = element_rect(fill="transparent", color =NA),
              panel.border = element_rect(color="black", fill=NA)) +
        theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
        scale_fill_manual(values=c('#F8CD9C','#EA7580'))+
        theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))+
        geom_text(aes(label=Family), position = position_stack(vjust=0.5), size=2)+
        ggtitle("Bacteria Y1 - Lm", subtitle = "All facilities Y1")

#Fig 6C
Logfold_ITSY1_Lm<-ggplot(stat_ITSY1_Lm_mean, aes(x=logchange, y=reorder(OTU,logchange), fill=L..monocytogenes))+
        geom_bar(stat='identity', color='black')+
        geom_text(aes(label=Family), position = position_stack(vjust=0.5), size=2)+
        theme(axis.title.y = element_blank(), axis.ticks=element_line(color='black'),
              axis.text=element_text(size=8, color='black')) + 
        guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
        xlab("Log fold change (log2 Negative/Positive)") + 
        scale_x_continuous(limits=c(-9,8), breaks = c(-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8))+
        theme(panel.background = element_rect(fill=NA, color =NA),
              plot.background = element_rect(fill="transparent", color =NA),
              panel.border = element_rect(color="black", fill=NA)) +
        theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
        scale_fill_manual(values=c('#F8CD9C','#EA7580'))+
        theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))+
        geom_text(aes(label=Family), position = position_stack(vjust=0.5), size=2)+
        ggtitle("Fungi Y1 - Lm", subtitle = "All Facilities Y1")


#Combine plots

Sig_Lm_all = plot_grid(Logfold_16sY1_Lm, Logfold_ITSY1_Lm,  
                          ncol=2, nrow=1, labels = c("A","B"), label_size = 20, vjust = 2, hjust = -1.5)
Sig_Lm_all

ggsave("Aldex_Logfold_Lm_by year.png", plot=Sig_Lm_all, device="png", width=13, height=11, units="in", dpi=600)
ggsave("Aldex_Logfold_Lm_by year.svg", plot=Sig_Lm_all, device="svg", width=13, height=11, units="in", dpi=600)


#### Differential abundance by LM PRESENCE at OTU level - Split by season and excluding F2 ####


#Compare between samples that were positive and negative
Aldex_16sY1_F1Lm.clr<-aldex.clr(F116sY1_otu, mc.samples = 128, conds = metadata_F116sY1$L..monocytogenes)
Aldex_16sY1_F3Lm.clr<-aldex.clr(F316sY1_otu, mc.samples = 128, conds = metadata_F316sY1$L..monocytogenes)

Aldex_16sY2_F1Lm.clr<-aldex.clr(F116sY2_otu, mc.samples = 128, conds = metadata_F116sY2$L..monocytogenes)
Aldex_16sY2_F3Lm.clr<-aldex.clr(F316sY2_otu, mc.samples = 128, conds = metadata_F316sY2$L..monocytogenes)

Aldex_ITSY1_F1Lm.clr<-aldex.clr(F1ITSY1_otu, mc.samples = 128, conds = metadata_F1ITSY1$L..monocytogenes)
Aldex_ITSY1_F3Lm.clr<-aldex.clr(F3ITSY1_otu, mc.samples = 128, conds = metadata_F3ITSY1$L..monocytogenes)

Aldex_ITSY2_F1Lm.clr<-aldex.clr(F1ITSY2_otu, mc.samples = 128, conds = metadata_F1ITSY2$L..monocytogenes)
Aldex_ITSY2_F3Lm.clr<-aldex.clr(F3ITSY2_otu, mc.samples = 128, conds = metadata_F3ITSY2$L..monocytogenes)

#Calculate the expected effect size
Aldex_16sY1_F1Lm.e<-aldex.effect(Aldex_16sY1_F1Lm.clr)
Aldex_16sY1_F3Lm.e<-aldex.effect(Aldex_16sY1_F3Lm.clr)

Aldex_16sY2_F1Lm.e<-aldex.effect(Aldex_16sY2_F1Lm.clr)
Aldex_16sY2_F3Lm.e<-aldex.effect(Aldex_16sY2_F3Lm.clr)

Aldex_ITSY1_F1Lm.e<-aldex.effect(Aldex_ITSY1_F1Lm.clr)
Aldex_ITSY1_F3Lm.e<-aldex.effect(Aldex_ITSY1_F3Lm.clr)

Aldex_ITSY2_F1Lm.e<-aldex.effect(Aldex_ITSY2_F1Lm.clr)
Aldex_ITSY2_F3Lm.e<-aldex.effect(Aldex_ITSY2_F3Lm.clr)

#Generate the P and Benjamini-Hochberg corrected P values.
Aldex_16sY1_F1Lm.t<-aldex.ttest(Aldex_16sY1_F1Lm.clr)
Aldex_16sY1_F3Lm.t<-aldex.ttest(Aldex_16sY1_F3Lm.clr)

Aldex_16sY2_F1Lm.t<-aldex.ttest(Aldex_16sY2_F1Lm.clr)
Aldex_16sY2_F3Lm.t<-aldex.ttest(Aldex_16sY2_F3Lm.clr)

Aldex_ITSY1_F1Lm.t<-aldex.ttest(Aldex_ITSY1_F1Lm.clr)
Aldex_ITSY1_F3Lm.t<-aldex.ttest(Aldex_ITSY1_F3Lm.clr)

Aldex_ITSY2_F1Lm.t<-aldex.ttest(Aldex_ITSY2_F1Lm.clr)
Aldex_ITSY2_F3Lm.t<-aldex.ttest(Aldex_ITSY2_F3Lm.clr)

#Merge data frames
Aldex_16sY1_F1Lm.all<-data.frame(Aldex_16sY1_F1Lm.e,Aldex_16sY1_F1Lm.t)
Aldex_16sY1_F3Lm.all<-data.frame(Aldex_16sY1_F3Lm.e,Aldex_16sY1_F3Lm.t)

Aldex_16sY2_F1Lm.all<-data.frame(Aldex_16sY2_F1Lm.e,Aldex_16sY2_F1Lm.t)
Aldex_16sY2_F3Lm.all<-data.frame(Aldex_16sY2_F3Lm.e,Aldex_16sY2_F3Lm.t)

Aldex_ITSY1_F1Lm.all<-data.frame(Aldex_ITSY1_F1Lm.e,Aldex_ITSY1_F1Lm.t)
Aldex_ITSY1_F3Lm.all<-data.frame(Aldex_ITSY1_F3Lm.e,Aldex_ITSY1_F3Lm.t)

Aldex_ITSY2_F1Lm.all<-data.frame(Aldex_ITSY2_F1Lm.e,Aldex_ITSY2_F1Lm.t)
Aldex_ITSY2_F3Lm.all<-data.frame(Aldex_ITSY2_F3Lm.e,Aldex_ITSY2_F3Lm.t)

#Determine which corrected values fall below a threshold
Aldex_16sY1_F1Lm.sig<-which(Aldex_16sY1_F1Lm.all$wi.eBH <=0.05)
Aldex_16sY1_F3Lm.sig<-which(Aldex_16sY1_F1Lm.all$wi.eBH <=0.05)

Aldex_16sY2_F1Lm.sig<-which(Aldex_16sY2_F1Lm.all$wi.eBH <=0.05)
Aldex_16sY2_F3Lm.sig<-which(Aldex_16sY2_F1Lm.all$wi.eBH <=0.05)

Aldex_ITSY1_F1Lm.sig<-which(Aldex_ITSY1_F1Lm.all$wi.eBH <=0.05)
Aldex_ITSY1_F3Lm.sig<-which(Aldex_ITSY1_F1Lm.all$wi.eBH <=0.05)

Aldex_ITSY2_F1Lm.sig<-which(Aldex_ITSY2_F1Lm.all$wi.eBH <=0.05)
Aldex_ITSY2_F3Lm.sig<-which(Aldex_ITSY2_F1Lm.all$wi.eBH <=0.05)

#Plot the results - Effect plots described by the documentation

png("Aldex - by Lm presence F1 and F3 - effect plot.png", width = 10, height = 8, units = 'in', res=600)

#16sY1
par(mfrow=c(2,4))
plot(Aldex_16sY1_F1Lm.all$diff.win, Aldex_16sY1_F1Lm.all$diff.btw, pch=19, cex=0.2,
     col='grey', xlab='within group dispersion', ylab='between group difference', main="Bacteria Y1 - F1")
points(Aldex_16sY1_F1Lm.all$diff.win[Aldex_16sY1_F1Lm.sig], 
       Aldex_16sY1_F1Lm.all$diff.btw[Aldex_16sY1_F1Lm.sig], pch=19, cex=0.2,col='red')#color significant ones in reg
abline(0,1, col='grey', lty=1)#Add approximate effect size lines
abline(0,-1, col='grey', lty=1)
abline(0,0.5, col='grey', lty=2)
abline(0,-0.5, col='grey', lty=2)

plot(Aldex_16sY1_F3Lm.all$diff.win, Aldex_16sY1_F3Lm.all$diff.btw, pch=19, cex=0.2,
     col='grey', xlab="within group dispersion", ylab='between group difference', main="Bacteria Y1 - F3 ")
points(Aldex_16sY1_F3Lm.all$diff.win[Aldex_16sY1_F3Lm.sig], 
       Aldex_16sY1_F3Lm.all$diff.btw[Aldex_16sY1_F3Lm.sig], pch=19, cex=0.2,col='red')#color significant ones in reg
abline(0,1, col='grey', lty=1)#Add approximate effect size lines
abline(0,-1, col='grey', lty=1)
abline(0,0.5, col='grey', lty=2)
abline(0,-0.5, col='grey', lty=2)

#16sY2
plot(Aldex_16sY2_F1Lm.all$diff.win, Aldex_16sY2_F1Lm.all$diff.btw, pch=19, cex=0.2,
     col='grey', xlab='within group dispersion', ylab='between group difference', main="Bacteria Y2 - F1")
points(Aldex_16sY2_F1Lm.all$diff.win[Aldex_16sY2_F1Lm.sig], 
       Aldex_16sY2_F1Lm.all$diff.btw[Aldex_16sY2_F1Lm.sig], pch=19, cex=0.2,col='red')#color significant ones in reg
abline(0,1, col='grey', lty=1)#Add approximate effect size lines
abline(0,-1, col='grey', lty=1)
abline(0,0.5, col='grey', lty=2)
abline(0,-0.5, col='grey', lty=2)

plot(Aldex_16sY2_F3Lm.all$diff.win, Aldex_16sY2_F3Lm.all$diff.btw, pch=19, cex=0.2,
     col='grey', xlab="within group dispersion", ylab='between group difference', main="Bacteria Y2 - F3 ")
points(Aldex_16sY2_F3Lm.all$diff.win[Aldex_16sY2_F3Lm.sig], 
       Aldex_16sY2_F3Lm.all$diff.btw[Aldex_16sY2_F3Lm.sig], pch=19, cex=0.2,col='red')#color significant ones in reg
abline(0,1, col='grey', lty=1)#Add approximate effect size lines
abline(0,-1, col='grey', lty=1)
abline(0,0.5, col='grey', lty=2)
abline(0,-0.5, col='grey', lty=2)

#ITSY1
plot(Aldex_ITSY1_F1Lm.all$diff.win, Aldex_ITSY1_F1Lm.all$diff.btw, pch=19, cex=0.2,
     col='grey', xlab='within group dispersion', ylab='between group difference', main="Fungi Y1 - F1")
points(Aldex_ITSY1_F1Lm.all$diff.win[Aldex_ITSY1_F1Lm.sig], 
       Aldex_ITSY1_F1Lm.all$diff.btw[Aldex_ITSY1_F1Lm.sig], pch=19, cex=0.2,col='red')#color significant ones in reg
abline(0,1, col='grey', lty=1)#Add approximate effect size lines
abline(0,-1, col='grey', lty=1)
abline(0,0.5, col='grey', lty=2)
abline(0,-0.5, col='grey', lty=2)

plot(Aldex_ITSY1_F3Lm.all$diff.win, Aldex_ITSY1_F3Lm.all$diff.btw, pch=19, cex=0.2,
     col='grey', xlab="within group dispersion", ylab='between group difference', main="Fungi Y1 - F3 ")
points(Aldex_ITSY1_F3Lm.all$diff.win[Aldex_ITSY1_F3Lm.sig], 
       Aldex_ITSY1_F3Lm.all$diff.btw[Aldex_ITSY1_F3Lm.sig], pch=19, cex=0.2,col='red')#color significant ones in reg
abline(0,1, col='grey', lty=1)#Add approximate effect size lines
abline(0,-1, col='grey', lty=1)
abline(0,0.5, col='grey', lty=2)
abline(0,-0.5, col='grey', lty=2)

#ITSY2
plot(Aldex_ITSY2_F1Lm.all$diff.win, Aldex_ITSY2_F1Lm.all$diff.btw, pch=19, cex=0.2,
     col='grey', xlab='within group dispersion', ylab='between group difference', main="Fungi Y2 - F1")
points(Aldex_ITSY2_F1Lm.all$diff.win[Aldex_ITSY2_F1Lm.sig], 
       Aldex_ITSY2_F1Lm.all$diff.btw[Aldex_ITSY2_F1Lm.sig], pch=19, cex=0.2,col='red')#color significant ones in reg
abline(0,1, col='grey', lty=1)#Add approximate effect size lines
abline(0,-1, col='grey', lty=1)
abline(0,0.5, col='grey', lty=2)
abline(0,-0.5, col='grey', lty=2)

plot(Aldex_ITSY2_F3Lm.all$diff.win, Aldex_ITSY2_F3Lm.all$diff.btw, pch=19, cex=0.2,
     col='grey', xlab="within group dispersion", ylab='between group difference', main="Fungi Y2 - F3 ")
points(Aldex_ITSY2_F3Lm.all$diff.win[Aldex_ITSY2_F3Lm.sig], 
       Aldex_ITSY2_F3Lm.all$diff.btw[Aldex_ITSY2_F3Lm.sig], pch=19, cex=0.2,col='red')#color significant ones in reg
abline(0,1, col='grey', lty=1)#Add approximate effect size lines
abline(0,-1, col='grey', lty=1)
abline(0,0.5, col='grey', lty=2)
abline(0,-0.5, col='grey', lty=2)

dev.off()









#### Differential abundance for each facility by year at the OTU level ####

#Split OTU table by year
#Make Phyloseq with count data
OTU_16s <- otu_table(otus_16s, taxa_are_rows = TRUE)
phyloseq_16s = phyloseq(OTU_16s, TAX_16s, META_16s)
TREE_16s = rtree(ntaxa(phyloseq_16s), rooted=TRUE, tip.label = taxa_names(phyloseq_16s))  
phyloseq16s <- phyloseq(OTU_16s, TAX_16s, TREE_16s, META_16s)

OTU_ITS <- otu_table(otus_ITS, taxa_are_rows = TRUE)
phyloseq_ITS = phyloseq(OTU_ITS, TAX_ITS, META_ITS)
TREE_ITS = rtree(ntaxa(phyloseq_ITS), rooted=TRUE, tip.label = taxa_names(phyloseq_ITS))  
phyloseqITS <- phyloseq(OTU_ITS, TAX_ITS, TREE_ITS, META_ITS)

# Subset Phyloseq for each year
physeq_16sF1 <- subset_samples(phyloseq16s, Facility == "F1") 
physeq_16sF2 <- subset_samples(phyloseq16s, Facility == "F2") 
physeq_16sF3 <- subset_samples(phyloseq16s, Facility == "F3") 

physeq_ITSF1 <- subset_samples(phyloseqITS, Facility == "F1") 
physeq_ITSF2 <- subset_samples(phyloseqITS, Facility == "F2") 
physeq_ITSF3 <- subset_samples(phyloseqITS, Facility == "F3") 

#Take otu table for each facility from Phyloseq
otus_16sF1<-as.data.frame(as(otu_table(physeq_16sF1), "matrix"))
otus_16sF2<-as.data.frame(as(otu_table(physeq_16sF2), "matrix"))
otus_16sF3<-as.data.frame(as(otu_table(physeq_16sF3), "matrix"))

otus_ITSF1<-as.data.frame(as(otu_table(physeq_ITSF1), "matrix"))
otus_ITSF2<-as.data.frame(as(otu_table(physeq_ITSF2), "matrix"))
otus_ITSF3<-as.data.frame(as(otu_table(physeq_ITSF3), "matrix"))

#Remove OTUs that have count zero in all samples - Necessary step for the zero imputation function
otus_16sF1<-otus_16sF1[ which(rowSums(otus_16sF1)>0),]
otus_16sF2<-otus_16sF2[ which(rowSums(otus_16sF2)>0),]
otus_16sF3<-otus_16sF3[ which(rowSums(otus_16sF3)>0),]

otus_ITSF1<-otus_ITSF1[ which(rowSums(otus_ITSF1)>0),]
otus_ITSF2<-otus_ITSF2[ which(rowSums(otus_ITSF2)>0),]
otus_ITSF3<-otus_ITSF3[ which(rowSums(otus_ITSF3)>0),]

#Subset metadata
metadata_16sF1<-subset(metadata_16s, Facility == "F1")
metadata_16sF2<-subset(metadata_16s, Facility == "F2")
metadata_16sF3<-subset(metadata_16s, Facility == "F3")

metadata_ITSF1<-subset(metadata_ITS, Facility == "F1")
metadata_ITSF2<-subset(metadata_ITS, Facility == "F2")
metadata_ITSF3<-subset(metadata_ITS, Facility == "F3")


#Change Year to character vector - as.factor affects aldex.effect function
metadata_16sF1$Year <- as.character(metadata_16sF1$Year)
metadata_16sF2$Year <- as.character(metadata_16sF2$Year)
metadata_16sF3$Year <- as.character(metadata_16sF3$Year)

metadata_ITSF1$Year <- as.character(metadata_ITSF1$Year)
metadata_ITSF2$Year <- as.character(metadata_ITSF2$Year)
metadata_ITSF3$Year <- as.character(metadata_ITSF3$Year)


#Aldex2
#Generate 128 Dirichlet distributed Monte Carlo instances and center-log ratio transform them
#OTU table needs to be as OTU in rows. Only two conditions can be compared at a time.

Aldex_16sF1.clr<-aldex.clr(otus_16sF1, mc.samples = 128, conds = metadata_16sF1$Year)
Aldex_16sF2.clr<-aldex.clr(otus_16sF2, mc.samples = 128, conds = metadata_16sF2$Year)
Aldex_16sF3.clr<-aldex.clr(otus_16sF3, mc.samples = 128, conds = metadata_16sF3$Year)

Aldex_ITSF1.clr<-aldex.clr(otus_ITSF1, mc.samples = 128, conds = metadata_ITSF1$Year)
Aldex_ITSF2.clr<-aldex.clr(otus_ITSF2, mc.samples = 128, conds = metadata_ITSF2$Year)
Aldex_ITSF3.clr<-aldex.clr(otus_ITSF3, mc.samples = 128, conds = metadata_ITSF3$Year)


#Calculate the expected effect size
Aldex_16sF1.e<-aldex.effect(Aldex_16sF1.clr)
Aldex_16sF2.e<-aldex.effect(Aldex_16sF2.clr)
Aldex_16sF3.e<-aldex.effect(Aldex_16sF3.clr)

Aldex_ITSF1.e<-aldex.effect(Aldex_ITSF1.clr)
Aldex_ITSF2.e<-aldex.effect(Aldex_ITSF2.clr)
Aldex_ITSF3.e<-aldex.effect(Aldex_ITSF3.clr)

#Generate the P and Benjamini-Hochberg corrected P values.
Aldex_16sF1.t<-aldex.ttest(Aldex_16sF1.clr)
Aldex_16sF2.t<-aldex.ttest(Aldex_16sF2.clr)
Aldex_16sF3.t<-aldex.ttest(Aldex_16sF3.clr)

Aldex_ITSF1.t<-aldex.ttest(Aldex_ITSF1.clr)
Aldex_ITSF2.t<-aldex.ttest(Aldex_ITSF2.clr)
Aldex_ITSF3.t<-aldex.ttest(Aldex_ITSF3.clr)


#Merge data frames
Aldex_16sF1.all<-data.frame(Aldex_16sF1.e,Aldex_16sF1.t)
Aldex_16sF2.all<-data.frame(Aldex_16sF2.e,Aldex_16sF2.t)
Aldex_16sF3.all<-data.frame(Aldex_16sF3.e,Aldex_16sF3.t)

Aldex_ITSF1.all<-data.frame(Aldex_ITSF1.e,Aldex_ITSF1.t)
Aldex_ITSF2.all<-data.frame(Aldex_ITSF2.e,Aldex_ITSF2.t)
Aldex_ITSF3.all<-data.frame(Aldex_ITSF3.e,Aldex_ITSF3.t)


#Determine which corrected values fall below a threshold
Aldex_16sF1.sig<-which(Aldex_16sF1.all$wi.eBH <=0.05)
Aldex_16sF2.sig<-which(Aldex_16sF2.all$wi.eBH <=0.05)
Aldex_16sF3.sig<-which(Aldex_16sF3.all$wi.eBH <=0.05)

Aldex_ITSF1.sig<-which(Aldex_ITSF1.all$wi.eBH <=0.05)
Aldex_ITSF2.sig<-which(Aldex_ITSF2.all$wi.eBH <=0.05)
Aldex_ITSF3.sig<-which(Aldex_ITSF3.all$wi.eBH <=0.05)


#Plot the results - Effect plots described by the documentation

png("Aldex by Year - effect plot.png", width = 8, height= 10, units = 'in', res=600)

#16sY1
par(mar=c(4,6,3,4))
par(mfrow=c(2,3))
plot(Aldex_16sF1.all$diff.win, Aldex_16sF1.all$diff.btw, pch=19, cex=0.2,
     col='grey', xlab='within group dispersion', ylab='between group difference', main="Bacteria F1 - Y1 v Y2")
points(Aldex_16sF1.all$diff.win[Aldex_16sF1.sig], 
       Aldex_16sF1.all$diff.btw[Aldex_16sF1.sig], pch=19, cex=0.2,col='red')#color significant ones in reg
abline(0,1, col='grey', lty=1)#Add approximate effect size lines
abline(0,-1, col='grey', lty=1)
abline(0,0.5, col='grey', lty=2)
abline(0,-0.5, col='grey', lty=2)

plot(Aldex_16sF2.all$diff.win, Aldex_16sF2.all$diff.btw, pch=19, cex=0.2,
     col='grey', xlab='within group dispersion', ylab='between group difference', main="Bacteria F2 - Y1 v Y2")
points(Aldex_16sF2.all$diff.win[Aldex_16sF2.sig], 
       Aldex_16sF2.all$diff.btw[Aldex_16sF2.sig], pch=19, cex=0.2,col='red')#color significant ones in reg
abline(0,1, col='grey', lty=1)#Add approximate effect size lines
abline(0,-1, col='grey', lty=1)
abline(0,0.5, col='grey', lty=2)
abline(0,-0.5, col='grey', lty=2)

plot(Aldex_16sF3.all$diff.win, Aldex_16sF3.all$diff.btw, pch=19, cex=0.2,
     col='grey', xlab='within group dispersion', ylab='between group difference', main="Bacteria F3 - Y1 v Y2")
points(Aldex_16sF3.all$diff.win[Aldex_16sF3.sig], 
       Aldex_16sF3.all$diff.btw[Aldex_16sF3.sig], pch=19, cex=0.2,col='red')#color significant ones in reg
abline(0,1, col='grey', lty=1)#Add approximate effect size lines
abline(0,-1, col='grey', lty=1)
abline(0,0.5, col='grey', lty=2)
abline(0,-0.5, col='grey', lty=2)

plot(Aldex_ITSF1.all$diff.win, Aldex_ITSF1.all$diff.btw, pch=19, cex=0.2,
     col='grey', xlab='within group dispersion', ylab='between group difference', main="Fungi F1 - Y1 v Y2")
points(Aldex_ITSF1.all$diff.win[Aldex_ITSF1.sig], 
       Aldex_ITSF1.all$diff.btw[Aldex_ITSF1.sig], pch=19, cex=0.2,col='red')#color significant ones in reg
abline(0,1, col='grey', lty=1)#Add approximate effect size lines
abline(0,-1, col='grey', lty=1)
abline(0,0.5, col='grey', lty=2)
abline(0,-0.5, col='grey', lty=2)

plot(Aldex_ITSF2.all$diff.win, Aldex_ITSF2.all$diff.btw, pch=19, cex=0.2,
     col='grey', xlab='within group dispersion', ylab='between group difference', main="Fungi F2 - Y1 v Y2")
points(Aldex_ITSF2.all$diff.win[Aldex_ITSF2.sig], 
       Aldex_ITSF2.all$diff.btw[Aldex_ITSF2.sig], pch=19, cex=0.2,col='red')#color significant ones in reg
abline(0,1, col='grey', lty=1)#Add approximate effect size lines
abline(0,-1, col='grey', lty=1)
abline(0,0.5, col='grey', lty=2)
abline(0,-0.5, col='grey', lty=2)

plot(Aldex_ITSF3.all$diff.win, Aldex_ITSF3.all$diff.btw, pch=19, cex=0.2,
     col='grey', xlab='within group dispersion', ylab='between group difference', main="Fungi F3 - Y1 v Y2")
points(Aldex_ITSF3.all$diff.win[Aldex_ITSF3.sig], 
       Aldex_ITSF3.all$diff.btw[Aldex_ITSF3.sig], pch=19, cex=0.2,col='red')#color significant ones in reg
abline(0,1, col='grey', lty=1)#Add approximate effect size lines
abline(0,-1, col='grey', lty=1)
abline(0,0.5, col='grey', lty=2)
abline(0,-0.5, col='grey', lty=2)

dev.off()

#Plots of significant families relative abundances
#Extract significant OTU data from ALDEx2 output
# 16s
Aldex_16sF1.sig.row<-rownames(Aldex_16sF1.all)[which(Aldex_16sF1.all$wi.eBH <=0.05)] #select the rownames of significant families
Aldex_16sF1.sig.table<-subset(Aldex_16sF1.all, rownames(Aldex_16sF1.all) %in% Aldex_16sF1.sig.row) #Subset significant families
# Aldex_16sF1.sig.otu<-as.data.frame(subset((t(otu.n0.acomp_16sY1)), rownames(t(otu.n0.acomp_16sY1)) %in% Aldex_16sF1.sig.row)) #Subset table with relative abundance per sample, only for significant fammilies
Aldex_16sF1.sig.taxon<-subset(taxon_16s, rownames(taxon_16s) %in% Aldex_16sF1.sig.row) #Subset taxon table to get the taxonomy of significant families
Aldex_16sF1.sig.table.all<-bind_cols(Aldex_16sF1.sig.taxon, Aldex_16sF1.sig.table) #combine tables
Aldex_16sF1.sig.table.all$OTU<-rownames(Aldex_16sF1.sig.table.all)

Aldex_16sF2.sig.row<-rownames(Aldex_16sF2.all)[which(Aldex_16sF2.all$wi.eBH <=0.05)] #select the rownames of significant families
Aldex_16sF2.sig.table<-subset(Aldex_16sF2.all, rownames(Aldex_16sF2.all) %in% Aldex_16sF2.sig.row) #Subset significant families
# Aldex_16sF2.sig.otu<-as.data.frame(subset((t(otu.n0.acomp_16sY1)), rownames(t(otu.n0.acomp_16sY1)) %in% Aldex_16sF2.sig.row)) #Subset table with relative abundance per sample, only for significant fammilies
Aldex_16sF2.sig.taxon<-subset(taxon_16s, rownames(taxon_16s) %in% Aldex_16sF2.sig.row) #Subset taxon table to get the taxonomy of significant families
Aldex_16sF2.sig.table.all<-bind_cols(Aldex_16sF2.sig.taxon, Aldex_16sF2.sig.table) #combine tables
Aldex_16sF2.sig.table.all$OTU<-rownames(Aldex_16sF2.sig.table.all)

Aldex_16sF3.sig.row<-rownames(Aldex_16sF3.all)[which(Aldex_16sF3.all$wi.eBH <=0.05)] #select the rownames of significant families
Aldex_16sF3.sig.table<-subset(Aldex_16sF3.all, rownames(Aldex_16sF3.all) %in% Aldex_16sF3.sig.row) #Subset significant families
# Aldex_16sF3.sig.otu<-as.data.frame(subset((t(otu.n0.acomp_16sY1)), rownames(t(otu.n0.acomp_16sY1)) %in% Aldex_16sF3.sig.row)) #Subset table with relative abundance per sample, only for significant fammilies
Aldex_16sF3.sig.taxon<-subset(taxon_16s, rownames(taxon_16s) %in% Aldex_16sF3.sig.row) #Subset taxon table to get the taxonomy of significant families
Aldex_16sF3.sig.table.all<-bind_cols(Aldex_16sF3.sig.taxon, Aldex_16sF3.sig.table) #combine tables
Aldex_16sF3.sig.table.all$OTU<-rownames(Aldex_16sF3.sig.table.all)

# ITS
Aldex_ITSF1.sig.row<-rownames(Aldex_ITSF1.all)[which(Aldex_ITSF1.all$wi.eBH <=0.05)] #select the rownames of significant families
Aldex_ITSF1.sig.table<-subset(Aldex_ITSF1.all, rownames(Aldex_ITSF1.all) %in% Aldex_ITSF1.sig.row) #Subset significant families
# Aldex_ITSF1.sig.otu<-as.data.frame(subset((t(otu.n0.acomp_ITSY1)), rownames(t(otu.n0.acomp_ITSY1)) %in% Aldex_ITSF1.sig.row)) #Subset table with relative abundance per sample, only for significant fammilies
Aldex_ITSF1.sig.taxon<-subset(taxon_ITS, rownames(taxon_ITS) %in% Aldex_ITSF1.sig.row) #Subset taxon table to get the taxonomy of significant families
Aldex_ITSF1.sig.table.all<-bind_cols(Aldex_ITSF1.sig.taxon, Aldex_ITSF1.sig.table) #combine tables
Aldex_ITSF1.sig.table.all$OTU<-rownames(Aldex_ITSF1.sig.table.all)

Aldex_ITSF2.sig.row<-rownames(Aldex_ITSF2.all)[which(Aldex_ITSF2.all$wi.eBH <=0.05)] #select the rownames of significant families
Aldex_ITSF2.sig.table<-subset(Aldex_ITSF2.all, rownames(Aldex_ITSF2.all) %in% Aldex_ITSF2.sig.row) #Subset significant families
# Aldex_ITSF2.sig.otu<-as.data.frame(subset((t(otu.n0.acomp_ITSY1)), rownames(t(otu.n0.acomp_ITSY1)) %in% Aldex_ITSF2.sig.row)) #Subset table with relative abundance per sample, only for significant fammilies
Aldex_ITSF2.sig.taxon<-subset(taxon_ITS, rownames(taxon_ITS) %in% Aldex_ITSF2.sig.row) #Subset taxon table to get the taxonomy of significant families
Aldex_ITSF2.sig.table.all<-bind_cols(Aldex_ITSF2.sig.taxon, Aldex_ITSF2.sig.table) #combine tables
Aldex_ITSF2.sig.table.all$OTU<-rownames(Aldex_ITSF2.sig.table.all)

Aldex_ITSF3.sig.row<-rownames(Aldex_ITSF3.all)[which(Aldex_ITSF3.all$wi.eBH <=0.05)] #select the rownames of significant families
Aldex_ITSF3.sig.table<-subset(Aldex_ITSF3.all, rownames(Aldex_ITSF3.all) %in% Aldex_ITSF3.sig.row) #Subset significant families
# Aldex_ITSF3.sig.otu<-as.data.frame(subset((t(otu.n0.acomp_ITSY1)), rownames(t(otu.n0.acomp_ITSY1)) %in% Aldex_ITSF3.sig.row)) #Subset table with relative abundance per sample, only for significant fammilies
Aldex_ITSF3.sig.taxon<-subset(taxon_ITS, rownames(taxon_ITS) %in% Aldex_ITSF3.sig.row) #Subset taxon table to get the taxonomy of significant families
Aldex_ITSF3.sig.table.all<-bind_cols(Aldex_ITSF3.sig.taxon, Aldex_ITSF3.sig.table) #combine tables
Aldex_ITSF3.sig.table.all$OTU<-rownames(Aldex_ITSF3.sig.table.all)



#Effect size plots for significant OTUs (p<0.05)

#16s
Sig_16sF1<-ggplot(Aldex_16sF1.sig.table.all, aes(x=effect, y=reorder(OTU,effect)))+
        geom_bar(stat='identity', color='black', fill="grey80")+
        geom_vline(xintercept = c(-1,1), color="grey50", linetype='dashed')+
        theme(axis.title.y = element_blank(), axis.ticks=element_line(color='black'),
              axis.text=element_text(size=5, color='black')) + 
        guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
        xlab("Effect size") + 
        theme(panel.background = element_rect(fill=NA, color =NA),
              plot.background = element_rect(fill="transparent", color =NA)) +
        theme(strip.background = element_rect(color='black', fill='white', size=1, linetype = 'solid'))+
        theme(strip.text = element_text(size=10, face='italic', angle=0),
              panel.border = element_rect(color="black", fill=NA))+
        theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
        theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))+
        ggtitle("Bacteria F1 - Y1 v Y2", subtitle = "Differential abundant OTUs")

Sig_16sF2<-ggplot(Aldex_16sF2.sig.table.all, aes(x=effect, y=reorder(OTU,effect)))+
        geom_bar(stat='identity', color='black', fill="grey80")+
        geom_vline(xintercept = c(-1,1), color="grey50", linetype='dashed')+
        theme(axis.title.y = element_blank(), axis.ticks=element_line(color='black'),
              axis.text=element_text(size=5, color='black')) + 
        guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
        xlab("Effect size") + 
        theme(panel.background = element_rect(fill=NA, color =NA),
              plot.background = element_rect(fill="transparent", color =NA)) +
        theme(strip.background = element_rect(color='black', fill='white', size=1, linetype = 'solid'))+
        theme(strip.text = element_text(size=10, face='italic', angle=0),
              panel.border = element_rect(color="black", fill=NA))+
        theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
        theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))+
        ggtitle("Bacteria F2 - Y1 v Y2", subtitle = "Differential abundant OTUs")

Sig_16sF3<-ggplot(Aldex_16sF3.sig.table.all, aes(x=effect, y=reorder(OTU,effect)))+
        geom_bar(stat='identity', color='black', fill="grey80")+
        geom_vline(xintercept = c(-1,1), color="grey50", linetype='dashed')+
        theme(axis.title.y = element_blank(), axis.ticks=element_line(color='black'),
              axis.text=element_text(size=5, color='black')) + 
        guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
        xlab("Effect size") + 
        theme(panel.background = element_rect(fill=NA, color =NA),
              plot.background = element_rect(fill="transparent", color =NA)) +
        theme(strip.background = element_rect(color='black', fill='white', size=1, linetype = 'solid'))+
        theme(strip.text = element_text(size=10, face='italic', angle=0),
              panel.border = element_rect(color="black", fill=NA))+
        theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
        theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))+
        ggtitle("Bacteria F3 - Y1 v Y2", subtitle = "Differential abundant OTUs")

#ITS
Sig_ITSF1<-ggplot(Aldex_ITSF1.sig.table.all, aes(x=effect, y=reorder(OTU,effect)))+
        geom_bar(stat='identity', color='black', fill="grey80")+
        geom_vline(xintercept = c(-1,1), color="grey50", linetype='dashed')+
        theme(axis.title.y = element_blank(), axis.ticks=element_line(color='black'),
              axis.text=element_text(size=5, color='black')) + 
        guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
        xlab("Effect size") + 
        theme(panel.background = element_rect(fill=NA, color =NA),
              plot.background = element_rect(fill="transparent", color =NA)) +
        theme(strip.background = element_rect(color='black', fill='white', size=1, linetype = 'solid'))+
        theme(strip.text = element_text(size=10, face='italic', angle=0),
              panel.border = element_rect(color="black", fill=NA))+
        theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
        theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))+
        ggtitle("Fungi F1 - Y1 v Y2", subtitle = "Differential abundant OTUs")

Sig_ITSF2<-ggplot(Aldex_ITSF2.sig.table.all, aes(x=effect, y=reorder(OTU,effect)))+
        geom_bar(stat='identity', color='black', fill="grey80")+
        geom_vline(xintercept = c(-1,1), color="grey50", linetype='dashed')+
        theme(axis.title.y = element_blank(), axis.ticks=element_line(color='black'),
              axis.text=element_text(size=5, color='black')) + 
        guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
        xlab("Effect size") + 
        theme(panel.background = element_rect(fill=NA, color =NA),
              plot.background = element_rect(fill="transparent", color =NA)) +
        theme(strip.background = element_rect(color='black', fill='white', size=1, linetype = 'solid'))+
        theme(strip.text = element_text(size=10, face='italic', angle=0),
              panel.border = element_rect(color="black", fill=NA))+
        theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
        theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))+
        ggtitle("Fungi F2 - Y1 v Y2", subtitle = "Differential abundant OTUs")

Sig_ITSF3<-ggplot(Aldex_ITSF3.sig.table.all, aes(x=effect, y=reorder(OTU,effect)))+
        geom_bar(stat='identity', color='black', fill="grey80")+
        geom_vline(xintercept = c(-1,1), color="grey50", linetype='dashed')+
        theme(axis.title.y = element_blank(), axis.ticks=element_line(color='black'),
              axis.text=element_text(size=5, color='black')) + 
        guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
        xlab("Effect size") + 
        theme(panel.background = element_rect(fill=NA, color =NA),
              plot.background = element_rect(fill="transparent", color =NA)) +
        theme(strip.background = element_rect(color='black', fill='white', size=1, linetype = 'solid'))+
        theme(strip.text = element_text(size=10, face='italic', angle=0),
              panel.border = element_rect(color="black", fill=NA))+
        theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
        theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))+
        ggtitle("Fungi F3 - Y1 v Y2", subtitle = "Differential abundant OTUs")



#Effect plots for those OTU with effect size over 1 and less than -1
Aldex_16sF1.eff.table.all<-bind_rows(subset(Aldex_16sF1.sig.table.all, effect >=1),subset(Aldex_16sF1.sig.table.all, effect <=-1))
Aldex_16sF2.eff.table.all<-bind_rows(subset(Aldex_16sF2.sig.table.all, effect >=1),subset(Aldex_16sF2.sig.table.all, effect <=-1))
Aldex_16sF3.eff.table.all<-bind_rows(subset(Aldex_16sF3.sig.table.all, effect >=1),subset(Aldex_16sF3.sig.table.all, effect <=-1))

Aldex_ITSF1.eff.table.all<-bind_rows(subset(Aldex_ITSF1.sig.table.all, effect >=1),subset(Aldex_ITSF1.sig.table.all, effect <=-1))
Aldex_ITSF2.eff.table.all<-bind_rows(subset(Aldex_ITSF2.sig.table.all, effect >=1),subset(Aldex_ITSF2.sig.table.all, effect <=-1))
Aldex_ITSF3.eff.table.all<-bind_rows(subset(Aldex_ITSF3.sig.table.all, effect >=1),subset(Aldex_ITSF3.sig.table.all, effect <=-1))


## Boxplots of OTUs with effect size grater than 1

#Remove OTUs that have count zero in all samples - Necessary step for the zero imputation function
otus_16sF1<-otus_16sF1[ which(rowSums(otus_16sF1)>0),]
otus_16sF2<-otus_16sF2[ which(rowSums(otus_16sF2)>0),]
otus_16sF3<-otus_16sF3[ which(rowSums(otus_16sF3)>0),]

otus_ITSF1<-otus_ITSF1[ which(rowSums(otus_ITSF1)>0),]
otus_ITSF2<-otus_ITSF2[ which(rowSums(otus_ITSF2)>0),]
otus_ITSF3<-otus_ITSF3[ which(rowSums(otus_ITSF3)>0),]

#Step 1: Convert OTU table to appropriate format
#Following step requires samples on rows and OTUs in columns
head(t(otus_16sF1))
head(t(otus_16sF2))
head(t(otus_16sF3))

head(t(otus_ITSF1))
head(t(otus_ITSF2))
head(t(otus_ITSF3))

#Step 2: Replace zero values with a small non- zero value before clr transformation. 
#Uses Count Zero Multiplicative method to replace zeros and outputs pseudo-counts\
library(zCompositions)
otu.n0_16sF1<-t(cmultRepl(t(otus_16sF1), label=0, method="CZM", output="p-counts"))
otu.n0_16sF2<-t(cmultRepl(t(otus_16sF2), label=0, method="CZM", output="p-counts"))
otu.n0_16sF3<-t(cmultRepl(t(otus_16sF3), label=0, method="CZM", output="p-counts"))

otu.n0_ITSF1<-t(cmultRepl(t(otus_ITSF1), label=0, method="CZM", output="p-counts"))
otu.n0_ITSF2<-t(cmultRepl(t(otus_ITSF2), label=0, method="CZM", output="p-counts"))
otu.n0_ITSF3<-t(cmultRepl(t(otus_ITSF3), label=0, method="CZM", output="p-counts"))

#Transform metadata into Phyloseq format
META_16sF1 = sample_data(metadata_16sF1)
META_16sF2 = sample_data(metadata_16sF2)
META_16sF3 = sample_data(metadata_16sF3)

META_ITSF1 = sample_data(metadata_ITSF1)
META_ITSF2 = sample_data(metadata_ITSF2)
META_ITSF3 = sample_data(metadata_ITSF3)

#Note: used compositional approach to transform the sample counts to compositions. 

# #Transform sample counts into compositions
otu.n0.acomp_16sF1<-as.data.frame(acomp(t(otu.n0_16sF1)), total=1)
otu.n0.acomp_16sF2<-as.data.frame(acomp(t(otu.n0_16sF2)), total=1)
otu.n0.acomp_16sF3<-as.data.frame(acomp(t(otu.n0_16sF3)), total=1)

otu.n0.acomp_ITSF1<-as.data.frame(acomp(t(otu.n0_ITSF1)), total=1)
otu.n0.acomp_ITSF2<-as.data.frame(acomp(t(otu.n0_ITSF2)), total=1)
otu.n0.acomp_ITSF3<-as.data.frame(acomp(t(otu.n0_ITSF3)), total=1)

#OTU level plots

#Make Phyloseq object
OTU_16sF1 <- otu_table(otu.n0.acomp_16sF1, taxa_are_rows = FALSE)
phyloseq_16sF1 = phyloseq(OTU_16sF1, TAX_16s, META_16sF1)
TREE_16sF1 = rtree(ntaxa(phyloseq_16sF1), rooted=TRUE, tip.label = taxa_names(phyloseq_16sF1))  
phyloseq16sF1 <- phyloseq(OTU_16sF1, TAX_16s, TREE_16sF1, META_16sF1)

OTU_16sF2 <- otu_table(otu.n0.acomp_16sF2, taxa_are_rows = FALSE)
phyloseq_16sF2 = phyloseq(OTU_16sF2, TAX_16s, META_16sF2)
TREE_16sF2 = rtree(ntaxa(phyloseq_16sF2), rooted=TRUE, tip.label = taxa_names(phyloseq_16sF2))  
phyloseq16sF2 <- phyloseq(OTU_16sF2, TAX_16s, TREE_16sF2, META_16sF2)

OTU_16sF3 <- otu_table(otu.n0.acomp_16sF3, taxa_are_rows = FALSE)
phyloseq_16sF3 = phyloseq(OTU_16sF3, TAX_16s, META_16sF3)
TREE_16sF3 = rtree(ntaxa(phyloseq_16sF3), rooted=TRUE, tip.label = taxa_names(phyloseq_16sF3))  
phyloseq16sF3 <- phyloseq(OTU_16sF3, TAX_16s, TREE_16sF3, META_16sF3)

OTU_ITSF1 <- otu_table(otu.n0.acomp_ITSF1, taxa_are_rows = FALSE)
phyloseq_ITSF1 = phyloseq(OTU_ITSF1, TAX_ITS, META_ITSF1)
TREE_ITSF1 = rtree(ntaxa(phyloseq_ITSF1), rooted=TRUE, tip.label = taxa_names(phyloseq_ITSF1))  
phyloseqITSF1 <- phyloseq(OTU_ITSF1, TAX_ITS, TREE_ITSF1, META_ITSF1)

OTU_ITSF2 <- otu_table(otu.n0.acomp_ITSF2, taxa_are_rows = FALSE)
phyloseq_ITSF2 = phyloseq(OTU_ITSF2, TAX_ITS, META_ITSF2)
TREE_ITSF2 = rtree(ntaxa(phyloseq_ITSF2), rooted=TRUE, tip.label = taxa_names(phyloseq_ITSF2))  
phyloseqITSF2 <- phyloseq(OTU_ITSF2, TAX_ITS, TREE_ITSF2, META_ITSF2)

OTU_ITSF3 <- otu_table(otu.n0.acomp_ITSF3, taxa_are_rows = FALSE)
phyloseq_ITSF3 = phyloseq(OTU_ITSF3, TAX_ITS, META_ITSF3)
TREE_ITSF3 = rtree(ntaxa(phyloseq_ITSF3), rooted=TRUE, tip.label = taxa_names(phyloseq_ITSF3))  
phyloseqITSF3 <- phyloseq(OTU_ITSF3, TAX_ITS, TREE_ITSF3, META_ITSF3)


#OTU level Plots
#Make long format table from Phyloseq object
otu_16sF1_long <- phyloseq16sF1 %>%  
        transform_sample_counts(function(x) {x * 100} ) %>% #Recalculates the relative abundance 
        psmelt() %>%  #Melts to long format
        arrange(desc(Abundance))

otu_16sF2_long <- phyloseq16sF2 %>%  
        transform_sample_counts(function(x) {x * 100} ) %>% #Recalculates the relative abundance 
        psmelt() %>%  #Melts to long format
        arrange(desc(Abundance))

otu_16sF3_long <- phyloseq16sF3 %>%  
        transform_sample_counts(function(x) {x * 100} ) %>% #Recalculates the relative abundance 
        psmelt() %>%  #Melts to long format
        arrange(desc(Abundance))

otu_ITSF1_long <- phyloseqITSF1 %>%  
        transform_sample_counts(function(x) {x * 100} ) %>% #Recalculates the relative abundance 
        psmelt() %>%  #Melts to long format
        arrange(desc(Abundance))

otu_ITSF2_long <- phyloseqITSF2 %>%  
        transform_sample_counts(function(x) {x * 100} ) %>% #Recalculates the relative abundance 
        psmelt() %>%  #Melts to long format
        arrange(desc(Abundance))

otu_ITSF3_long <- phyloseqITSF3 %>%  
        transform_sample_counts(function(x) {x * 100} ) %>% #Recalculates the relative abundance 
        psmelt() %>%  #Melts to long format
        arrange(desc(Abundance))


#Make vector with OTU names for differential abundant OTUs with effect size over 1 for each pair of comparisons
otu_16sF1_effect<-Aldex_16sF1.eff.table.all$OTU
otu_16sF2_effect<-Aldex_16sF2.eff.table.all$OTU
otu_16sF3_effect<-Aldex_16sF3.eff.table.all$OTU

otu_ITSF1_effect<-Aldex_ITSF1.eff.table.all$OTU
otu_ITSF2_effect<-Aldex_ITSF2.eff.table.all$OTU
otu_ITSF3_effect<-Aldex_ITSF3.eff.table.all$OTU


#Filter table to obtain only OTUs that with effect size greater than absolute 1
otu_16sF1_filter_effect <- filter(otu_16sF1_long, OTU %in% otu_16sF1_effect)
otu_16sF2_filter_effect <- filter(otu_16sF2_long, OTU %in% otu_16sF2_effect)
otu_16sF3_filter_effect <- filter(otu_16sF3_long, OTU %in% otu_16sF3_effect)

otu_ITSF1_filter_effect <- filter(otu_ITSF1_long, OTU %in% otu_ITSF1_effect)
otu_ITSF2_filter_effect <- filter(otu_ITSF2_long, OTU %in% otu_ITSF2_effect)
otu_ITSF3_filter_effect <- filter(otu_ITSF3_long, OTU %in% otu_ITSF3_effect)

otu_ITSF1_filter_effect$Family<-gsub("k__","", otu_ITSF1_filter_effect$Family)
otu_ITSF1_filter_effect$Family<-gsub("p__","", otu_ITSF1_filter_effect$Family)
otu_ITSF1_filter_effect$Family<-gsub("o__","", otu_ITSF1_filter_effect$Family)
otu_ITSF1_filter_effect$Family<-gsub("c__","", otu_ITSF1_filter_effect$Family)
otu_ITSF1_filter_effect$Family<-gsub("f__","", otu_ITSF1_filter_effect$Family)

otu_ITSF2_filter_effect$Family<-gsub("k__","", otu_ITSF2_filter_effect$Family)
otu_ITSF2_filter_effect$Family<-gsub("p__","", otu_ITSF2_filter_effect$Family)
otu_ITSF2_filter_effect$Family<-gsub("o__","", otu_ITSF2_filter_effect$Family)
otu_ITSF2_filter_effect$Family<-gsub("c__","", otu_ITSF2_filter_effect$Family)
otu_ITSF2_filter_effect$Family<-gsub("f__","", otu_ITSF2_filter_effect$Family)

otu_ITSF3_filter_effect$Family<-gsub("k__","", otu_ITSF3_filter_effect$Family)
otu_ITSF3_filter_effect$Family<-gsub("p__","", otu_ITSF3_filter_effect$Family)
otu_ITSF3_filter_effect$Family<-gsub("o__","", otu_ITSF3_filter_effect$Family)
otu_ITSF3_filter_effect$Family<-gsub("c__","", otu_ITSF3_filter_effect$Family)
otu_ITSF3_filter_effect$Family<-gsub("f__","", otu_ITSF3_filter_effect$Family)

#calculate statistics by OTU and Year
#16s
stat_16sF1<-otu_16sF1_filter_effect %>%
        group_by(OTU, Year, Family)%>%
        summarise(Mean=mean(Abundance, na.rm=TRUE))

stat_16sF2<-otu_16sF2_filter_effect %>%
        group_by(OTU, Year, Family)%>%
        summarise(Mean=mean(Abundance, na.rm=TRUE))

stat_16sF3<-otu_16sF3_filter_effect %>%
        group_by(OTU, Year, Family)%>%
        summarise(Mean=mean(Abundance, na.rm=TRUE))

#ITS
stat_ITSF1<-otu_ITSF1_filter_effect %>%
        group_by(OTU, Year, Family)%>%
        summarise(Mean=mean(Abundance, na.rm=TRUE))

stat_ITSF2<-otu_ITSF2_filter_effect %>%
        group_by(OTU, Year, Family)%>%
        summarise(Mean=mean(Abundance, na.rm=TRUE))

stat_ITSF3<-otu_ITSF3_filter_effect %>%
        group_by(OTU, Year, Family)%>%
        summarise(Mean=mean(Abundance, na.rm=TRUE))


#Reshape table to wide format
stat_16sF1_mean<-dcast(stat_16sF1, OTU + Family~ Year , value.var='Mean')
stat_16sF2_mean<-dcast(stat_16sF2, OTU + Family~ Year , value.var='Mean')
stat_16sF3_mean<-dcast(stat_16sF3, OTU + Family~ Year , value.var='Mean')

stat_ITSF1_mean<-dcast(stat_ITSF1, OTU + Family~ Year , value.var='Mean')
stat_ITSF2_mean<-dcast(stat_ITSF2, OTU + Family~ Year , value.var='Mean')
stat_ITSF3_mean<-dcast(stat_ITSF3, OTU + Family~ Year , value.var='Mean')

#Calculate log fold chage - log2 of the ratio of mean RA of first facility by second facility in comparison
stat_16sF1_mean$logchange<-log2(stat_16sF1_mean$Y1/stat_16sF1_mean$Y2)
stat_16sF2_mean$logchange<-log2(stat_16sF2_mean$Y1/stat_16sF2_mean$Y2)
stat_16sF3_mean$logchange<-log2(stat_16sF3_mean$Y1/stat_16sF3_mean$Y2)

stat_ITSF1_mean$logchange<-log2(stat_ITSF1_mean$Y1/stat_ITSF1_mean$Y2)
stat_ITSF2_mean$logchange<-log2(stat_ITSF2_mean$Y1/stat_ITSF2_mean$Y2)
stat_ITSF3_mean$logchange<-log2(stat_ITSF3_mean$Y1/stat_ITSF3_mean$Y2)



#Add facility identifier - Year with higher RA 
stat_16sF1_mean$Year<-ifelse(stat_16sF1_mean$logchange <0 , "Y2", "Y1")
stat_16sF2_mean$Year<-ifelse(stat_16sF2_mean$logchange <0 , "Y2", "Y1")
stat_16sF3_mean$Year<-ifelse(stat_16sF3_mean$logchange <0 , "Y2", "Y1")

stat_ITSF1_mean$Year<-ifelse(stat_ITSF1_mean$logchange <0 , "Y2", "Y1")
stat_ITSF2_mean$Year<-ifelse(stat_ITSF2_mean$logchange <0 , "Y2", "Y1")
stat_ITSF3_mean$Year<-ifelse(stat_ITSF3_mean$logchange <0 , "Y2", "Y1")


#Plot 
#Fig 4
#16s
Logfold_16sF1<-ggplot(stat_16sF1_mean, aes(x=logchange, y=reorder(OTU,logchange), fill=Year))+
        geom_bar(stat='identity', color='black')+
        geom_text(aes(label=Family), position = position_stack(vjust=0.5), size=2)+
        theme(axis.title.y = element_blank(), axis.ticks=element_line(color='black'),
              axis.text=element_text(size=8, color='black')) + 
        guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
        xlab("Log fold change (log2 Y1/Y2)") + 
        scale_x_continuous(limits=c(-9,8), breaks = c(-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8))+
        theme(panel.background = element_rect(fill=NA, color =NA),
              plot.background = element_rect(fill="transparent", color =NA),
              panel.border = element_rect(color="black", fill=NA)) +
        theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
        scale_fill_manual(values=c("#8D6CA480"))+
        theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))+
        geom_text(aes(label=Family), position = position_stack(vjust=0.5), size=2)+
        ggtitle("Bacteria F1 - Y1 v Y2", subtitle = "Differential abundant OTUs")

Logfold_16sF2<-ggplot(stat_16sF2_mean, aes(x=logchange, y=reorder(OTU,logchange), fill=Year))+
        geom_bar(stat='identity', color='black')+
        geom_text(aes(label=Family), position = position_stack(vjust=0.5), size=2)+
        theme(axis.title.y = element_blank(), axis.ticks=element_line(color='black'),
              axis.text=element_text(size=8, color='black')) + 
        guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
        xlab("Log fold change (log2 Y1/Y2)") + 
        scale_x_continuous(limits=c(-9,8), breaks = c(-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8))+
        theme(panel.background = element_rect(fill=NA, color =NA),
              plot.background = element_rect(fill="transparent", color =NA),
              panel.border = element_rect(color="black", fill=NA)) +
        theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
        scale_fill_manual(values=c("#D6879880"))+
        theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))+
        geom_text(aes(label=Family), position = position_stack(vjust=0.5), size=2)+
        ggtitle("Bacteria F2 - Y1 v Y2", subtitle = "Differential abundant OTUs")

Logfold_16sF3<-ggplot(stat_16sF3_mean, aes(x=logchange, y=reorder(OTU,logchange), fill=Year))+
        geom_bar(stat='identity', color='black')+
        geom_text(aes(label=Family), position = position_stack(vjust=0.5), size=2)+
        theme(axis.title.y = element_blank(), axis.ticks=element_line(color='black'),
              axis.text=element_text(size=8, color='black')) + 
        guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
        xlab("Log fold change (log2 Y1/Y2)") + 
        scale_x_continuous(limits=c(-9,8), breaks = c(-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8))+
        theme(panel.background = element_rect(fill=NA, color =NA),
              plot.background = element_rect(fill="transparent", color =NA),
              panel.border = element_rect(color="black", fill=NA)) +
        theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
        scale_fill_manual(values=c("#FDC96C80"))+
        theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))+
        geom_text(aes(label=Family), position = position_stack(vjust=0.5), size=2)+
        ggtitle("Bacteria F3 - Y1 v Y2", subtitle = "Differential abundant OTUs")

#ITS
Logfold_ITSF1<-ggplot(stat_ITSF1_mean, aes(x=logchange, y=reorder(OTU,logchange), fill=Year))+
        geom_bar(stat='identity', color='black')+
        geom_text(aes(label=Family), position = position_stack(vjust=0.5), size=2)+
        theme(axis.title.y = element_blank(), axis.ticks=element_line(color='black'),
              axis.text=element_text(size=8, color='black')) + 
        guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
        xlab("Log fold change (log2 Y1/Y2)") + 
        scale_x_continuous(limits=c(-9,8), breaks = c(-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8))+
        theme(panel.background = element_rect(fill=NA, color =NA),
              plot.background = element_rect(fill="transparent", color =NA),
              panel.border = element_rect(color="black", fill=NA)) +
        theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
        scale_fill_manual(values=c("#8D6CA480"))+
        theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))+
        geom_text(aes(label=Family), position = position_stack(vjust=0.5), size=2)+
        ggtitle("Fungi F1 - Y1 v Y2", subtitle = "Differential abundant OTUs")

Logfold_ITSF2<-ggplot(stat_ITSF2_mean, aes(x=logchange, y=reorder(OTU,logchange), fill=Year))+
        geom_bar(stat='identity', color='black')+
        geom_text(aes(label=Family), position = position_stack(vjust=0.5), size=2)+
        theme(axis.title.y = element_blank(), axis.ticks=element_line(color='black'),
              axis.text=element_text(size=8, color='black')) + 
        guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
        xlab("Log fold change (log2 Y1/Y2)") + 
        scale_x_continuous(limits=c(-9,8), breaks = c(-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8))+
        theme(panel.background = element_rect(fill=NA, color =NA),
              plot.background = element_rect(fill="transparent", color =NA),
              panel.border = element_rect(color="black", fill=NA)) +
        theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
        scale_fill_manual(values=c("#D6879880"))+
        theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))+
        geom_text(aes(label=Family), position = position_stack(vjust=0.5), size=2)+
        ggtitle("Fungi F2 - Y1 v Y2", subtitle = "Differential abundant OTUs")

Logfold_ITSF3<-ggplot(stat_ITSF3_mean, aes(x=logchange, y=reorder(OTU,logchange), fill=Year))+
        geom_bar(stat='identity', color='black')+
        geom_text(aes(label=Family), position = position_stack(vjust=0.5), size=2)+
        theme(axis.title.y = element_blank(), axis.ticks=element_line(color='black'),
              axis.text=element_text(size=8, color='black')) + 
        guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
        xlab("Log fold change (log2 Y1/Y2)") + 
        scale_x_continuous(limits=c(-9,8), breaks = c(-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8))+
        theme(panel.background = element_rect(fill=NA, color =NA),
              plot.background = element_rect(fill="transparent", color =NA),
              panel.border = element_rect(color="black", fill=NA)) +
        theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
        scale_fill_manual(values=c("#FDC96C80"))+
        theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))+
        geom_text(aes(label=Family), position = position_stack(vjust=0.5), size=2)+
        ggtitle("Fungi F3 - Y1 v Y2", subtitle = "Differential abundant OTUs")



#Combine plots
Logfold_byYear <- plot_grid(Logfold_16sF1, Logfold_16sF2, Logfold_16sF3, Logfold_ITSF1, Logfold_ITSF2, Logfold_ITSF3,
                        ncol=3, nrow=2, labels = c("A","B","C","D","E","F"), label_size = 20, vjust = 2, hjust = -1.5)
Logfold_byYear
ggsave("Aldex_LogFold_byYear.png", plot=Logfold_byYear, device="png", width=13, height=11, units="in", dpi=600)
ggsave("Aldex_LogFold_byYear.svg", plot=Logfold_byYear, device="svg", width=13, height=11, units="in", dpi=600)


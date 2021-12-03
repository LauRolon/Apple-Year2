#Seasonal analysis of two-year sampling of tree fruit packing environments
#Analysis of Nanopore data obtained by Narjol 

#Last updated MLR 09/15/2021

#Set working directory
setwd("/storage/work/m/mlr355/Apple/Nanopore")

#Attach libraries
library(taxize)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

#### Import data from Nanopore WIMP ####
M1<-read.csv('Taxonomy_M1.csv', header = TRUE, stringsAsFactors = FALSE)
M4<-read.csv('Taxonomy_M4.csv', header = TRUE, stringsAsFactors = FALSE)
M7<-read.csv('Taxonomy_M7.csv', header = TRUE, stringsAsFactors = FALSE)

M1<-as.data.frame(M1)
M4<-as.data.frame(M4)
M7<-as.data.frame(M7)

#Calculate number of taxonomonic ranks
M1$taxranks<-str_count(M1$lineage, pattern = ":")
M4$taxranks<-str_count(M4$lineage, pattern = ":")
M7$taxranks<-str_count(M7$lineage, pattern = ":")

max(M1$taxranks)
max(M4$taxranks)
max(M7$taxranks)

#Segregate lineage column into taxonomy ranks
Taxranks<-sprintf("Tax%s",seq(1:32)) #Column header for lineage

M1<-separate(M1, col=9, Taxranks, sep=":")
M1[Taxranks]<-lapply(M1[Taxranks], factor)
str(M1)

M4<-separate(M4, col=9, Taxranks, sep=":")
M4[taxranks]<-lapply(M4[Taxranks], factor)
str(M4)

M7<-separate(M7, col=9, Taxranks, sep=":")
M7[Taxranks]<-lapply(M7[Taxranks], factor)
str(M7)

#Subset by taxonomic rank
#M1
M1_classified<-subset(M1, exit_status=="Classified")
M1_unclassified<-subset(M1, exit_status=="Unclassified")
M1_virus<-subset(M1_classified, Tax2=="10239")
M1_cellorg<-subset(M1_classified, Tax2=="131567")
M1_bacteria<-subset(M1_cellorg, Tax3=="2")
M1_archaea<-subset(M1_cellorg, Tax3=="2157")
M1_eukaryota<-subset(M1_cellorg, Tax3=="2759")
M1_kingunclassified<-subset(M1_cellorg, is.na(M1_cellorg$Tax3))
M1_fungi<-subset(M1_eukaryota, Tax5=="4751")
M1_homosapiens<-subset(M1_eukaryota,taxID=="9606")
M1_listeria<-subset(M1_bacteria, Tax9=="1637")
M1_lmono<-subset(M1_bacteria, taxID=="1639")


M4_classified<-subset(M4, exit_status=="Classified")
M4_unclassified<-subset(M4, exit_status=="Unclassified")
M4_virus<-subset(M4_classified, Tax2=="10239")
M4_cellorg<-subset(M4_classified, Tax2=="131567")
M4_bacteria<-subset(M4_cellorg, Tax3=="2")
M4_archaea<-subset(M4_cellorg, Tax3=="2157")
M4_eukaryota<-subset(M4_cellorg, Tax3=="2759")
M4_kingunclassified<-subset(M4_cellorg, is.na(M1_cellorg$Tax3))
M4_fungi<-subset(M4_eukaryota, Tax5=="4751")
M4_homosapiens<-subset(M4_eukaryota,taxID=="9606")
M4_listeria<-subset(M4_bacteria, Tax9=="1637")
M4_lmono<-subset(M4_bacteria, taxID=="1639")


M7_classified<-subset(M7, exit_status=="Classified")
M7_unclassified<-subset(M7, exit_status=="Unclassified")
M7_virus<-subset(M7_classified, Tax2=="10239")
M7_cellorg<-subset(M7_classified, Tax2=="131567")
M7_bacteria<-subset(M7_cellorg, Tax3=="2")
M7_archaea<-subset(M7_cellorg, Tax3=="2157")
M7_eukaryota<-subset(M7_cellorg, Tax3=="2759")
M7_kingunclassified<-subset(M7_cellorg, is.na(M1_cellorg$Tax3))
M7_fungi<-subset(M7_eukaryota, Tax5=="4751")
M7_homosapiens<-subset(M7_eukaryota,taxID=="9606")
M7_listeria<-subset(M7_bacteria, Tax9=="1637")
M7_lmono<-subset(M7_bacteria, taxID=="1639")


####Bacteria####
#Calculate relative abundance by species
#M1
M1_bacteria_abundance<-subset(M1_bacteria[6:17]) %>% #Select columns with lineage and taxonomy information
  group_by(taxID, name) %>% #Remove duplicates
  tally() #Count reads in each group

sum(M1_bacteria_abundance$n) #Verify that all reads were counted

#Function to calculate relative abundance in percentage 
RA_M1<-(function(x) {x*100/(sum(M1_bacteria_abundance$n))}) 

#Calculate relative abundace
M1_bacteria_abundance$RA<-RA_M1(M1_bacteria_abundance$n)

#Exclude samples with less than 1%
M1_bacteria_abundance1<-subset(M1_bacteria_abundance, RA>=1)
M1_bacteria_abundance1$name<-as.character(M1_bacteria_abundance1$name)

#Get full taxonomy from NCBI - Adapted from Copyright 2021 Simone Maestri.
taxonomy_M1<-c()
raw_classification_M1 <- classification(M1_bacteria_abundance1$taxID, db = 'ncbi')
for (i in 1:length(M1_bacteria_abundance1$taxID)) {
  if (!is.null(ncol(M1_bacteria_abundance1$taxID[[i]]))) {
    if (nrow(M1_bacteria_abundance1$taxID[[i]]) > 1) {
      subspecies_name_curr <- M1_bacteria_abundance1$taxID[[i]][which(M1_bacteria_abundance1$taxID[[i]][, 2] == "subspecies"), 1]
      species_name_curr <- M1_bacteria_abundance1$taxID[[i]][which(M1_bacteria_abundance1$taxID[[i]][, 2] == "species"), 1]
      genus_name_curr <- M1_bacteria_abundance1$taxID[[i]][which(M1_bacteria_abundance1$taxID[[i]][, 2] == "genus"), 1]
      family_name_curr <- M1_bacteria_abundance1$taxID[[i]][which(M1_bacteria_abundance1$taxID[[i]][, 2] == "family"), 1]
      order_name_curr <- M1_bacteria_abundance1$taxID[[i]][which(M1_bacteria_abundance1$taxID[[i]][, 2] == "order"), 1]
      class_name_curr <- M1_bacteria_abundance1$taxID[[i]][which(M1_bacteria_abundance1$taxID[[i]][, 2] == "class"), 1]
      phylum_name_curr <- M1_bacteria_abundance1$taxID[[i]][which(M1_bacteria_abundance1$taxID[[i]][, 2] == "phylum"), 1]
      kingdom_name_curr <- M1_bacteria_abundance1$taxID[[i]][which(M1_bacteria_abundance1$taxID[[i]][, 2] == "superkingdom"), 1]
      taxonomy_M1[i] <- paste(kingdom_name_curr, phylum_name_curr, class_name_curr, order_name_curr, family_name_curr, genus_name_curr, species_name_curr, subspecies_name_curr, sep = "\t")
    } else {
      taxonomy_M1[i] <- "Unclassified"
    }
  } else {
    taxonomy_M1[i] <- "Unclassified"
  }
}

M1_bacteria_taxfam1<-tax_name(M1_bacteria_abundance1$name[1:10], get=c('Superkingdom','phylum','class','order','family','genus','species'), db='ncbi')
M1_bacteria_taxfam2<-tax_name(M1_bacteria_abundance1$name[11:17], get=c('Superkingdom','phylum','class','order','family','genus','species'), db='ncbi')

M1_bacteria_tax<-rbind(M1_bacteria_taxfam1,M1_bacteria_taxfam2)
M1_bacteria_tax_RA<-bind_cols(M1_bacteria_tax,M1_bacteria_abundance1)

write.csv(M1_bacteria_tax_RA,"M1_species_RA.csv")

#M4
M4_bacteria_abundance<-subset(M4_bacteria[6:17]) %>% #Select columns with lineage and taxonomy information
  group_by(taxID, name) %>% #Remove duplicates
  tally() #Count reads in each group

sum(M4_bacteria_abundance$n) #Verify that all reads were counted

RA_M4<-(function(x) {x/(sum(M4_bacteria_abundance$n))}) #Function to calculate realtive abundace

M4_bacteria_abundance$RA<-RA_M4(M4_bacteria_abundance$n) #Calculate relative abundance per group

M4_bacteria_abundance1<-subset(M4_bacteria_abundance, RA>=0.01) #Select all rows were RA>=0.01

sum(M4_bacteria_abundance$n)-sum(M4_bacteria_abundance1$n) #Number of reads below 1% RA (For "Others" category)
1-(sum(M4_bacteria_abundance1$RA)) # RA for "Others" category

M4_bacteria_abundance1$name<-as.character(M4_bacteria_abundance1$name)

M4_bacteria_taxfam1<-tax_name(M4_bacteria_abundance1$name[1:9], get=c('Superkingdom','phylum','class','order','family','genus','species'), db='ncbi') #Can only query NCBI for 10 results per time >:(
M4_bacteria_taxfam2<-tax_name(M4_bacteria_abundance1$name[10:16], get=c('Superkingdom','phylum','class','order','family','genus','species'), db='ncbi')

M4_bacteria_tax<-rbind(M4_bacteria_taxfam1,M4_bacteria_taxfam2)
M4_bacteria_tax_RA<-bind_cols(M4_bacteria_tax,M4_bacteria_abundance1)

write.csv(M4_bacteria_tax_RA,"M4_species_RA.csv")

#M7
M7_bacteria_abundance<-subset(M7_bacteria[6:17]) %>% #Select columns with lineage and taxonomy information
  group_by(taxID, name) %>% #Remove duplicates
  tally() #Count reads in each group

sum(M7_bacteria_abundance$n) #Verify that all reads were counted

RA_M7<-(function(x) {x/(sum(M7_bacteria_abundance$n))}) #Function to calculate realtive abundace

M7_bacteria_abundance$RA<-RA_M7(M7_bacteria_abundance$n) #Calculate relative abundance per group

M7_bacteria_abundance1<-subset(M7_bacteria_abundance, RA>=0.01) #Select all rows were RA>=0.01

sum(M7_bacteria_abundance$n)-sum(M7_bacteria_abundance1$n) #Number of reads below 1% RA (For "Others" category)
1-(sum(M7_bacteria_abundance1$RA)) # RA for "Others" category

M7_bacteria_abundance1$name<-as.character(M7_bacteria_abundance1$name)

M7_bacteria_taxfam1<-tax_name(M7_bacteria_abundance1$name[1:10], get=c('Superkingdom','phylum','class','order','family','genus','species'), db='ncbi') #Can only query NCBI for 10 results per time >:(

M7_bacteria_tax_RA<-bind_cols(M7_bacteria_taxfam1,M7_bacteria_abundance1)

write.csv(M7_bacteria_tax_RA,"M7_species_RA.csv")

#Open .csv file in Excel. Add "Others" row with n and RA calculated above. 
#Add metadata information for the sample
#In cells were there is no taxonomic information, copy from prior tax rank.

M1_bacteria_species_RA<-read.csv("M1_species_RA.csv", header = TRUE, row.names=1)
M4_bacteria_species_RA<-read.csv("M4_species_RA.csv", header = TRUE, row.names=1)
M7_bacteria_species_RA<-read.csv("M7_species_RA.csv", header = TRUE, row.names=1)

Bacteria_species_all<-rbind(M1_bacteria_species_RA,M4_bacteria_species_RA,M7_bacteria_species_RA)
write.csv(Bacteria_species_all, "Bacteria_all.csv")

#Add number label to names

Bacteria_species<-read.csv("Bacteria_all.csv", header = TRUE, row.names = 1)

#Plot bacterial species by facility
Barplot_species<- ggplot(Bacteria_species, aes(x = Facility, y = RA , fill = L.name))  + 
  geom_bar(stat = "identity", color='black') + 
  theme(legend.text=element_text(size=7), legend.title= element_text(size=7)) +
  theme(axis.title.x = element_text(size=15), axis.ticks.x=element_blank(),
        axis.title.y=element_text(size=15)) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  ylab("Relative Abundance") + theme(panel.background = element_rect(fill="transparent", color =NA),
                                     plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_blank(), strip.text = element_text(size=15),
        panel.border = element_rect(color="black", fill=NA))+
  geom_text(aes(label=Label), position = position_stack(vjust=0.5), size=3)

ggsave("Barplot_bacteria_species.png", plot=Barplot_species, device="png", width=13, height=8, units="in", dpi=600)

#Calculate relative abundance by family
M1_bacteria_abundance_family<-subset(M1_bacteria[6:17]) %>% #Select columns with lineage and taxonomy information
  group_by(taxID, name, Kingdom, Phylum, Class, Order, Family) %>% #Remove duplicates
  tally() #Count reads in each group
M1_bacteria_abundance_family$RA<-RA_M1(M1_bacteria_abundance_family$n)

M4_bacteria_abundance_family<-subset(M4_bacteria[6:17]) %>% #Select columns with lineage and taxonomy information
  group_by(taxID, name, Kingdom, Phylum, Class, Order, Family) %>% #Remove duplicates
  tally() #Count reads in each group
M4_bacteria_abundance_family$RA<-RA_M4(M4_bacteria_abundance_family$n)

M7_bacteria_abundance_family<-subset(M7_bacteria[6:17]) %>% #Select columns with lineage and taxonomy information
  group_by(taxID, name, Kingdom, Phylum, Class, Order, Family) %>% #Remove duplicates
  tally() #Count reads in each group
M7_bacteria_abundance_family$RA<-RA_M7(M7_bacteria_abundance_family$n)

write.csv(M1_bacteria_abundance_family,"M1_family_RA.csv")
write.csv(M4_bacteria_abundance_family,"M4_family_RA.csv")
write.csv(M7_bacteria_abundance_family,"M7_family_RA.csv")

#Open in Excel and replace NA with prior tax rank
M1_bacteria_family<-read.csv("M1_family_RA.csv", header = TRUE,row.names = 1)

M1_bacteria_RA_family<-M1_bacteria_family %>% #Select columns with lineage and taxonomy information
  group_by(Family) %>% #Remove duplicates
  tally(RA) #Count RA in each group

M1_bacteria_family1<-subset(M1_bacteria_RA_family, n>=0.01) #Select all rows were RA>=0.01

M4_bacteria_family<-read.csv("M4_family_RA.csv", header = TRUE,row.names = 1)

M4_bacteria_RA_family<-M4_bacteria_family %>% #Select columns with lineage and taxonomy information
  group_by(Family) %>% #Remove duplicates
  tally(RA) #Count RA in each group

M4_bacteria_family1<-subset(M4_bacteria_RA_family, n>=0.01) #Select all rows were RA>=0.01

M7_bacteria_family<-read.csv("M7_family_RA.csv", header = TRUE,row.names = 1)

M7_bacteria_RA_family<-M7_bacteria_family %>% #Select columns with lineage and taxonomy information
  group_by(Family) %>% #Remove duplicates
  tally(RA) #Count RA in each group

M7_bacteria_family1<-subset(M7_bacteria_RA_family, n>=0.01) #Select all rows were RA>=0.01

write.csv(M1_bacteria_family1, "M1_family.csv")
write.csv(M4_bacteria_family1, "M4_family.csv")
write.csv(M7_bacteria_family1, "M7_family.csv")

#Open files in Excel. Combine all tables and add metadata information.
#Look on NCBI the taxonomy names for each taxID
#Add column for type of sequencing, sample names, and labels.
#Plot families by type of sequencing

Bacteria_family<-read.csv("Bacteria_family_all_RLE.csv", header = TRUE)

Barplot_family<- ggplot(Bacteria_family, aes(x = Seq, y = RA , fill = L.family))  + 
  facet_grid(.~Facility)+
  geom_bar(stat = "identity", color='black') + 
  theme(legend.text=element_text(size=10), legend.title= element_blank()) +
  theme(axis.title.y = element_text(size=20,color = 'black'),axis.title.x = element_blank(),
        axis.ticks.x=element_blank(), axis.text = element_text(size=15, color='black'), axis.text.x = element_text(angle=90)) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  ylab("Relative Abundance") + theme(panel.background = element_rect(fill="transparent", color =NA),
                                     plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_blank(), strip.text = element_text(size=18),
        panel.border = element_rect(color="black", fill=NA))+
  geom_text(aes(label=Label), position = position_stack(vjust=0.5), size=3, color='black')+
  theme(panel.background = element_rect(fill='transparent', color = NA),
        plot.background = element_rect(fill = 'transparent',color = NA), 
        panel.border = element_rect(color='black',fill = NA,size=1))+
  theme(strip.background = element_rect(color='black', fill='white', size=1, linetype = 'solid'))+
  scale_fill_viridis_d(option='inferno', begin=0.2, end=0.8, alpha = 0.85)

ggsave("Barplot_bacteria_family_poster.png", plot=Barplot_family, device="png", width=13, height=10, units="in", dpi=600)


####Fungi####
#Calculate relative abundace by species
#M1
M1_fungi_abundance<-subset(M1_fungi[6:21]) %>% #Select columns with lineage and taxonomy information
  group_by(taxID, name) %>% #Remove duplicates
  tally() #Count reads in each group

sum(M1_fungi_abundance$n) #Verify that all reads were counted

RA_fungi_M1<-(function(x) {x/(sum(M1_fungi_abundance$n))}) #Function to calculate realtive abundace

sum(M1_fungi_abundance$n)-sum(M1_fungi_abundance1$n) #Number of reads below 1% RA (For "Others" category)
1-(sum(M1_fungi_abundance1$RA)) # RA for "Others" category

M1_fungi_abundance$RA<-RA_fungi_M1(M1_fungi_abundance$n) #Calculate relative abundance per group
M1_fungi_abundance1<-subset(M1_fungi_abundance, RA>=0.01)

M1_fungi_taxfam1<-tax_name(M1_fungi_abundance1$name, get=c('Superkingdom','phylum','class','order','family','genus','species'), db='ncbi')

M1_fungi_tax_RA<-bind_cols(M1_fungi_taxfam1,M1_fungi_abundance1)

write.csv(M1_fungi_tax_RA,"M1_fungi_species_RA.csv")



#M4
M4_fungi_abundance<-subset(M4_fungi[6:21]) %>% #Select columns with lineage and taxonomy information
  group_by(taxID, name) %>% #Remove duplicates
  tally() #Count reads in each group

sum(M4_fungi_abundance$n) #Verify that all reads were counted

RA_fungi_M4<-(function(x) {x/(sum(M4_fungi_abundance$n))}) #Function to calculate realtive abundace

sum(M4_fungi_abundance$n)-sum(M4_fungi_abundance1$n) #Number of reads below 1% RA (For "Others" category)
1-(sum(M4_fungi_abundance1$RA)) # RA for "Others" category

M4_fungi_abundance$RA<-RA_fungi_M4(M4_fungi_abundance$n) #Calculate relative abundance per group

M4_fungi_abundance1<-subset(M4_fungi_abundance, RA>=0.01)

write.csv(M4_fungi_abundance1,"M4_fungi_species_RA.csv")

#M7
M7_fungi_abundance<-subset(M7_fungi[6:21]) %>% #Select columns with lineage and taxonomy information
  group_by(taxID, name) %>% #Remove duplicates
  tally() #Count reads in each group

sum(M7_fungi_abundance$n) #Verify that all reads were counted

RA_fungi_M7<-(function(x) {x/(sum(M7_fungi_abundance$n))}) #Function to calculate realtive abundace

sum(M7_fungi_abundance$n)-sum(M7_fungi_abundance1$n) #Number of reads below 1% RA (For "Others" category)
1-(sum(M7_fungi_abundance1$RA)) # RA for "Others" category

M7_fungi_abundance$RA<-RA_fungi_M7(M7_fungi_abundance$n) #Calculate relative abundance per group

M7_fungi_abundance1<-subset(M7_fungi_abundance, RA>=0.01)

write.csv(M7_fungi_abundance1,"M7_fungi_species_RA.csv")

#Open .csv file in Excel. Add "Others" row with n and RA calculated above. 
#Add metadata information for the sample
#In cells were there is no taxonomic information, copy from prior tax rank.

fungi_species<-read.csv("Fungi_all.csv", header = TRUE, row.names = 1)

#Plot fungal species by facility
Barplot_fungi_species<- ggplot(fungi_species, aes(x = Facility, y = RA , fill = L.Species))  + 
  geom_bar(stat = "identity", color='black') + 
  theme(legend.text=element_text(size=7), legend.title= element_text(size=7)) +
  theme(axis.title.x = element_text(size=15), axis.ticks.x=element_blank(),
        axis.title.y=element_text(size=15)) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  ylab("Relative Abundance") + theme(panel.background = element_rect(fill="transparent", color =NA),
                                     plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_blank(), strip.text = element_text(size=15),
        panel.border = element_rect(color="black", fill=NA))+
  geom_text(aes(label=Label), position = position_stack(vjust=0.5), size=3)

ggsave("Barplot_fungi_species.png", plot=Barplot_fungi_species, device="png", width=13, height=8, units="in", dpi=600)

#Calculate relative abundance by family
M1_fungi_abundance_family<-subset(M1_fungi[6:21]) %>% #Select columns with lineage and taxonomy information
  group_by(taxID, name, Tax5, Tax6, Tax7, Tax8, Tax9, Tax10, Tax11) %>% #Remove duplicates
  tally() #Count reads in each group

M1_fungi_abundance_family$RA<-RA_fungi_M1(M1_fungi_abundance_family$n)

M4_fungi_abundance_family<-subset(M4_fungi[6:22]) %>% #Select columns with lineage and taxonomy information
  group_by(taxID, name, Tax5, Tax6, Tax7, Tax8, Tax9, Tax10, Tax11) %>% #Remove duplicates
  tally() #Count reads in each group
M4_fungi_abundance_family$RA<-RA_fungi_M4(M4_fungi_abundance_family$n)

M7_fungi_abundance_family<-subset(M7_fungi[6:22]) %>% #Select columns with lineage and taxonomy information
  group_by(taxID, name, Tax5, Tax6, Tax7, Tax8, Tax9, Tax10, Tax11) %>% #Remove duplicates
  tally() #Count reads in each group
M7_fungi_abundance_family$RA<-RA_fungi_M7(M7_fungi_abundance_family$n)

write.csv(M1_fungi_abundance_family,"M1_fungi_family_RA.csv")
write.csv(M4_fungi_abundance_family,"M4_fungi_family_RA.csv")
write.csv(M7_fungi_abundance_family,"M7_fungi_family_RA.csv")

#Open in Excel and replace NA with prior tax rank
M1_fungi_family<-read.csv("M1_fungi_family_RA.csv", header = TRUE,row.names = 1)

M1_fungi_RA_family<-M1_fungi_family %>% #Select columns with lineage and taxonomy information
  group_by(Tax11) %>% #Remove duplicates
  tally(RA) #Count RA in each group

M1_fungi_family1<-subset(M1_fungi_RA_family, n>=0.01) #Select all rows were RA>=0.01

M4_fungi_family<-read.csv("M4_fungi_family_RA.csv", header = TRUE,row.names = 1)

M4_fungi_RA_family<-M4_fungi_family %>% #Select columns with lineage and taxonomy information
  group_by(Tax11) %>% #Remove duplicates
  tally(RA) #Count RA in each group

M4_fungi_family1<-subset(M4_fungi_RA_family, n>=0.01) #Select all rows were RA>=0.01

M7_fungi_family<-read.csv("M7_fungi_family_RA.csv", header = TRUE,row.names = 1)

M7_fungi_RA_family<-M7_fungi_family %>% #Select columns with lineage and taxonomy information
  group_by(Tax11) %>% #Remove duplicates
  tally(RA) #Count RA in each group

M7_fungi_family1<-subset(M7_fungi_RA_family, n>=0.01) #Select all rows were RA>=0.01

write.csv(M1_fungi_family1, "M1_fungi_family.csv")
write.csv(M4_fungi_family1, "M4_fungi_family.csv")
write.csv(M7_fungi_family1, "M7_fungi_family.csv")

#Open files in Excel. Combine all tables and add metadata information.
#Look on NCBI the taxonomy names for each taxID
#Add column for type of sequencing, sample names, and labels.
#Plot families by type of sequencing

Fungi_family<-read.csv("All_fungi_family.csv", header = TRUE)

Barplot_fungi_family<- ggplot(Fungi_family, aes(x = Seq, y = RA , fill = L.Family))  + 
  facet_grid(Facility~.)+
  geom_bar(stat = "identity", color='black') + 
  theme(legend.text=element_text(size=9), legend.title= element_text(size=9)) +
  theme(axis.title.x = element_text(size=15), axis.ticks.x=element_blank(),
        axis.title.y=element_text(size=15)) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  ylab("Relative Abundance") + theme(panel.background = element_rect(fill="transparent", color =NA),
                                     plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_blank(), strip.text = element_text(size=15),
        panel.border = element_rect(color="black", fill=NA))+
  geom_text(aes(label=Label), position = position_stack(vjust=0.5), size=3)

ggsave("Barplot_fungi_family.png", plot=Barplot_fungi_family, device="png", width=13, height=10, units="in", dpi=600)



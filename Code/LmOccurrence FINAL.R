#Two-year monitoring of Lm in apple packing houses
#Analysis of Lm occurrence
#Laura Rolon
#Last updated: 05/26/21

#Set working directory
setwd()

#### FIG 2. Occurrence of Lm in tree fruit packing houses ####
#Attach libraries
library(readxl)
library(ggplot2)
library(cowplot)
library(LaCroixColoR)
library(svglite)

#Load the file containing data on L.monocytogenes occurrence
lm_facility <- read_excel('Lm_occurance.xlsx', sheet=1, col_names = TRUE) #by facility
lm_section <- read_excel('Lm_occurance.xlsx', sheet=2, col_names = TRUE) #by section
#lm_facsec <- read_excel('Lm_occurance.xlsx', sheet=3, col_names = TRUE) #by section and facility

# Panel A - Barplot of Lm occurrance by facility per year
lmoccurrence_facility <- ggplot(data = lm_facility, aes(x=Year, y=Percentage, fill=L.monocytogenes))+
  geom_bar(stat = "identity",width = 0.8)+
  facet_grid(.~Facility) +
  geom_text(aes(label=Percentage), position = position_stack(vjust=0.5), size=4) +
  theme(panel.background = element_rect(fill='transparent', color = NA),
        plot.background = element_rect(fill = 'transparent',color = NA), 
        panel.border = element_rect(color='black',fill = NA, size=0.8))+
  theme(strip.text = element_text(size=11, color='black'))+
  theme(legend.text=element_text(size=11), legend.title= element_text(size=11, face = 'italic')) +
  theme(axis.text.x = element_text(size=13, color='black')) +
  theme(axis.title = element_text(size=13, color = 'black')) +
  theme(strip.background = element_rect(color='black', fill=NA, size=0.8, linetype = 'solid'))+
  scale_x_discrete(name = "") +
  scale_y_discrete(name = "Percentage (%)")+
  scale_fill_manual(values = c('#F8CD9C','#EA7580'))+
  theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(legend.position = 'bottom')
lmoccurrence_facility
ggsave("lm_facility_Y1Y2percentage.png", plot = lmoccurrence_facility, device="png", width=8, height=5, units="in",dpi=600)


#Panel B - Barplotof Lm occurrance by section per year
lmoccurrence_section <- ggplot(data = lm_section, aes(x=Year, y=Percentage, fill=L.monocytogenes))+
  geom_bar(stat = "identity",width = 0.8)+
  facet_grid(.~reorder(Section,Order)) +
  scale_x_discrete(name = "") +
  scale_y_discrete(name = "Percentage (%)")+
  geom_text(aes(label=Percentage), position = position_stack(vjust=0.5), size=4) +
  theme(strip.text = element_text(size=11, color='black'))+
  theme(panel.background = element_rect(fill='transparent', color = NA),
        plot.background = element_rect(fill = 'transparent',color = NA), 
        panel.border = element_rect(color='black',fill = NA, size = 0.8))+
  theme(legend.text=element_text(size=11, color = 'black'), legend.title= element_text(size=11, face = 'italic', color='black')) +
  theme(axis.text.x = element_text(size=13, color='black'), axis.text.y = element_text(size = 17,color='black')) +
  theme(axis.title = element_text(size=13, color = 'black')) +
  theme(legend.position = 'bottom')+
  theme(strip.background = element_rect(color='black', fill='white', size=1, linetype = 'solid'))+
  scale_fill_manual(values = c('#F8CD9C','#EA7580'))+
  theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))
lmoccurrence_section
ggsave("lm_section_Y1Y2percetage.png", plot = lmoccurrence_section, device="png", width=8, height=5, units="in",dpi=600)


#Combine panels
Lmplot_top = plot_grid(lmoccurrence_facility+ theme(legend.position = "none"),
                         lmoccurrence_section,
                         ncol=1, nrow=2, labels = c("A","B"), label_size = 20, rel_heights = c(0.85,1))
Lmplot_top

# Lmplot_all = plot_grid(Lmplot_top,lmoccurrence_facsec, ncol=1, nrow = 2, rel_heights = c(2,2.5), labels = c("","C"), 
#                        label_size = 20)
# Lmplot_all

ggsave('LmocurranceAB_percentage.png', plot=Lmplot_top, device='png', width = 8, height=8, units = 'in', dpi = 600)

ggsave('LmocurranceAB_percentage.eps', plot=Lmplot_top, device='eps', width = 8, height=8, units = 'in', dpi = 600)


#### Statistical analysis ####
#Attach libraries
library(MASS) #for Chi-quare and Likelihood test
library(fmsb) #for pairwise Fisher's test.

## Y1: By facility
#Make contingency table
FacY1 <- c(28,11,0,39,23,16)
FacY1.table <- as.table(matrix(FacY1, nrow = 2, byrow = FALSE, dimnames = list(Lm = c('Absent', 'Present'), 
                                                                               Facility = c('F1', 'F2','F3'))))
addmargins(FacY1.table)

#Chi-square table
loglm( ~ Lm + Facility, data = FacY1.table)

# Fisher's pairwise comparison of proportions with Bonferroni correction
pairwise.fisher.test(c(11,39,16), c(39,39,39), p.adjust.method = 'bonferroni') 

## Y1: By Section 
#Make contingency table
SecY1 <- c(18,21,12,27,21,18)
SecY1.table <- as.table(matrix(SecY1, nrow = 2, byrow = FALSE, dimnames = list(Lm = c('Absent', 'Present'), 
                                                                               Section = c('Dry', 'Wash','Wax'))))
addmargins(SecY1.table)

#Chi-square table
loglm( ~ Lm + Section, data = SecY1.table)


## Y2: By facility
#Make contingency table
LmY2 <- c(21,15,0,36,5,31)
LmY2.table <- as.table(matrix(LmY2, nrow = 2, byrow = FALSE, dimnames = list(Lm = c('Absent', 'Present'), 
                                                                             Facility = c('F1', 'F2','F3'))))
addmargins(LmY2.table)

#Chi-square table
loglm( ~ Lm + Facility, data = LmY2.table) #prints Chi-square table

# Fisher's pairwise comparison of proportions with Bonferroni correction
pairwise.fisher.test(c(15,36,31), c(36,36,36), p.adjust.method = 'bonferroni') 

##Y2: By Section 
#Make contingency table
SecY2 <- c(25,11,29,7,28,8)
SecY2.table <- as.table(matrix(SecY2, nrow = 2, byrow = FALSE, dimnames = list(Lm = c('Absent', 'Present'), 
                                                                               Section = c('Dry', 'Wash','Wax'))))
addmargins(SecY2.table)

#Chi-square table
loglm( ~ Lm + Section, data = SecY2.table)

## Two-year comparison: by facility
pos.F1<-c(11,15)
tot.F1<-c(39,36)
prop.test(pos.F1, tot.F1) #X-squared = 0.96239, df = 1, p-value = 0.3266

pos.F3<-c(16,31)
tot.F3<-c(39,36)
prop.test(pos.F3, tot.F3) #X-squared = 14.395, df = 1, p-value = 0.0001482

## Two-year comparison: by section
pos.wash<-c(27,29)
tot.wash<-c(39,36)
prop.test(pos.wash, tot.wash) #X-squared = 0.74115, df = 1, p-value = 0.3893

pos.dry<-c(21,25)
tot.dry<-c(39,36)
prop.test(pos.dry, tot.dry) #X-squared = 1.3191, df = 1, p-value = 0.2507

pos.wax<-c(18,28)
tot.wax<-c(39,36)
prop.test(pos.wax, tot.wax) #X-squared = 6.617, df = 1, p-value = 0.0101

##Two-year comparison: by facility and section
#F1-wash
pos.F1wash<-c(7,7)
tot.F1wash<-c(13,12)
prop.test(pos.F1wash, tot.F1wash) #X-squared = 5.388e-31, df = 1, p-value = 1

#F1-dry
pos.F1dry<-c(3,3)
tot.F1dry<-c(13,12)
prop.test(pos.F1dry, tot.F1dry) #X-squared = 1.9764e-31, df = 1, p-value = 1

#F1-wax
pos.F1wax<-c(1,5)
tot.F1wax<-c(13,12)
prop.test(pos.F1wax, tot.F1wax) #X-squared = 2.3058, df = 1, p-value = 0.1289

#F3-wash
pos.F3wash<-c(7,10)
tot.F3wash<-c(13,12)
prop.test(pos.F3wash, tot.F3wash) #X-squared = 1.3224, df = 1, p-value = 0.2502

#F3-dry
pos.F3dry<-c(5,10)
tot.F3dry<-c(13,12)
prop.test(pos.F3dry, tot.F3dry) #X-squared = 3.5323, df = 1, p-value = 0.06018

#F3-wax
pos.F3wax<-c(4,11)
tot.F3wax<-c(13,12)
prop.test(pos.F3wax, tot.F3wax) #X-squared = 7.2716, df = 1, p-value = 0.007005

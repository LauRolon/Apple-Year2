#Two-year monitoring of Lm in apple packing houses
#Analysis of read sequencing quality 
#Laura Rolon
#Last updated: 05/26/21

#Set working directory
setwd()

#Attach libraries
library(ggplot2)
library(viridis)
library(LaCroixColoR)
library(psych)
library(readxl)
library(cowplot)

#### Analysis of sequencing reads prior to bioinformatics ####
#Import data - number of reads after sequencing - 
#16s
SEQreads_16s<-read_excel('Seq_Reads.xlsx', col_names = TRUE, sheet=1)
str(SEQreads_16s)
max(SEQreads_16s$Readsx1000) # 567.957

#ITS
SEQreads_ITS<-read_excel('Seq_Reads.xlsx', col_names = TRUE, sheet=2)
str(SEQreads_ITS)
max(SEQreads_ITS$Readsx1000) # 400.204

#Summary statistics
describeBy(SEQreads_16s$Reads, group=SEQreads_16s$Year, mat = TRUE)
describeBy(SEQreads_ITS$Reads, group=SEQreads_ITS$Year, mat = TRUE)

#Statistical analysis
t.test_16s<-t.test(Reads~Year, data=SEQreads_16s) #t = 1.3722, df = 276.68, p-value = 0.1711
t.test_ITS<-t.test(Reads~Year, data=SEQreads_ITS) # t = -11.752, df = 424.3, p-value < 2.2e-16

#Plot reads by year
SEQreads16s_plot<-ggplot(SEQreads_16s, aes(x=Year, y=Readsx1000, fill=Year))+
  geom_boxplot(color='black')+
  scale_fill_manual(values = lacroix_palette('Pamplemousse', 2, type='discrete'))+
  theme(panel.background = element_rect(fill='transparent', color = NA),
        plot.background = element_rect(fill = 'transparent',color = NA), 
        panel.border = element_rect(color='black',fill = NA, size=0.8))+
  theme(legend.text=element_text(size=13), legend.title= element_blank()) +
  theme(axis.text.x = element_text(size=13, color='black'), axis.ticks.x = element_blank()) +
  theme(axis.ticks.y = element_line(color='black'), axis.text.y = element_text(color='black', size=11))+
  theme(axis.title = element_text(size=13, color = 'black')) +
  theme(strip.background = element_rect(color='black', fill=NA, size=0.8, linetype = 'solid'))+
  scale_y_continuous(name="Number of reads (x 1,000)",  breaks = seq(0,600, 50), limits = c(0,600))+
  theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))

SEQreadsITS_plot<-ggplot(SEQreads_ITS, aes(x=Year, y=Readsx1000, fill=Year))+
  geom_boxplot(color='black')+
  scale_fill_manual(values = lacroix_palette('Pamplemousse', 2, type='discrete'))+
  theme(panel.background = element_rect(fill='transparent', color = NA),
        plot.background = element_rect(fill = 'transparent',color = NA), 
        panel.border = element_rect(color='black',fill = NA, size=0.8))+
  theme(legend.text=element_text(size=13), legend.title= element_blank()) +
  theme(axis.text.x = element_text(size=13, color='black'), axis.ticks.x = element_blank()) +
  theme(axis.ticks.y = element_line(color='black'), axis.text.y = element_text(color='black', size=11))+
  theme(axis.title = element_text(size=13, color = 'black')) +
  theme(strip.background = element_rect(color='black', fill=NA, size=0.8, linetype = 'solid'))+
  scale_y_continuous(name="Number of reads (x 1,000)",  breaks = seq(0,600, 50), limits = c(0,600))+
  theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))




#### Analysis of OTU reads after to bioinformatics ####
OTUreads_16s<-read_excel('OTUreads.xlsx', col_names = TRUE, sheet=1)
OTUreads_ITS<-read_excel('OTUreads.xlsx', col_names = TRUE, sheet=2)

#Summary statistics
describeBy(OTUreads_16s$Reads, group=OTUreads_16s$Year, mat = TRUE)
describeBy(OTUreads_ITS$Reads, group=OTUreads_ITS$Year, mat = TRUE)

#Statistical analysis
t.test_otu16s<-t.test(Reads~Year, data=OTUreads_16s) #t = -0.96598, df = 126.6, p-value = 0.3359
t.test_otuITS<-t.test(Reads~Year, data=OTUreads_ITS) # t = -9.4707, df = 155.56, p-value < 2.2e-16

#Plot reads by year
OTUreads16s_plot<-ggplot(OTUreads_16s, aes(x=Year, y=Readsx1000, fill=Year))+
  geom_boxplot(color='black')+
  scale_fill_manual(values = lacroix_palette('Pamplemousse', 2, type='discrete'))+
  theme(panel.background = element_rect(fill='transparent', color = NA),
        plot.background = element_rect(fill = 'transparent',color = NA), 
        panel.border = element_rect(color='black',fill = NA, size=0.8))+
  theme(legend.text=element_text(size=13), legend.title= element_blank()) +
  theme(axis.text.x = element_text(size=13, color='black'), axis.ticks.x = element_blank()) +
  theme(axis.ticks.y = element_line(color='black'), axis.text.y = element_text(color='black', size=11))+
  theme(axis.title = element_text(size=13, color = 'black')) +
  theme(strip.background = element_rect(color='black', fill=NA, size=0.8, linetype = 'solid'))+
  scale_y_continuous(name="Number of reads (x 1,000)",  breaks = seq(0,500, 50), limits = c(0,500))+
  theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))

OTUreadsITS_plot<-ggplot(OTUreads_ITS, aes(x=Year, y=Readsx1000, fill=Year))+
  geom_boxplot(color='black')+
  scale_fill_manual(values = lacroix_palette('Pamplemousse', 2, type='discrete'))+
  theme(panel.background = element_rect(fill='transparent', color = NA),
        plot.background = element_rect(fill = 'transparent',color = NA), 
        panel.border = element_rect(color='black',fill = NA, size=0.8))+
  theme(legend.text=element_text(size=13), legend.title= element_blank()) +
  theme(axis.text.x = element_text(size=13, color='black'), axis.ticks.x = element_blank()) +
  theme(axis.ticks.y = element_line(color='black'), axis.text.y = element_text(color='black', size=11))+
  theme(axis.title = element_text(size=13, color = 'black')) +
  theme(strip.background = element_rect(color='black', fill=NA, size=0.8, linetype = 'solid'))+
  scale_y_continuous(name="Number of reads (x 1,000)",  breaks = seq(0,500, 50), limits = c(0,500))+
  theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))


#Combine plots
FigS1 = plot_grid(SEQreads16s_plot+ theme(legend.position = "none"),
                      SEQreadsITS_plot+ theme(legend.position = "none"),
                      OTUreads16s_plot+ theme(legend.position = "none"),
                      OTUreadsITS_plot+ theme(legend.position = "none"),
                      ncol=2, nrow=2, labels = c("A","B","C","D"), label_size = 20)
FigS1

ggsave('FigS1.png', plot=FigS1, device='png', width = 8, height=10, units = 'in', dpi = 600)

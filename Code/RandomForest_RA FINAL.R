#Random Forest Analysis of microbiome data

#Looking for ASVs that classify samples into Lm positive or negative
#Written by Taejung Chung
#Updated by Laura Rolon / Taejung Chung
#Last updated: 6/24/22

#Load packages
library(ape) 
library(phyloseq)
library(randomForest) 
library(caret) 
library(svglite)
library(zCompositions)
library(compositions)

#Set working directory to where files are located
setwd("/storage/work/m/mlr355/Apple/Downstream")

#### IMPORT DATA ####

##16s
asv_16s<-read.csv("ASV_16s_clean.csv", header=TRUE, row.names = 1)
taxon_16s<-as.data.frame(read.csv("Taxon_16s_clean.csv", header = TRUE, row.names = 1))
metadata.16s<-read.csv("metadata_16s_clean.csv", header = TRUE, row.names = 1)

##ITS
asv_ITS<-read.csv("ASV_ITS_clean.csv", header=TRUE, row.names = 1)
taxon_ITS<-as.data.frame(read.csv("Taxon_ITS_clean.csv", header = TRUE, row.names = 1))
metadata.ITS<-read.csv("metadata_ITS_clean.csv", header = TRUE, row.names = 1)


#Make phyloseq with count data
phyloseq16s<-phyloseq(otu_table(asv_16s, taxa_are_rows = FALSE), tax_table(taxon_16s), sample_data(metadata.16s))
phyloseqITS<-phyloseq(otu_table(asv_ITS, taxa_are_rows = FALSE), tax_table(taxon_ITS), sample_data(metadata.ITS))

#Calculate relative abundances
#Replace zero values before clr transformation: Uses Count Zero Multiplicative method to replace zeros and outputs pseudo-counts
asv.n0_16s<-t(cmultRepl(asv_16s, label=0, method="CZM", output="p-counts")) #1740394 corrected values
asv.n0_ITS<-t(cmultRepl(asv_ITS, label=0, method="CZM", output="p-counts")) #989240  corrected values

#Note: Check the output to make sure there are no negative numbers. If samples or asv are sparce, the CZM method will 
#add a negative number that interferes with the log normalization. If necessary use function elow to convert negative values
#into positives
asv_n0_16s<-ifelse(asv.n0_16s < 0, asv.n0_16s*(-1), asv.n0_16s)
asv_n0_ITS<-ifelse(asv.n0_ITS < 0, asv.n0_ITS*(-1), asv.n0_ITS)


#Transform sample counts into relative abundances
asv.n0.acomp_16s<-as.data.frame(acomp(t(asv_n0_16s)), total=1)
asv.n0.acomp_ITS<-as.data.frame(acomp(t(asv_n0_ITS)), total=1)

#Make Phyloseq with relative abundances
phyloseq16s_rel = phyloseq(otu_table(asv.n0.acomp_16s, taxa_are_rows = FALSE), tax_table(taxon_16s), sample_data(metadata.16s))
phyloseqITS_rel = phyloseq(otu_table(asv.n0.acomp_ITS, taxa_are_rows = FALSE), tax_table(taxon_ITS), sample_data(metadata.ITS))

#prune taxa to remove ASVS that have less than 5 reads in the ASV table
phyloseq16s_prune<-prune_taxa(taxa_sums(phyloseq16s)>29, phyloseq16s_rel)
phyloseqITS_prune<-prune_taxa(taxa_sums(phyloseqITS)>29, phyloseqITS_rel)


# Identify bacterial families associated with the presence of L. monocytogenes samples.
model_lm_16s <- as.data.frame(cbind(otu_table(phyloseq16s_prune), sample_data(phyloseq16s_prune)[,5])) #Number in brackets is the column for metadata info with Lm presence
model_lm_16s$L.monocytogenes <- as.factor(model_lm_16s$L.monocytogenes)

model_lm_ITS <- as.data.frame(cbind(otu_table(phyloseqITS_prune), sample_data(phyloseqITS_prune)[,5])) #Number in brackets is the column for metadata info with Lm presence
model_lm_ITS$L.monocytogenes <- as.factor(model_lm_ITS$L.monocytogenes)

# Load Dataset
#16s
x_16s <- model_lm_16s[,1:ncol(model_lm_16s)-1] #ASV table
y_16s <- model_lm_16s[,ncol(model_lm_16s)] #Lm presence column is used as classifier

dataset_rf_16s <- as.data.frame(cbind(x_16s,y_16s))

#ITS
x_ITS <- model_lm_ITS[,1:ncol(model_lm_ITS)-1] #ASV table
y_ITS <- model_lm_ITS[,ncol(model_lm_ITS)] #Here, you can identify which column will be used as classifier. Lm presence column

dataset_rf_ITS <- as.data.frame(cbind(x_ITS,y_ITS))


#Randomforest parameter tuning script
#ntree is the number of trees that RF is going to make
#mtry is the number of trees that are going to be used in each split. If mtry is too low the tree might not be reliable, if too high you are overfitting

customRF <- list(type = "Classification", library = "randomForest", loop = NULL)
customRF$parameters <- data.frame(parameter = c("mtry", "ntree"), class = rep("numeric", 2), label = c("mtry", "ntree"))
customRF$grid <- function(x, y, len = NULL, search = "grid") {}
customRF$fit <- function(x, y, wts, param, lev, last, weights, classProbs, ...) {
  randomForest(x, y, mtry = param$mtry, ntree=param$ntree, ...) 
  #You can change other randomforest (i.e., conditional forest) here
}
customRF$predict <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
  predict(modelFit, newdata)
customRF$prob <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
  predict(modelFit, newdata, type = "prob")
customRF$sort <- function(x) x[order(x[,1]),]
customRF$levels <- function(x) x$classes


# train model
control <- trainControl(method="repeatedcv", number=10, repeats=3) #Change parameters (i.e., number, repeats)
#Repeating the cross validation to verify the accuracy of the model. number= will take that percentage of samples and reclassify to see if they fit in the same place
#For actual model, number=10, repeats=10
#number=10, repeats=3 is the standard


#16s
tunegrid_16s <- expand.grid(.mtry=c(105:155), .ntree=c(1000, 1500, 2000)) #mtry, ntree tuning parameters
#mtry sqroot of #of ASV is mean, add 20 to each side to establish parameter
custom_16s <- train(y_16s~., data=dataset_rf_16s, method=customRF, tuneGrid=tunegrid_16s, trControl=control)
summary(custom_16s)
png(file="RandomForest_train16s2.png")
plot(custom_16s)
dev.off()

#ITS
tunegrid_ITS <- expand.grid(.mtry=c(21:61), .ntree=c(1000, 1500, 2000))
custom_ITS <- train(y_ITS~., data=dataset_rf_ITS, method=customRF, tuneGrid=tunegrid_ITS, trControl=control)
summary(custom_ITS)
png(file="RandomForest_trainITS2.png")
plot(custom_ITS)
dev.off()


#See the best tuned paramter mtry and ntree
custom_16s$bestTune 
custom_ITS$bestTune 
#Accuracy of the model - over 80% is a stable/ reliable model

#Get accuracy of the RF models
max(custom_16s$results$Accuracy)
max(custom_ITS$results$Accuracy)

#16s
rf_classifier_16s <- randomForest(y_16s ~ ., data= dataset_rf_16s, ntree=1500, mtry=153, importance=TRUE) #mtry and mtree based on previous line
varimp_16s <- varImpPlot(rf_classifier_16s) #shows variable importance - which ASV has the highest classification power
varimp_16s <- as.data.frame(varimp_16s) 

varimp_16s = varimp_16s[order(varimp_16s$MeanDecreaseGini, na.last=NA, decreasing = FALSE),] #purity of samples - the higher the better

row_varimp_16s<-rownames(varimp_16s)

#ITS
rf_classifier_ITS <- randomForest(y_ITS ~ ., data= dataset_rf_ITS, ntree=1000, mtry=26, importance=TRUE) #mtry and mtree based on previous line
varimp_ITS <- varImpPlot(rf_classifier_ITS) #shows variable importance - which ASV has the highest classification power
varimp_ITS <- as.data.frame(varimp_ITS) 

f = varimp_ITS[order(varimp_ITS$MeanDecreaseGini, na.last=NA, decreasing = FALSE),] #purity of samples - the higher the better

row_varimp_ITS<-rownames(varimp_ITS)

#Replot accuracy and gini (aka purity) results
#16s
#Subset 30 ASVs with highest mean decrease accuracy
library(dplyr)
accuracy_16s<-varimp_16s[,1, drop=FALSE]

accuracy_16s_top30<-accuracy_16s %>%
  top_n(30)

accuracy_16s_top30$OTU<-rownames(accuracy_16s_top30)  
accuracy_16s_top30<-arrange(accuracy_16s_top30, OTU)
taxon_16s<-arrange(taxon_16s, rownames(taxon_16s))
accuracy_16s_top30$Genus<-subset(taxon_16s$Genus, rownames(taxon_16s) %in% rownames(accuracy_16s_top30))
accuracy_16s_top30$Seq<-rep("Bacteria",30)


library(ggplot2)
library(svglite)

#Fig 5D
accuracy_16s_top30_plot<-ggplot(accuracy_16s_top30, aes(x=MeanDecreaseAccuracy, y=reorder(OTU,MeanDecreaseAccuracy), color=Genus))+
  geom_point(stat='identity', size=5)+geom_text(aes(label=Genus, size=5), nudge_x = 0.7)+
  theme(axis.title.y = element_blank(), axis.text=element_text(color='black', size=13), axis.ticks=element_line(color='black')) +
  xlab("Mean Deacrease Accuracy")+
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  theme(panel.background = element_rect(fill="grey99", color =NA),
        plot.background = element_rect(fill="white", color =NA)) +
  theme(strip.background= element_rect(fill="white", color = "black"), strip.text.y = element_text(size=6, angle=0), 
        panel.border = element_rect(color="black", fill=NA))+
  scale_color_viridis_d(begin=0.1, end=0.9, option='magma')+
  theme(legend.position = 'none')+scale_x_continuous(limits = c(2,8), breaks = c(2,4,6,8))+
  ggtitle("Bacteria")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=15))
ggsave("RandomForest_16s_top30RA2.png", plot =accuracy_16s_top30_plot, device="png", width=10, height=8, units="in",dpi=600)
ggsave("RandomForest_16s_top30RA2.svg", plot =accuracy_16s_top30_plot, device="svg", width=10, height=8, units="in",dpi=600)



#ITS
#Subset 30 ASVs with highest mean decrease accuracy
accuracy_ITS<-varimp_ITS[,1, drop=FALSE]

accuracy_ITS_top30<-accuracy_ITS %>%
  top_n(30)

accuracy_ITS_top30$OTU<-rownames(accuracy_ITS_top30)  
accuracy_ITS_top30<-arrange(accuracy_ITS_top30, OTU)
taxon_ITS<-arrange(taxon_ITS, rownames(taxon_ITS))
accuracy_ITS_top30$Genus<-subset(taxon_ITS$Genus, rownames(taxon_ITS) %in% rownames(accuracy_ITS_top30))
accuracy_ITS_top30$Seq<-rep("Fungi",30)

#Fig 6D
accuracy_ITS_top30_plot<-ggplot(accuracy_ITS_top30, aes(x=MeanDecreaseAccuracy, y=reorder(OTU,MeanDecreaseAccuracy), color=Genus))+
  geom_point(stat='identity', size=5)+geom_text(aes(label=Genus, size=5), nudge_x = 0.5)+
  theme(axis.title.y = element_blank(), axis.text=element_text(color='black', size=13), axis.ticks=element_line(color='black')) +
  xlab("Mean Deacrease Accuracy")+
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  theme(panel.background = element_rect(fill="grey99", color =NA),
        plot.background = element_rect(fill="white", color =NA)) +
  theme(strip.background= element_rect(fill="white", color = "black"), strip.text.y = element_text(size=6, angle=0), 
        panel.border = element_rect(color="black", fill=NA))+
  scale_color_viridis_d(begin=0.1, end=0.9, option='viridis')+
  theme(legend.position = 'none')+scale_x_continuous(limits = c(2,10), breaks = c(2,4,6,8,10))+
  ggtitle("Fungi")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=15))
ggsave("RandomForest_ITS_top30RA2.png", plot =accuracy_ITS_top30_plot, device="png", width=10, height=8, units="in",dpi=600)
ggsave("RandomForest_ITS_top30RA2.svg", plot =accuracy_ITS_top30_plot, device="svg", width=10, height=8, units="in",dpi=600)



## ROC curve and AUC analysis
## Out of Bag (OOB) sample tree classification votes for each data point. These votes roughly represent a probability, and therefore can be used to create a ROC and AUC measure
library(ROCR)
# 16S rRNA
predictions_16s=as.vector(rf_classifier_16s$votes[,2])
pred_16s=prediction(predictions_16s,dataset_rf_16s$y_16s)

perf_AUC_16s=performance(pred_16s,"auc") #Calculate the AUC value
AUC_16s=perf_AUC_16s@y.values[[1]]
kappa_16s = caret::confusionMatrix(rf_classifier_16s$predicted, rf_classifier_16s$y)$overall[["Kappa"]]
perf_ROC_16s=performance(pred_16s,"tpr","fpr") #plot the actual ROC curve
performanceAUC_16s<-data.frame("FalsePositive"=unlist(perf_ROC_16s@x.values),"TruePositive"=unlist(perf_ROC_16s@y.values))

#Fig 5D insert
AUCplot_16s<-ggplot(performanceAUC_16s, aes(x=FalsePositive, y=TruePositive))+geom_line()+
  theme(axis.text=element_text(color='black', size=13), axis.ticks=element_line(color='black')) +
  xlab("False Positive Rate")+ylab("True Positive Rate")+
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  theme(panel.background = element_rect(fill="grey99", color ="black"),
        plot.background = element_rect(fill="white", color =NA),panel.border = element_rect(color="black", fill=NA))+
  geom_abline(intercept=0, linetype="dashed", color="red")+
  annotate("text", x = 0.8, y=0.2, label = paste("AUC = ",format(AUC_16s, digits=4, scientific=FALSE)))+
  annotate("text", x = 0.8, y=0.1, label = paste("Kappa = ",format(kappa_16s, digits=4, scientific=FALSE)))+
  ggtitle("Bacteria")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))


#ITS
predictions_ITS = as.vector(rf_classifier_ITS$votes[,2])
pred_ITS=prediction(predictions_ITS,dataset_rf_ITS$y_ITS)

perf_AUC_ITS=performance(pred_ITS,"auc") #Calculate the AUC value
AUC_ITS=perf_AUC_ITS@y.values[[1]]
kappa_ITS= caret::confusionMatrix(rf_classifier_ITS$predicted, rf_classifier_ITS$y)$overall[["Kappa"]]
perf_ROC_ITS=performance(pred_ITS,"tpr","fpr") #plot the actual ROC curve
performanceAUC_ITS<-data.frame("FalsePositive"=unlist(perf_ROC_ITS@x.values),"TruePositive"=unlist(perf_ROC_ITS@y.values))

#Fig 6D insert
AUCplot_ITS<-ggplot(performanceAUC_ITS, aes(x=FalsePositive, y=TruePositive))+geom_line()+
  theme(axis.text=element_text(color='black', size=13), axis.ticks=element_line(color='black')) +
  xlab("False Positive Rate")+ylab("True Positive Rate")+
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  theme(panel.background = element_rect(fill="grey99", color ="black"),
        plot.background = element_rect(fill="white", color =NA),panel.border = element_rect(color="black", fill=NA))+
  geom_abline(intercept=0, linetype="dashed", color="red")+
  annotate("text", x = 0.8, y=0.2, label = paste("AUC = ",format(AUC_ITS, digits=4, scientific=FALSE)))+
  annotate("text", x = 0.8, y=0.1, label = paste("Kappa = ",format(kappa_ITS, digits=4, scientific=FALSE)))+
  ggtitle("Fungi")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))

library(cowplot)
RandomForest_AUC = plot_grid(AUCplot_16s, AUCplot_ITS, 
                             ncol=2, nrow=1)
RandomForest_AUC
ggsave("RandomForest_AUC2.png", plot =RandomForest_AUC, device="png", width=16, height=8, units="in",dpi=600)
ggsave("RandomForest_AUC2.svg", plot =RandomForest_AUC, device="svg", width=16, height=8, units="in",dpi=600)


#Check RA of most important predictors
#Import asv table with RA
asv_16s_long<-read.csv("asv_16s.csv", header = TRUE, row.names = 1)
asv_ITS_long<-read.csv("asv_ITS.csv", header = TRUE, row.names = 1)

#Make string with 30 top predictor ASVs
asv_16s_pred<-accuracy_16s_top30$OTU
asv_ITS_pred<-accuracy_ITS_top30$OTU

#Filter ASV table to keep only ASV predictors
asv_16s_pred_RA<-filter(asv_16s_long, OTU %in% asv_16s_pred)
asv_ITS_pred_RA<-filter(asv_ITS_long, OTU %in% asv_ITS_pred)

#Calculate mean RA by L. monocytogenes presence/absence
asv_16s_pred_mean<-asv_16s_pred_RA%>%
  group_by(OTU,Genus,L.monocytogenes)%>%
  summarize(Mean=mean(Abundance))%>%
  arrange(desc(Mean))

asv_ITS_pred_mean<-asv_ITS_pred_RA%>%
  group_by(OTU,Genus,L.monocytogenes)%>%
  summarize(Mean=mean(Abundance))%>%
  arrange(desc(Mean))

#Make plot by mean RA
asv_16s_Lm_mean<-ggplot(asv_16s_pred_mean, aes(x=L.monocytogenes, y=Mean, fill=L.monocytogenes))+
  geom_bar(stat='identity', color='black')+facet_wrap(~OTU, scales = "free")
ggsave(plot=asv_16s_Lm_mean, "RFpredictors_16s_Lm_mean.png", device="png", width=16, height=8, units="in",dpi=600)

asv_ITS_Lm_mean<-ggplot(asv_ITS_pred_mean, aes(x=L.monocytogenes, y=Mean, fill=L.monocytogenes))+
  geom_bar(stat='identity', color='black')+facet_wrap(~OTU, scales = "free")
ggsave(plot=asv_ITS_Lm_mean, "RFpredictors_ITS_Lm_mean.png", device="png", width=16, height=8, units="in",dpi=600)

#Make boxplot including all data
asv_16s_Lm<-ggplot(asv_16s_pred_RA, aes(x=L.monocytogenes, y=Abundance, fill=L.monocytogenes))+
  geom_boxplot(color='black')+facet_wrap(~OTU, scales = "free")
ggsave(plot=asv_16s_Lm, "RFpredictors_16s_Lm.png", device="png", width=16, height=8, units="in",dpi=600)

asv_ITS_Lm<-ggplot(asv_ITS_pred_RA, aes(x=L.monocytogenes, y=Abundance, fill=L.monocytogenes))+
  geom_boxplot(color='black')+facet_wrap(~OTU, scales = "free")
ggsave(plot=asv_ITS_Lm, "RFpredictors_ITS_Lm.png", device="png", width=16, height=8, units="in",dpi=600)


#### Random Forest model without F2 ####
#Prune samples to remove F2
phyloseq16s_F1F3<-subset_samples(phyloseq16s_prune, Facility != "F2")
phyloseqITS_F1F3<-subset_samples(phyloseqITS_prune, Facility != "F2")

metadata_16s_F1F3<-subset(metadata.16s, Facility!="F2")
metadata_ITS_F1F3<-subset(metadata.ITS, Facility!="F2")


# Identify bacterial families associated with the presence of Salmonella in water samples.
model_lm_16s_F1F3 <- as.data.frame(cbind(otu_table(phyloseq16s_F1F3), sample_data(phyloseq16s_F1F3)[,5])) #Number in brackets is the column for metadata info with Lm presence
model_lm_16s_F1F3$L.monocytogenes <- as.factor(model_lm_16s_F1F3$L.monocytogenes)

model_lm_ITS_F1F3 <- as.data.frame(cbind(otu_table(phyloseqITS_F1F3), sample_data(phyloseqITS_F1F3)[,5])) #Number in brackets is the column for metadata info with Lm presence
model_lm_ITS_F1F3$L.monocytogenes <- as.factor(model_lm_ITS_F1F3$L.monocytogenes)

# Load Dataset
#16s
x_16s_F1F3 <- model_lm_16s_F1F3[,1:ncol(model_lm_16s_F1F3)-1] #ASV table
y_16s_F1F3 <- model_lm_16s_F1F3[,ncol(model_lm_16s_F1F3)] #Here, you can identify which column will be used as classifier. Lm presence column

dataset_rf_16s_F1F3 <- as.data.frame(cbind(x_16s_F1F3,y_16s_F1F3))

#ITS
x_ITS_F1F3 <- model_lm_ITS_F1F3[,1:ncol(model_lm_ITS_F1F3)-1] #ASV table
y_ITS_F1F3 <- model_lm_ITS_F1F3[,ncol(model_lm_ITS_F1F3)] #Here, you can identify which column will be used as classifier. Lm presence column

dataset_rf_ITS_F1F3 <- as.data.frame(cbind(x_ITS_F1F3,y_ITS_F1F3))


# train model
#16s
tunegrid_16s_F1F3 <- expand.grid(.mtry=c(105:155), .ntree=c(1000, 1500, 2000)) #mtry, ntree tuning parameters
#mtry sqroot of #of ASV is mean, add 20 to each side to establish parameter
set.seed(1001)
custom_16s_F1F3 <- train(y_16s_F1F3~., data=dataset_rf_16s_F1F3, method=customRF, tuneGrid=tunegrid_16s_F1F3, trControl=control)
summary(custom_16s_F1F3)

#ITS
tunegrid_ITS_F1F3 <- expand.grid(.mtry=c(21:61), .ntree=c(1000, 1500, 2000))
custom_ITS_F1F3 <- train(y_ITS_F1F3~., data=dataset_rf_ITS_F1F3, method=customRF, tuneGrid=tunegrid_ITS_F1F3, trControl=control)
summary(custom_ITS_F1F3)


#See the best tuned paramter mtry and ntree
custom_16s_F1F3$bestTune 
custom_ITS_F1F3$bestTune 
#Accuracy of the model - over 80% is a stable/ reliable model

#Get accuracy of the RF models
max(custom_16s_F1F3$results$Accuracy)
max(custom_ITS_F1F3$results$Accuracy)


#16s
rf_classifier_16s_F1F3 <- randomForest(y_16s_F1F3 ~ ., data= dataset_rf_16s_F1F3, ntree=1000, mtry=114, importance=TRUE) #mtry and mtree based on previous line
varimp_16s_F1F3 <- varImpPlot(rf_classifier_16s_F1F3) #shows variable importance - which ASV has the highest classification power
varimp_16s_F1F3 <- as.data.frame(varimp_16s_F1F3) 

f = varimp_16s_F1F3[order(varimp_16s_F1F3$MeanDecreaseGini, na.last=NA),] #purity of samples - the higher the better

row_varimp_16s_F1F3<-rownames(varimp_16s_F1F3)

#ITS
rf_classifier_ITS_F1F3 <- randomForest(y_ITS_F1F3 ~ ., data= dataset_rf_ITS_F1F3, ntree=1500, mtry=29, importance=TRUE) #mtry and mtree based on previous line
varimp_ITS_F1F3 <- varImpPlot(rf_classifier_ITS_F1F3) #shows variable importance - which ASV has the highest classification power
varimp_ITS_F1F3 <- as.data.frame(varimp_ITS_F1F3) 

f = varimp_ITS_F1F3[order(varimp_ITS_F1F3$MeanDecreaseGini, na.last=NA),] #purity of samples - the higher the better

row_varimp_ITS_F1F3<-rownames(varimp_ITS_F1F3)

#Replot accuracy results
#16s
#Subset 10 ASVs with highest mean decrease accuracy
library(dplyr)
accuracy_16s_F1F3<-varimp_16s_F1F3[,1, drop=FALSE]

accuracy_16s_F1F3_top30<-accuracy_16s_F1F3 %>%
  top_n(30)

accuracy_16s_F1F3_top30$OTU<-rownames(accuracy_16s_F1F3_top30)  
accuracy_16s_F1F3_top30<-arrange(accuracy_16s_F1F3_top30, OTU)
accuracy_16s_F1F3_top30$Genus<-subset(taxon_16s$Genus, rownames(taxon_16s) %in% rownames(accuracy_16s_F1F3_top30))
accuracy_16s_F1F3_top30$Seq<-rep("Bacteria",30)


library(ggplot2)
accuracy_16s_F1F3_top30_plot<-ggplot(accuracy_16s_F1F3_top30, aes(x=MeanDecreaseAccuracy, y=reorder(OTU,MeanDecreaseAccuracy), color=Genus))+
  geom_point(stat='identity', size=5)+geom_text(aes(label=Genus, size=5), nudge_x = 0.7)+
  theme(axis.title.y = element_blank(), axis.text=element_text(color='black', size=13), axis.ticks=element_line(color='black')) +
  xlab("Mean Deacrease Accuracy")+
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  theme(panel.background = element_rect(fill="grey99", color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_rect(fill="white", color = "black"), strip.text.y = element_text(size=6, angle=0), 
        panel.border = element_rect(color="black", fill=NA))+
  scale_color_viridis_d(begin=0.1, end=0.9, option='magma')+
  theme(legend.position = 'none')+scale_x_continuous(limits = c(0,10), breaks = c(0,2,4,6,8,10))+
  ggtitle("Bacteria")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=15))


#Subset 30 samples with highest mean decrease gini
gini_16s_F1F3<-varimp_16s_F1F3[,2, drop=FALSE]

gini_16s_F1F3_top30<-gini_16s_F1F3 %>%
  top_n(30)

gini_16s_F1F3_top30$OTU<-rownames(gini_16s_F1F3_top30)  
gini_16s_F1F3_top30<-arrange(gini_16s_F1F3_top30, OTU)
gini_16s_F1F3_top30$Genus<-subset(taxon_16s$Genus, rownames(taxon_16s) %in% rownames(gini_16s_F1F3_top30))


gini_16s_F1F3_top30_plot<-ggplot(gini_16s_F1F3_top30, aes(x=MeanDecreaseGini, y=reorder(OTU,MeanDecreaseGini), color=Genus))+
  geom_point(stat='identity')+
  geom_point(stat='identity', size=5)+geom_text(aes(label=Genus, size=5), nudge_x = 0.7)+
  theme(axis.title.y = element_blank(), axis.text=element_text(color='black', size=13), axis.ticks=element_line(color='black')) +
  xlab("Mean Deacrease Gini")+
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  theme(panel.background = element_rect(fill="grey99", color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_rect(fill="white", color = "black"), strip.text.y = element_text(size=6, angle=0), 
        panel.border = element_rect(color="black", fill=NA))+
  scale_color_viridis_d(begin=0.1, end=0.9, option='magma')+
  theme(legend.position = 'none')+scale_x_continuous(limits = c(0,10), breaks = c(0,2,4,6,8,10))+
  ggtitle("Bacteria")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=15))

#Combine plots
library(cowplot)
RandomForest_16s_F1F3 = plot_grid(accuracy_16s_F1F3_top30_plot, gini_16s_F1F3_top30_plot, 
                             ncol=2, nrow=1)
RandomForest_16s_F1F3
ggsave("RandomForest_16s_F1F3_top30.png", plot =RandomForest_16s_F1F3, device="png", width=16, height=8, units="in",dpi=600)

#ITS
#Subset 30 ASVs with highest mean decrease accuracy
accuracy_ITS_F1F3<-varimp_ITS_F1F3[,1, drop=FALSE]

accuracy_ITS_F1F3_top30<-accuracy_ITS_F1F3 %>%
  top_n(30)

accuracy_ITS_F1F3_top30$OTU<-rownames(accuracy_ITS_F1F3_top30)  
accuracy_ITS_F1F3_top30<-arrange(accuracy_ITS_F1F3_top30, OTU)
accuracy_ITS_F1F3_top30$Genus<-subset(taxon_ITS$Genus, rownames(taxon_ITS) %in% rownames(accuracy_ITS_F1F3_top30))
accuracy_ITS_F1F3_top30$Seq<-rep("Fungi",30)

accuracy_ITS_F1F3_top30_plot<-ggplot(accuracy_ITS_F1F3_top30, aes(x=MeanDecreaseAccuracy, y=reorder(OTU,MeanDecreaseAccuracy), color=Genus))+
  geom_point(stat='identity', size=5)+geom_text(aes(label=Genus, size=5), nudge_x = 0.5)+
  theme(axis.title.y = element_blank(), axis.text=element_text(color='black', size=13), axis.ticks=element_line(color='black')) +
  xlab("Mean Deacrease Accuracy")+
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  theme(panel.background = element_rect(fill="grey99", color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_rect(fill="white", color = "black"), strip.text.y = element_text(size=6, angle=0), 
        panel.border = element_rect(color="black", fill=NA))+
  scale_color_viridis_d(begin=0.1, end=0.9, option='viridis')+
  theme(legend.position = 'none')+scale_x_continuous(limits = c(0,10), breaks = c(0,2,4,6,8,10))+
  ggtitle("Fungi")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=15))


#Subset 30 samples with highest mean decrease gini
gini_ITS_F1F3<-varimp_ITS_F1F3[,2, drop=FALSE]

gini_ITS_F1F3_top30<-gini_ITS_F1F3 %>%
  top_n(30)

gini_ITS_F1F3_top30$OTU<-rownames(gini_ITS_F1F3_top30)  
gini_ITS_F1F3_top30<-arrange(gini_ITS_F1F3_top30, OTU)
gini_ITS_F1F3_top30$Genus<-subset(taxon_ITS$Genus, rownames(taxon_ITS) %in% rownames(gini_ITS_F1F3_top30))
gini_ITS_F1F3_top30$Seq<-rep("Fungi",30)


gini_ITS_F1F3_top30_plot<-ggplot(gini_ITS_F1F3_top30, aes(x=MeanDecreaseGini, y=reorder(OTU,MeanDecreaseGini), color=Genus))+
  geom_point(stat='identity')+
  geom_point(stat='identity', size=5)+geom_text(aes(label=Genus, size=5), nudge_x = 0.7)+
  theme(axis.title.y = element_blank(), axis.text=element_text(color='black', size=13), axis.ticks=element_line(color='black')) +
  xlab("Mean Deacrease Gini")+
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  theme(panel.background = element_rect(fill="grey99", color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_rect(fill="white", color = "black"), strip.text.y = element_text(size=6, angle=0), 
        panel.border = element_rect(color="black", fill=NA))+
  scale_color_viridis_d(begin=0.1, end=0.9, option='viridis')+
  theme(legend.position = 'none')+scale_x_continuous(limits = c(0,10), breaks = c(0,2,4,6,8,10))+
  ggtitle("Fungi")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=15))

#Combine plots
RandomForest_ITS_F1F3 = plot_grid(accuracy_ITS_F1F3_top30_plot, gini_ITS_F1F3_top30_plot, 
                             ncol=2, nrow=1)
RandomForest_ITS_F1F3
ggsave("RandomForest_ITS_top30.png", plot =RandomForest_ITS_F1F3, device="png", width=16, height=8, units="in",dpi=600)

#Combine accuracy for 16s and ITS


RandomForest_accuracy_F1F3 = plot_grid(accuracy_16s_F1F3_top30_plot, accuracy_ITS_F1F3_top30_plot, 
                                   ncol=2, nrow=1)
RandomForest_accuracy_F1F3
ggsave("RandomForest_top30_F1F3_FINAL.png", plot =RandomForest_accuracy_F1F3, device="png", width=16, height=8, units="in",dpi=600)
ggsave("RandomForest_top30_F1F3_FINAL.svg", plot =RandomForest_accuracy_F1F3, device="svg", width=16, height=8, units="in",dpi=600)



## ROC curve and AUC analysis
## Out of Bag (OOB) sample tree classification votes for each data point. These votes roughly represent a probability, and therefore can be used to create a ROC and AUC measure
library(ROCR)
# 16S rRNA
predictions_16s_F1F3=as.vector(rf_classifier_16s_F1F3$votes[,2])
pred_16s_F1F3=prediction(predictions_16s_F1F3,dataset_rf_16s_F1F3$y_16s_F1F3)

perf_AUC_16s_F1F3=performance(pred_16s_F1F3,"auc") #Calculate the AUC value
AUC_16s_F1F3=perf_AUC_16s_F1F3@y.values[[1]]
kappa_16s_F1F3 = caret::confusionMatrix(rf_classifier_16s_F1F3$predicted, rf_classifier_16s_F1F3$y)$overall[["Kappa"]]
perf_ROC_16s_F1F3=performance(pred_16s_F1F3,"tpr","fpr") #plot the actual ROC curve
performanceAUC_16s_F1F3<-data.frame("FalsePositive"=unlist(perf_ROC_16s_F1F3@x.values),"TruePositive"=unlist(perf_ROC_16s_F1F3@y.values))


AUCplot_16s_F1F3<-ggplot(performanceAUC_16s_F1F3, aes(x=FalsePositive, y=TruePositive))+geom_line()+
  theme(axis.text=element_text(color='black', size=13), axis.ticks=element_line(color='black')) +
  xlab("False Positive Rate")+ylab("True Positive Rate")+
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  theme(panel.background = element_rect(fill="grey99", color ="black"),
        plot.background = element_rect(fill="transparent", color =NA),panel.border = element_rect(color="black", fill=NA))+
  geom_abline(intercept=0, linetype="dashed", color="red")+
  annotate("text", x = 0.8, y=0.2, label = paste("AUC = ",format(AUC_16s_F1F3, digits=4, scientific=FALSE)))+
  annotate("text", x = 0.8, y=0.1, label = paste("Kappa = ",format(kappa_16s_F1F3, digits=4, scientific=FALSE)))+
  ggtitle("Bacteria")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))

#ITS
predictions_ITS_F1F3 = as.vector(rf_classifier_ITS_F1F3$votes[,2])
pred_ITS_F1F3=prediction(predictions_ITS_F1F3,dataset_rf_ITS_F1F3$y_ITS)

perf_AUC_ITS_F1F3=performance(pred_ITS_F1F3,"auc") #Calculate the AUC value
AUC_ITS_F1F3=perf_AUC_ITS_F1F3@y.values[[1]]
kappa_ITS_F1F3= caret::confusionMatrix(rf_classifier_ITS_F1F3$predicted, rf_classifier_ITS_F1F3$y)$overall[["Kappa"]]
perf_ROC_ITS_F1F3=performance(pred_ITS_F1F3,"tpr","fpr") #plot the actual ROC curve
performanceAUC_ITS_F1F3<-data.frame("FalsePositive"=unlist(perf_ROC_ITS_F1F3@x.values),"TruePositive"=unlist(perf_ROC_ITS_F1F3@y.values))


AUCplot_ITS_F1F3<-ggplot(performanceAUC_ITS_F1F3, aes(x=FalsePositive, y=TruePositive))+geom_line()+
  theme(axis.text=element_text(color='black', size=13), axis.ticks=element_line(color='black')) +
  xlab("False Positive Rate")+ylab("True Positive Rate")+
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  theme(panel.background = element_rect(fill="grey99", color ="black"),
        plot.background = element_rect(fill="transparent", color =NA),panel.border = element_rect(color="black", fill=NA))+
  geom_abline(intercept=0, linetype="dashed", color="red")+
  annotate("text", x = 0.8, y=0.2, label = paste("AUC = ",format(AUC_ITS_F1F3, digits=4, scientific=FALSE)))+
  annotate("text", x = 0.8, y=0.1, label = paste("Kappa = ",format(kappa_ITS_F1F3, digits=4, scientific=FALSE)))+
  ggtitle("Fungi")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))

RandomForest_AUC_F1F3 = plot_grid(AUCplot_16s_F1F3, AUCplot_ITS_F1F3, 
                             ncol=2, nrow=1)
RandomForest_AUC_F1F3
ggsave("RandomForest_AUC_F1F3.png", plot =RandomForest_AUC_F1F3, device="png", width=16, height=8, units="in",dpi=600)
ggsave("RandomForest_AUC_F1F3.svg", plot =RandomForest_AUC_F1F3, device="svg", width=16, height=8, units="in",dpi=600)



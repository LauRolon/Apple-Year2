#Random Forest Analysis of microbiome data

#Looking for OTUs that classify samples into Lm positive or negative
#Written by Taejung Chung
#Updated by Laura Rolon / Taejung Chung
#Last updated: 10/19/21 by TCC

#Load packages
library(ape) 
library(phyloseq)
library(randomForest) 
library(caret) 
library(svglite)
library(zCompositions)
library(compositions)

#Set working directory to where files are located
setwd()

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

#Make phyloseq with count data
phyloseq16s<-phyloseq(otu_table(otus_16s, taxa_are_rows = TRUE), TAX_16s, META_16s)
phyloseqITS<-phyloseq(otu_table(otus_ITS, taxa_are_rows = TRUE), TAX_ITS, META_ITS)

#Calculate relative abundances
#Replace zero values before clr transformation: Uses Count Zero Multiplicative method to replace zeros and outputs pseudo-counts
otu.n0_16s<-t(cmultRepl(t(otus_16s), label=0, method="CZM", output="p-counts")) #1740394 corrected values
otu.n0_ITS<-t(cmultRepl(t(otus_ITS), label=0, method="CZM", output="p-counts")) #989240  corrected values

#Transform sample counts into relative abundances
otu.n0.acomp_16s<-as.data.frame(acomp(t(otu.n0_16s)), total=1)
otu.n0.acomp_ITS<-as.data.frame(acomp(t(otu.n0_ITS)), total=1)

#Make Phyloseq with relative abundances
OTU_16s <- otu_table(otu.n0.acomp_16s, taxa_are_rows = FALSE)
phyloseq16s_rel = phyloseq(OTU_16s, TAX_16s, META_16s)

OTU_ITS <- otu_table(otu.n0.acomp_ITS, taxa_are_rows = FALSE)
phyloseqITS_rel = phyloseq(OTU_ITS, TAX_ITS, META_ITS)


#prune taxa to remove OTUS that have less than 5 reads in the OTU table
phyloseq16s_prune<-prune_taxa(taxa_sums(phyloseq16s)>4, phyloseq16s_rel)
phyloseqITS_prune<-prune_taxa(taxa_sums(phyloseqITS)>4, phyloseqITS_rel)


# Identify bacterial families associated with the presence of /l. monocytogenes samples.
model_lm_16s <- as.data.frame(cbind(otu_table(phyloseq16s_prune), sample_data(phyloseq16s_prune)[,5])) #Number in brackets is the column for metadata info with Lm presence
model_lm_16s$L..monocytogenes <- as.factor(model_lm_16s$L..monocytogenes)

model_lm_ITS <- as.data.frame(cbind(otu_table(phyloseqITS_prune), sample_data(phyloseqITS_prune)[,5])) #Number in brackets is the column for metadata info with Lm presence
model_lm_ITS$L..monocytogenes <- as.factor(model_lm_ITS$L..monocytogenes)

# Load Dataset
#16s
x_16s <- model_lm_16s[,1:ncol(model_lm_16s)-1] #OTU table
y_16s <- model_lm_16s[,ncol(model_lm_16s)] #Lm presence column is used as classifier

dataset_rf_16s <- as.data.frame(cbind(x_16s,y_16s))

#ITS
x_ITS <- model_lm_ITS[,1:ncol(model_lm_ITS)-1] #OTU table
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
tunegrid_16s <- expand.grid(.mtry=c(60:110), .ntree=c(1000, 1500, 2000)) #mtry, ntree tuning parameters
#mtry sqroot of #of otu is mean, add 20 to each side to establish parameter
set.seed(1001)
custom_16s <- train(y_16s~., data=dataset_rf_16s, method=customRF, tuneGrid=tunegrid_16s, trControl=control)
summary(custom_16s)
plot(custom_16s)

#ITS
tunegrid_ITS <- expand.grid(.mtry=c(35:75), .ntree=c(1000, 1500, 2000))
custom_ITS <- train(y_ITS~., data=dataset_rf_ITS, method=customRF, tuneGrid=tunegrid_ITS, trControl=control)
summary(custom_ITS)
plot(custom_ITS)



#See the best tuned paramter mtry and ntree
custom_16s$bestTune 
custom_ITS$bestTune 
#Accuracy of the model - over 80% is a stable/ reliable model

#Get accuracy of the RF models
max(custom_16s$results$Accuracy)
max(custom_ITS$results$Accuracy)

#16s
rf_classifier_16s <- randomForest(y_16s ~ ., data= dataset_rf_16s, ntree=1000, mtry=88, importance=TRUE) #mtry and mtree based on previous line
varimp_16s <- varImpPlot(rf_classifier_16s) #shows variable importance - which OTU has the highest classification power
varimp_16s <- as.data.frame(varimp_16s) 

varimp_16s = varimp_16s[order(varimp_16s$MeanDecreaseGini, na.last=NA, decreasing = FALSE),] #purity of samples - the higher the better

row_varimp_16s<-rownames(varimp_16s)

#ITS
rf_classifier_ITS <- randomForest(y_ITS ~ ., data= dataset_rf_ITS, ntree=1500, mtry=69, importance=TRUE) #mtry and mtree based on previous line
varimp_ITS <- varImpPlot(rf_classifier_ITS) #shows variable importance - which OTU has the highest classification power
varimp_ITS <- as.data.frame(varimp_ITS) 

f = varimp_ITS[order(varimp_ITS$MeanDecreaseGini, na.last=NA, decreasing = FALSE),] #purity of samples - the higher the better

row_varimp_ITS<-rownames(varimp_ITS)

#Replot accuracy and gini (aka purity) results
#16s
#Subset 30 OTUs with highest mean decrease accuracy
library(dplyr)
accuracy_16s<-varimp_16s[,1, drop=FALSE]

accuracy_16s_top30<-accuracy_16s %>%
  top_n(30)

accuracy_16s_top30$OTU<-rownames(accuracy_16s_top30)  
accuracy_16s_top30<-arrange(accuracy_16s_top30, OTU)
accuracy_16s_top30$Family<-subset(taxon_16s$Family, rownames(taxon_16s) %in% rownames(accuracy_16s_top30))
accuracy_16s_top30$Seq<-rep("Bacteria",30)


library(ggplot2)
library(svglite)

#Fig 5D
accuracy_16s_top30_plot<-ggplot(accuracy_16s_top30, aes(x=MeanDecreaseAccuracy, y=reorder(OTU,MeanDecreaseAccuracy), color=Family))+
  geom_point(stat='identity', size=5)+geom_text(aes(label=Family, size=5), nudge_x = 0.7)+
  theme(axis.title.y = element_blank(), axis.text=element_text(color='black', size=13), axis.ticks=element_line(color='black')) +
  xlab("Mean Deacrease Accuracy")+
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  theme(panel.background = element_rect(fill="grey99", color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_rect(fill="white", color = "black"), strip.text.y = element_text(size=6, angle=0), 
        panel.border = element_rect(color="black", fill=NA))+
  scale_color_viridis_d(begin=0.1, end=0.9, option='magma')+
  theme(legend.position = 'none')+scale_x_continuous(limits = c(2,8), breaks = c(2,4,6,8))+
  ggtitle("Bacteria")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=15))
ggsave("RandomForest_16s_top30RA.png", plot =accuracy_16s_top30_plot, device="png", width=10, height=8, units="in",dpi=600)
ggsave("RandomForest_16s_top30RA.svg", plot =accuracy_16s_top30_plot, device="svg", width=10, height=8, units="in",dpi=600)



#ITS
#Subset 30 OTUs with highest mean decrease accuracy
accuracy_ITS<-varimp_ITS[,1, drop=FALSE]

accuracy_ITS_top30<-accuracy_ITS %>%
  top_n(30)

accuracy_ITS_top30$OTU<-rownames(accuracy_ITS_top30)  
accuracy_ITS_top30<-arrange(accuracy_ITS_top30, OTU)
accuracy_ITS_top30$Family<-subset(taxon_ITS$Family, rownames(taxon_ITS) %in% rownames(accuracy_ITS_top30))
accuracy_ITS_top30$Family<-gsub("f__","", accuracy_ITS_top30$Family)
accuracy_ITS_top30$Family<-gsub("o__","", accuracy_ITS_top30$Family)
accuracy_ITS_top30$Family<-gsub("p__","", accuracy_ITS_top30$Family)
accuracy_ITS_top30$Family<-gsub("c__","", accuracy_ITS_top30$Family)
accuracy_ITS_top30$Seq<-rep("Fungi",30)

#Fig 6D
accuracy_ITS_top30_plot<-ggplot(accuracy_ITS_top30, aes(x=MeanDecreaseAccuracy, y=reorder(OTU,MeanDecreaseAccuracy), color=Family))+
  geom_point(stat='identity', size=5)+geom_text(aes(label=Family, size=5), nudge_x = 0.5)+
  theme(axis.title.y = element_blank(), axis.text=element_text(color='black', size=13), axis.ticks=element_line(color='black')) +
  xlab("Mean Deacrease Accuracy")+
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  theme(panel.background = element_rect(fill="grey99", color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_rect(fill="white", color = "black"), strip.text.y = element_text(size=6, angle=0), 
        panel.border = element_rect(color="black", fill=NA))+
  scale_color_viridis_d(begin=0.1, end=0.9, option='viridis')+
  theme(legend.position = 'none')+scale_x_continuous(limits = c(2,8), breaks = c(2,4,6,8))+
  ggtitle("Fungi")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=15))
ggsave("RandomForest_ITS_top30RA.png", plot =accuracy_ITS_top30_plot, device="png", width=10, height=8, units="in",dpi=600)
ggsave("RandomForest_ITS_top30RA.svg", plot =accuracy_ITS_top30_plot, device="svg", width=10, height=8, units="in",dpi=600)



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
        plot.background = element_rect(fill="transparent", color =NA),panel.border = element_rect(color="black", fill=NA))+
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
        plot.background = element_rect(fill="transparent", color =NA),panel.border = element_rect(color="black", fill=NA))+
  geom_abline(intercept=0, linetype="dashed", color="red")+
  annotate("text", x = 0.8, y=0.2, label = paste("AUC = ",format(AUC_ITS, digits=4, scientific=FALSE)))+
  annotate("text", x = 0.8, y=0.1, label = paste("Kappa = ",format(kappa_ITS, digits=4, scientific=FALSE)))+
  ggtitle("Fungi")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))

library(cowplot)
RandomForest_AUC = plot_grid(AUCplot_16s, AUCplot_ITS, 
                             ncol=2, nrow=1)
RandomForest_AUC
ggsave("RandomForest_AUC.png", plot =RandomForest_AUC, device="png", width=16, height=8, units="in",dpi=600)
ggsave("RandomForest_AUC.svg", plot =RandomForest_AUC, device="svg", width=16, height=8, units="in",dpi=600)



#### Random Forest model without F2 ####
#Prune samples to remove F2
phyloseq16s_F1F3<-subset_samples(phyloseq16s_rel, Facility != "F2")
phyloseqITS_F1F3<-subset_samples(phyloseqITS_rel, Facility != "F2")

metadata_16s_F1F3<-subset(metadata_16s, Facility!="F2")
metadata_ITS_F1F3<-subset(metadata_ITS, Facility!="F2")

phyloseq16s_F1F3_RA<-phyloseq(otu_table(phyloseq16s_F1F3, taxa_are_rows = FALSE), TAX_16s, sample_data(metadata_16s_F1F3))
phyloseqITS_F1F3_RA<-phyloseq(otu_table(phyloseqITS_F1F3, taxa_are_rows = FALSE), TAX_ITS, sample_data(metadata_ITS_F1F3))

#prune taxa to remove OTUS that have less than 0.002 relative abundance in the OTU table
phyloseq16s_F1F3_RAprune<-prune_taxa(taxa_sums(phyloseq16s_F1F3_RA)>=0.002, phyloseq16s_F1F3_RA)
phyloseqITS_F1F3_RAprune<-prune_taxa(taxa_sums(phyloseqITS_F1F3_RA)>=0.002, phyloseqITS_F1F3_RA)

# Identify bacterial families associated with the presence of Salmonella in water samples.
model_lm_16s_F1F3 <- as.data.frame(cbind(otu_table(phyloseq16s_F1F3_RAprune), sample_data(phyloseq16s_F1F3_RAprune)[,5])) #Number in brackets is the column for metadata info with Lm presence
model_lm_16s_F1F3$L..monocytogenes <- as.factor(model_lm_16s_F1F3$L..monocytogenes)

model_lm_ITS_F1F3 <- as.data.frame(cbind(otu_table(phyloseqITS_F1F3_RAprune), sample_data(phyloseqITS_F1F3_RAprune)[,5])) #Number in brackets is the column for metadata info with Lm presence
model_lm_ITS_F1F3$L..monocytogenes <- as.factor(model_lm_ITS_F1F3$L..monocytogenes)

# Load Dataset
#16s
x_16s_F1F3 <- model_lm_16s_F1F3[,1:ncol(model_lm_16s_F1F3)-1] #OTU table
y_16s_F1F3 <- model_lm_16s_F1F3[,ncol(model_lm_16s_F1F3)] #Here, you can identify which column will be used as classifier. Lm presence column

dataset_rf_16s_F1F3 <- as.data.frame(cbind(x_16s_F1F3,y_16s_F1F3))

#ITS
x_ITS_F1F3 <- model_lm_ITS_F1F3[,1:ncol(model_lm_ITS_F1F3)-1] #OTU table
y_ITS_F1F3 <- model_lm_ITS_F1F3[,ncol(model_lm_ITS_F1F3)] #Here, you can identify which column will be used as classifier. Lm presence column

dataset_rf_ITS_F1F3 <- as.data.frame(cbind(x_ITS_F1F3,y_ITS_F1F3))


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
tunegrid_16s_F1F3 <- expand.grid(.mtry=c(40:90), .ntree=c(1000, 1500, 2000)) #mtry, ntree tuning parameters
#mtry sqroot of #of otu is mean, add 20 to each side to establish parameter
set.seed(1001)
custom_16s_F1F3 <- train(y_16s_F1F3~., data=dataset_rf_16s_F1F3, method=customRF, tuneGrid=tunegrid_16s_F1F3, trControl=control)
summary(custom_16s_F1F3)
plot(custom_16s_F1F3)

#ITS
tunegrid_ITS_F1F3 <- expand.grid(.mtry=c(35:75), .ntree=c(1000, 1500, 2000))
custom_ITS_F1F3 <- train(y_ITS_F1F3~., data=dataset_rf_ITS_F1F3, method=customRF, tuneGrid=tunegrid_ITS_F1F3, trControl=control)
summary(custom_ITS_F1F3)
plot(custom_ITS_F1F3)



#See the best tuned paramter mtry and ntree
custom_16s_F1F3$bestTune 
custom_ITS_F1F3$bestTune 
#Accuracy of the model - over 80% is a stable/ reliable model

#Get accuracy of the RF models
max(custom_16s_F1F3$results$Accuracy)
max(custom_ITS_F1F3$results$Accuracy)


#16s
rf_classifier_16s_F1F3 <- randomForest(y_16s_F1F3 ~ ., data= dataset_rf_16s_F1F3, ntree=1000, mtry=88, importance=TRUE) #mtry and mtree based on previous line
varimp_16s_F1F3 <- varImpPlot(rf_classifier_16s_F1F3) #shows variable importance - which OTU has the highest classification power
varimp_16s_F1F3 <- as.data.frame(varimp_16s_F1F3) 

f = varimp_16s_F1F3[order(varimp_16s_F1F3$MeanDecreaseGini, na.last=NA),] #purity of samples - the higher the better

row_varimp_16s_F1F3<-rownames(varimp_16s_F1F3)

#ITS
rf_classifier_ITS_F1F3 <- randomForest(y_ITS_F1F3 ~ ., data= dataset_rf_ITS_F1F3, ntree=1500, mtry=69, importance=TRUE) #mtry and mtree based on previous line
varimp_ITS_F1F3 <- varImpPlot(rf_classifier_ITS_F1F3) #shows variable importance - which OTU has the highest classification power
varimp_ITS_F1F3 <- as.data.frame(varimp_ITS_F1F3) 

f = varimp_ITS_F1F3[order(varimp_ITS_F1F3$MeanDecreaseGini, na.last=NA),] #purity of samples - the higher the better

row_varimp_ITS_F1F3<-rownames(varimp_ITS_F1F3)

#Replot accuracy and gini results
#16s
#Subset 10 OTUs with highest mean decrease accuracy
library(dplyr)
accuracy_16s_F1F3<-varimp_16s_F1F3[,1, drop=FALSE]

accuracy_16s_F1F3_top30<-accuracy_16s_F1F3 %>%
  top_n(30)

accuracy_16s_F1F3_top30$OTU<-rownames(accuracy_16s_F1F3_top30)  
accuracy_16s_F1F3_top30<-arrange(accuracy_16s_F1F3_top30, OTU)
accuracy_16s_F1F3_top30$Family<-subset(taxon_16s$Family, rownames(taxon_16s) %in% rownames(accuracy_16s_F1F3_top30))
accuracy_16s_F1F3_top30$Seq<-rep("Bacteria",30)


library(ggplot2)
accuracy_16s_F1F3_top30_plot<-ggplot(accuracy_16s_F1F3_top30, aes(x=MeanDecreaseAccuracy, y=reorder(OTU,MeanDecreaseAccuracy), color=Family))+
  geom_point(stat='identity', size=5)+geom_text(aes(label=Family, size=5), nudge_x = 0.7)+
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
gini_16s_F1F3_top30$Family<-subset(taxon_16s$Family, rownames(taxon_16s) %in% rownames(gini_16s_F1F3_top30))


gini_16s_F1F3_top30_plot<-ggplot(gini_16s_F1F3_top30, aes(x=MeanDecreaseGini, y=reorder(OTU,MeanDecreaseGini), color=Family))+
  geom_point(stat='identity')+
  geom_point(stat='identity', size=5)+geom_text(aes(label=Family, size=5), nudge_x = 0.7)+
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

#Combine plots
library(cowplot)
RandomForest_16s_F1F3 = plot_grid(accuracy_16s_F1F3_top30_plot, gini_16s_F1F3_top30_plot, 
                             ncol=2, nrow=1)
RandomForest_16s_F1F3
ggsave("RandomForest_16s_F1F3_top30.png", plot =RandomForest_16s_F1F3, device="png", width=16, height=8, units="in",dpi=600)

#ITS
#Subset 30 OTUs with highest mean decrease accuracy
accuracy_ITS_F1F3<-varimp_ITS_F1F3[,1, drop=FALSE]

accuracy_ITS_F1F3_top30<-accuracy_ITS_F1F3 %>%
  top_n(30)

accuracy_ITS_F1F3_top30$OTU<-rownames(accuracy_ITS_F1F3_top30)  
accuracy_ITS_F1F3_top30<-arrange(accuracy_ITS_F1F3_top30, OTU)
accuracy_ITS_F1F3_top30$Family<-subset(taxon_ITS$Family, rownames(taxon_ITS) %in% rownames(accuracy_ITS_F1F3_top30))
accuracy_ITS_F1F3_top30$Family<-gsub("f__","", accuracy_ITS_F1F3_top30$Family)
accuracy_ITS_F1F3_top30$Family<-gsub("o__","", accuracy_ITS_F1F3_top30$Family)
accuracy_ITS_F1F3_top30$Family<-gsub("p__","", accuracy_ITS_F1F3_top30$Family)
accuracy_ITS_F1F3_top30$Family<-gsub("c__","", accuracy_ITS_F1F3_top30$Family)
accuracy_ITS_F1F3_top30$Seq<-rep("Fungi",30)

accuracy_ITS_F1F3_top30_plot<-ggplot(accuracy_ITS_F1F3_top30, aes(x=MeanDecreaseAccuracy, y=reorder(OTU,MeanDecreaseAccuracy), color=Family))+
  geom_point(stat='identity', size=5)+geom_text(aes(label=Family, size=5), nudge_x = 0.5)+
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
gini_ITS_F1F3_top30$Family<-subset(taxon_ITS$Family, rownames(taxon_ITS) %in% rownames(gini_ITS_F1F3_top30))
gini_ITS_F1F3_top30$Family<-gsub("f__","", gini_ITS_F1F3_top30$Family)
gini_ITS_F1F3_top30$Family<-gsub("o__","", gini_ITS_F1F3_top30$Family)
gini_ITS_F1F3_top30$Family<-gsub("p__","", gini_ITS_F1F3_top30$Family)
gini_ITS_F1F3_top30$Family<-gsub("c__","", gini_ITS_F1F3_top30$Family)
gini_ITS_F1F3_top30$Seq<-rep("Fungi",30)


gini_ITS_F1F3_top30_plot<-ggplot(gini_ITS_F1F3_top30, aes(x=MeanDecreaseGini, y=reorder(OTU,MeanDecreaseGini), color=Family))+
  geom_point(stat='identity')+
  geom_point(stat='identity', size=5)+geom_text(aes(label=Family, size=5), nudge_x = 0.7)+
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




#outlier ranking by hclust
#Compared to dist whose input must be numeric variables, 
#the main feature of daisy is its ability to handle other variable types as well 
#(e.g. nominal, ordinal, (a)symmetric binary) even when different types occur in the same data set
install.packages("cluster")
install.packages("DMwR")
library(cluster)
library(DMwR)
dm<-daisy(admodule7)
o <- outliers.ranking(dm)

## Now let us check the outlier ranking factors ordered by decreasing score of outlyingness
o$prob.outliers[o$rank.outliers]
outlier<-cbind(o$prob.outliers,o$rank.outliers,amodule7)
write.table(outlier, "//data/ADI/Informatics/Informatics_Personnel/Ly Tran/ALE06/data/outlier.csv",sep=",",row.names=F)


outliers.ranking(data, test.data = NULL, method = "sizeDiff", method.pars = NULL, 
                 clus = list(dist = "euclidean",alg = "hclust", meth = "ward"), power = 1, verb = F)


install.packages("OutlierDetection")
OutlierDetection(iris[,-5])
#PCA outlier
PCOutlierDetection(iris[,-5])

#Distance-based outlier detection based on user-given neighborhood size
# Classify observations
cls_observations <- DB(dataset=X, d=1, fraction=0.05)$classification
# Remove outliers from dataset
X <- X[cls_observations=='Inlier',]


#find cluster:
#Unsupervised random forest_ find cluster:
library(readxl)
library(lattice)
require(caret)
library(cluster)
library(MASS)
library(randomForest)
library(fpc)
library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra) # clustering algorithms & visualization

all <- read_excel("//data/ADI/Informatics/Informatics_Personnel/Ly Tran/BLAST/data/module score/filtered_metadata.xlsx", sheet="all")
all <- read_excel("//data/ADI/Informatics/Informatics_Personnel/Ly Tran/BLAST/data/module score/filtered_metadata.xlsx", sheet="removeoutlier")
all <- read_excel("//data/ADI/Informatics/Informatics_Personnel/Ly Tran/BLAST/data/module score/filtered_metadata.xlsx", sheet="removeoutlier_hadldg12")
all <- read_excel("//data/ADI/Informatics/Informatics_Personnel/Ly Tran/BLAST/data/module score/filtered_metadata.xlsx", sheet="removeoutlier_hadldgab")
all <- read_excel("//data/ADI/Informatics/Informatics_Personnel/Ly Tran/BLAST/data/module score/filtered_metadata.xlsx", sheet="removeoutlier_noldg")

fm<-all[c("ID", "M1.1",	"M2.3",	"M3.1",	"M6.18",	"M1.2",	"M3.4",	"M5.12",	"M3.2",	"M4.2",	"M4.6",
          "M4.13",	"M5.1",	"M5.7",	"M7.1",	"M3.6",	"M8.46",	"M4.1",	"M4.15",	"M4.10",	"M4.11",	
          "M4.14",	"M5.15",	"LDG1.1",	"LDG2.1",	"M8.83",	"M6.6",	"M2.2",	"M3.3",	"M3.5",	"M4.7",
          "M6.11",	"M6.16",	"M9.42",	"M6.13",	"M4.3",	"M4.5",	"M5.9",	"M5.6",	"M5.10",	"M6.2",	"M6.12")]

fm<-all[c("ID", "M1.1",	"M2.3",	"M3.1",	"M6.18",	"M1.2",	"M3.4",	"M5.12",	"M3.2",	"M4.2",	"M4.6",
          "M4.13",	"M5.1",	"M5.7",	"M7.1",	"M3.6",	"M8.46",	"M4.1",	"M4.15",	"M4.10",	"M4.11",	
          "M4.14",	"M5.15",	"LDG_a",	"LDG_b",	"M8.83",	"M6.6",	"M2.2",	"M3.3",	"M3.5",	"M4.7",
          "M6.11",	"M6.16",	"M9.42",	"M6.13",	"M4.3",	"M4.5",	"M5.9",	"M5.6",	"M5.10",	"M6.2",	"M6.12")]

fm<-all[c("ID", "M1.1",	"M2.3",	"M3.1",	"M6.18",	"M1.2",	"M3.4",	"M5.12",	"M3.2",	"M4.2",	"M4.6",
          "M4.13",	"M5.1",	"M5.7",	"M7.1",	"M3.6",	"M8.46",	"M4.1",	"M4.15",	"M4.10",	"M4.11",	
          "M4.14",	"M5.15",	"LDG1.1",	"LDG2.1",	"LDG_a",	"LDG_b",	"M8.83",	"M6.6",	"M2.2",	"M3.3",	"M3.5",	"M4.7",
          "M6.11",	"M6.16",	"M9.42",	"M6.13",	"M4.3",	"M4.5",	"M5.9",	"M5.6",	"M5.10",	"M6.2",	"M6.12")]

fm<-all[c("ID", "M1.1",	"M2.3",	"M3.1",	"M6.18",	"M1.2",	"M3.4",	"M5.12",	"M3.2",	"M4.2",	"M4.6",
          "M4.13",	"M5.1",	"M5.7",	"M7.1",	"M3.6",	"M8.46",	"M4.1",	"M4.15",	"M4.10",	"M4.11",	
          "M4.14",	"M5.15",	"M8.83",	"M6.6",	"M2.2",	"M3.3",	"M3.5",	"M4.7",
          "M6.11",	"M6.16",	"M9.42",	"M6.13",	"M4.3",	"M4.5",	"M5.9",	"M5.6",	"M5.10",	"M6.2",	"M6.12")]

fm<-data.frame(fm, row.names=1)


set.seed(12345)
urfAll7a<- randomForest(x=fm, y=NULL,prox=T)

distance.matrix7a<-dist(1-urfAll7a$proximity)
mds.stuff7a <- cmdscale(distance.matrix7a, k=3)


k1 <- kmeans(mds.stuff7a, centers = 7)
centers <- k1$centers[k1$cluster, ]
distances <- sqrt(rowSums((mds.stuff7a - centers)^2))
outliers <- order(distances, decreasing=T)[1:20]
print(mds.stuff7a[outliers,])


#3-D 
mds.stuff7a=cbind(mds.stuff7a,k1[["cluster"]])
write.table(mds.stuff7a, "C:/Users/tranl/Desktop/Ly Tran/BLAST/data/mds.stuff7_191101.csv",sep=",",row.names=F)

#clustering data
fm_full=cbind(k1[["cluster"]],all$ID, fm)
write.table(fm_full, "C:/Users/tranl/Desktop/Ly Tran/BLAST/cluster_all.csv",sep=",",row.names=F)
write.table(fm_full, "C:/Users/tranl/Desktop/Ly Tran/BLAST/cluster_removeoulier.csv",sep=",",row.names=F)
write.table(fm_full, "C:/Users/tranl/Desktop/Ly Tran/BLAST/cluster_abremoveoulier.csv",sep=",",row.names=F)
write.table(fm_full, "C:/Users/tranl/Desktop/Ly Tran/BLAST/cluster_noldgremoveoulier.csv",sep=",",row.names=F)

###heatmap
amodule7sc<-read.csv(file="C:/Users/tranl/Desktop/Ly Tran/BLAST/cluster_all.csv", header=TRUE, sep=",")
amodule7sc<-read.csv(file="C:/Users/tranl/Desktop/Ly Tran/BLAST/cluster_removeoulier.csv", header=TRUE, sep=",")
amodule7sc<-read_excel("C:/Users/tranl/Desktop/Ly Tran/BLAST/cluster_abremoveoulier.xlsx")
amodule7sc<-read.csv(file="C:/Users/tranl/Desktop/Ly Tran/BLAST/cluster_noldgremoveoulier.csv", header=TRUE, sep=",")
amodule7sc<-read.csv(file="C:/Users/tranl/Desktop/Ly Tran/BLAST/cluster_removeoulierwithldgab.csv", header=TRUE, sep=",")

amodule7sc<-read_excel("//data/ADI/Informatics/Informatics_Personnel/Ly Tran/BLAST/data/module score/filtered_metadata.xlsx", sheet="BLAST")
amodule7sc<-read_excel("//data/ADI/Informatics/Informatics_Personnel/Ly Tran/BLAST/data/module score/filtered_metadata.xlsx", sheet="ABC")
amodule7sc<-read_excel("//data/ADI/Informatics/Informatics_Personnel/Ly Tran/BLAST/data/module score/filtered_metadata.xlsx", sheet="DxTerity")
amodule7sc<-read_excel("//data/ADI/Informatics/Informatics_Personnel/Ly Tran/BLAST/data/module score/filtered_metadata.xlsx", sheet="MESS")

admodule7sc<-amodule7sc[c("ID", "M1.1",	"M2.3",	"M3.1",	"M6.18",	"M1.2",	"M3.4",	"M5.12",	"M3.2",	"M4.2",	"M4.6",
                          "M4.13",	"M5.1",	"M5.7",	"M7.1",	"M3.6",	"M8.46",	"M4.1",	"M4.15",	"M4.10",	"M4.11",	
                          "M4.14",	"M5.15",	"LDG1.1",	"LDG2.1","LDG_a",	"LDG_b",	"M8.83",	"M6.6",	"M2.2",	"M3.3",	"M3.5",	"M4.7",
                          "M6.11",	"M6.16",	"M9.42",	"M6.13",	"M4.3",	"M4.5",	"M5.9",	"M5.6",	"M5.10",	"M6.2",	"M6.12")]

colnames(admodule7sc) <- c("ID","M1.1_Platelets",	"M2.3_Erythrocytes",	"M3.1_Erythrocytes",	
                           "M6.18_Erythrocytes",	"M1.2_Interferon",	"M3.4_Interferon",	"M5.12_Interferon",	
                           "M3.2_Inflammation",	"M4.2_Inflammation",	"M4.6_Inflammation",	"M4.13_Inflammation",	
                           "M5.1_Inflammation",	"M5.7_Inflammation",	"M7.1_Inflammation",	"M3.6_Cytotoxic/NK",
                           "M8.46_Cytotoxic/NK",	"M4.1_T cells",	"M4.15_T cells",	"M4.10_B cell",	
                           "M4.11_Plasma Cells",	"M4.14_Monocytes",	"M5.15_Neutrophils",	"LDG1.1",	"LDG2.1","LDG_a",	"LDG_b",
                           "M8.83_Immune Responses",	"M6.6_Apoptosis/Survival",	"M2.2_Cell Cycle",	
                           "M3.3_Cell Cycle",	"M3.5_Cell Cycle",	"M4.7_Cell Cycle",	"M6.11_Cell Cycle",	
                           "M6.16_Cell Cycle",	"M9.42_Cell Cycle",	"M6.13_Cell Death",	"M4.3_Protein Synthesis",	
                           "M4.5_Protein Synthesis",	"M5.9_Protein Synthesis",	"M5.6_Mitochondrial",	
                           "M5.10_Mitochondrial",	"M6.2_Mitochondrial",	"M6.12_Mitochondrial")


admodule7sc<-amodule7sc[c("ID", "M1.1",	"M2.3",	"M3.1",	"M6.18",	"M1.2",	"M3.4",	"M5.12",	"M3.2",	"M4.2",	"M4.6",
                          "M4.13",	"M5.1",	"M5.7",	"M7.1",	"M3.6",	"M8.46",	"M4.1",	"M4.15",	"M4.10",	"M4.11",	
                          "M4.14",	"M5.15",	"LDG1.1",	"LDG2.1",	"M8.83",	"M6.6",	"M2.2",	"M3.3",	"M3.5",	"M4.7",
                          "M6.11",	"M6.16",	"M9.42",	"M6.13",	"M4.3",	"M4.5",	"M5.9",	"M5.6",	"M5.10",	"M6.2",	"M6.12")]

colnames(admodule7sc) <- c("ID","M1.1_Platelets",	"M2.3_Erythrocytes",	"M3.1_Erythrocytes",	
                           "M6.18_Erythrocytes",	"M1.2_Interferon",	"M3.4_Interferon",	"M5.12_Interferon",	
                           "M3.2_Inflammation",	"M4.2_Inflammation",	"M4.6_Inflammation",	"M4.13_Inflammation",	
                           "M5.1_Inflammation",	"M5.7_Inflammation",	"M7.1_Inflammation",	"M3.6_Cytotoxic/NK",
                           "M8.46_Cytotoxic/NK",	"M4.1_T cells",	"M4.15_T cells",	"M4.10_B cell",	
                           "M4.11_Plasma Cells",	"M4.14_Monocytes",	"M5.15_Neutrophils",	"LDG1.1",	"LDG2.1",
                           "M8.83_Immune Responses",	"M6.6_Apoptosis/Survival",	"M2.2_Cell Cycle",	
                           "M3.3_Cell Cycle",	"M3.5_Cell Cycle",	"M4.7_Cell Cycle",	"M6.11_Cell Cycle",	
                           "M6.16_Cell Cycle",	"M9.42_Cell Cycle",	"M6.13_Cell Death",	"M4.3_Protein Synthesis",	
                           "M4.5_Protein Synthesis",	"M5.9_Protein Synthesis",	"M5.6_Mitochondrial",	
                           "M5.10_Mitochondrial",	"M6.2_Mitochondrial",	"M6.12_Mitochondrial")

admodule7sc<-amodule7sc[c("ID", "M1.1",	"M2.3",	"M3.1",	"M6.18",	"M1.2",	"M3.4",	"M5.12",	"M3.2",	"M4.2",	"M4.6",
                          "M4.13",	"M5.1",	"M5.7",	"M7.1",	"M3.6",	"M8.46",	"M4.1",	"M4.15",	"M4.10",	"M4.11",	
                          "M4.14",	"M5.15",	"LDG_a",	"LDG_b",	"M8.83",	"M6.6",	"M2.2",	"M3.3",	"M3.5",	"M4.7",
                          "M6.11",	"M6.16",	"M9.42",	"M6.13",	"M4.3",	"M4.5",	"M5.9",	"M5.6",	"M5.10",	"M6.2",	"M6.12")]

colnames(admodule7sc) <- c("ID","M1.1_Platelets",	"M2.3_Erythrocytes",	"M3.1_Erythrocytes",	
                           "M6.18_Erythrocytes",	"M1.2_Interferon",	"M3.4_Interferon",	"M5.12_Interferon",	
                           "M3.2_Inflammation",	"M4.2_Inflammation",	"M4.6_Inflammation",	"M4.13_Inflammation",	
                           "M5.1_Inflammation",	"M5.7_Inflammation",	"M7.1_Inflammation",	"M3.6_Cytotoxic/NK",
                           "M8.46_Cytotoxic/NK",	"M4.1_T cells",	"M4.15_T cells",	"M4.10_B cell",	
                           "M4.11_Plasma Cells",	"M4.14_Monocytes",	"M5.15_Neutrophils",	"LDG_a",	"LDG_b",
                           "M8.83_Immune Responses",	"M6.6_Apoptosis/Survival",	"M2.2_Cell Cycle",	
                           "M3.3_Cell Cycle",	"M3.5_Cell Cycle",	"M4.7_Cell Cycle",	"M6.11_Cell Cycle",	
                           "M6.16_Cell Cycle",	"M9.42_Cell Cycle",	"M6.13_Cell Death",	"M4.3_Protein Synthesis",	
                           "M4.5_Protein Synthesis",	"M5.9_Protein Synthesis",	"M5.6_Mitochondrial",	
                           "M5.10_Mitochondrial",	"M6.2_Mitochondrial",	"M6.12_Mitochondrial")

admodule7sc<-amodule7sc[c("ID", "M1.1",	"M2.3",	"M3.1",	"M6.18",	"M1.2",	"M3.4",	"M5.12",	"M3.2",	"M4.2",	"M4.6",
                          "M4.13",	"M5.1",	"M5.7",	"M7.1",	"M3.6",	"M8.46",	"M4.1",	"M4.15",	"M4.10",	"M4.11",	
                          "M4.14",	"M5.15",		"M8.83",	"M6.6",	"M2.2",	"M3.3",	"M3.5",	"M4.7",
                          "M6.11",	"M6.16",	"M9.42",	"M6.13",	"M4.3",	"M4.5",	"M5.9",	"M5.6",	"M5.10",	"M6.2",	"M6.12")]

colnames(admodule7sc) <- c("ID","M1.1_Platelets",	"M2.3_Erythrocytes",	"M3.1_Erythrocytes",	
                           "M6.18_Erythrocytes",	"M1.2_Interferon",	"M3.4_Interferon",	"M5.12_Interferon",	
                           "M3.2_Inflammation",	"M4.2_Inflammation",	"M4.6_Inflammation",	"M4.13_Inflammation",	
                           "M5.1_Inflammation",	"M5.7_Inflammation",	"M7.1_Inflammation",	"M3.6_Cytotoxic/NK",
                           "M8.46_Cytotoxic/NK",	"M4.1_T cells",	"M4.15_T cells",	"M4.10_B cell",	
                           "M4.11_Plasma Cells",	"M4.14_Monocytes",	"M5.15_Neutrophils",	
                           "M8.83_Immune Responses",	"M6.6_Apoptosis/Survival",	"M2.2_Cell Cycle",	
                           "M3.3_Cell Cycle",	"M3.5_Cell Cycle",	"M4.7_Cell Cycle",	"M6.11_Cell Cycle",	
                           "M6.16_Cell Cycle",	"M9.42_Cell Cycle",	"M6.13_Cell Death",	"M4.3_Protein Synthesis",	
                           "M4.5_Protein Synthesis",	"M5.9_Protein Synthesis",	"M5.6_Mitochondrial",	
                           "M5.10_Mitochondrial",	"M6.2_Mitochondrial",	"M6.12_Mitochondrial")

admodule7sc2 <- data.frame(admodule7sc, row.names=1)
admodulesct <- data.frame(t(admodule7sc2))

ano_col <- read_excel("//data/ADI/Informatics/Informatics_Personnel/Ly Tran/BLAST/data/module score/filtered_metadata.xlsx", sheet = "ano_col")
ano_col <- read_excel("//data/ADI/Informatics/Informatics_Personnel/Ly Tran/BLAST/data/module score/filtered_metadata.xlsx", sheet = "ano_colab")
ano_col <- read_excel("//data/ADI/Informatics/Informatics_Personnel/Ly Tran/BLAST/data/module score/filtered_metadata.xlsx", sheet = "ano_colnoldg")
ano_col <- read_excel("//data/ADI/Informatics/Informatics_Personnel/Ly Tran/BLAST/data/module score/filtered_metadata.xlsx", sheet = "ano_col12")

ano_col <- data.frame(ano_col, row.names=1)

col<-c("purple", "black", "yellow")
pal<-colorRampPalette(col)(100)
breaks<-seq(-3,3, length=100)

library("pheatmap")
pheatmap(admodulesct, 
         fontsize = 8,
         scale = "row",
         border_color = NA,
         cluster_rows = F,
         cluster_cols=F,
         breaks = breaks, show_colnames = F,
         angle_col = 315, color = pal, annotation_col=ano_col,
         gaps_row= c(1,4,7,14,20,25),
         gaps_col=c(231,371,634,932,1058,1190))

pheatmap(admodulesct, 
         fontsize = 8,
         scale = "row",
         border_color = NA,
         cluster_rows = F,
         cluster_cols=F,
         breaks = breaks, show_colnames = F,
         angle_col = 315, color = pal, annotation_col=ano_col,
         gaps_row= c(1,4,7,14,20,23),
         gaps_col=c(200,345,606,843,960,1113))

pheatmap(admodulesct, 
         fontsize = 8,
         scale = "row",
         border_color = NA,
         cluster_rows = F,
         cluster_cols=F,
         breaks = breaks, show_colnames = F,
         angle_col = 315, color = pal, 
         gaps_row= c(1,4,7,14,20,27),
         gaps_col=c(132,299,622,947,1071,1215))  

pheatmap(admodulesct, 
         fontsize = 8,
         scale = "row",
         border_color = NA,
         cluster_rows = F,
         cluster_cols=F,
         breaks = breaks, show_colnames = F,
         angle_col = 315, color = pal, 
         gaps_row= c(1,4,7,14,20,25),annotation_col=ano_col,
         gaps_col=c(243,384,635,958,1076,1255))

pheatmap(admodulesct, 
         fontsize = 8,
         scale = "row",
         border_color = NA,
         cluster_rows = F,
         cluster_cols=F,
         breaks = breaks, show_colnames = F,
         angle_col = 315, color = pal, annotation_col=ano_col,
         gaps_row= c(1,4,7,14,20,25),
         gaps_col=c(252,375,635,959,1087,1234))

pheatmap(admodulesct, 
         fontsize = 8,
         scale = "row",
         border_color = NA,
         cluster_rows = F,
         cluster_cols=F,
         main = "Heatmap for projects BLAST_removed outliers (266 samples)",
         breaks = breaks, show_colnames = F,
         angle_col = 315, color = pal, annotation_col=ano_col,
         gaps_row= c(1,4,7,14,20,25),
         gaps_col=c(40,88,121,168,189,218))

pheatmap(admodulesct, 
         fontsize = 8,
         scale = "row",
         border_color = NA,
         cluster_rows = F,
         cluster_cols=F,
         main = "Heatmap for projects ABC_removed outliers (232 samples)",
         breaks = breaks, show_colnames = F,
         angle_col = 315, color = pal, annotation_col=ano_col,
         gaps_row= c(1,4,7,14,20,25),
         gaps_col=c(23,59,98,139,160,188))

pheatmap(admodulesct, 
         fontsize = 8,
         scale = "row",
         border_color = NA,
         cluster_rows = F,
         cluster_cols=F,
         main = "Heatmap for projects DxTerity_removed outliers (495 samples)",
         breaks = breaks, show_colnames = F,
         angle_col = 315, color = pal, annotation_col=ano_col,
         gaps_row= c(1,4,7,14,20,25),
         gaps_col=c(42,139,229,342,380,427))

pheatmap(admodulesct, 
         fontsize = 8,
         scale = "row",
         border_color = NA,
         cluster_rows = F,
         cluster_cols=F,
         main = "Heatmap for projects MESS_removed outliers (197 samples)",
         breaks = breaks, show_colnames = F,
         angle_col = 315, color = pal, annotation_col=ano_col,
         gaps_row= c(1,4,7,14,20,25),
         gaps_col=c(14,50,92,133,154,167))

pheatmap(admodulesct, 
         fontsize = 8,
         scale = "row",
         border_color = NA,
         cluster_rows = F,
         cluster_cols=F,
         breaks = breaks, show_colnames = F,
         angle_col = 315, color = pal, annotation_col=ano_col,
         gaps_row= c(1,4,7,14,20,23),
         gaps_col=c(373,634,751,896,1049,1249))

##Test model
library(lattice)
require(caret)
library(cluster)
library(MASS)
library(randomForest)
library(fpc)
dmodule <- read.csv(file="C:/Users/tranl/Desktop/Ly Tran/BLAST/cluster.csv", header=TRUE, sep=",")

intrainc<-createDataPartition(dmodule$cluster,p=0.7, list=FALSE)
trainXc<-dmodule[intrainc,]
testXc<-dmodule[-intrainc,]

#find best mtry
set.seed(12345)
mtryc<-tuneRF(trainXc[,-1], factor(trainXc$cluster), ntreeTry=1000, stepFactor =1.5, improve=0.01, trace=TRUE, plot=TRUE)
bestmc<-mtryc[mtryc[,2]==min(mtryc[,2]),1]
print(mtryc)
print(bestmc)

#find best ntree
set.seed(12345)
modelRFc<-randomForest(factor(cluster)~., data=trainXc[,-1], mtry=12, ntree=90, importance=TRUE, prox=TRUE)
plot(modelRFc)

#final model
set.seed(12345)
modelRFc<-randomForest(factor(cluster)~., data=trainXc[,-1], mtry=12, ntree=90, importance=TRUE, prox=TRUE, keep.forest=TRUE)
print(modelRFc)
predc<-predict(modelRFc,testXc[,-1])
confusionMatrix(factor(predc), factor(testXc$cluster))

#post estimate_ROC curve
predc=predict(modelRFc,type = "prob")

library(ROCR)
perfc = prediction(predc, trainXc$cluster)
# 1. Area under curve
aucd = performance(perfd, "auc")
aucd
# 2. True Positive and Negative Rate
pred3 = performance(perf, "tpr","fpr")
# 3. Plot the ROC curve
plot(pred3,main="ROC Curve for Random Forest",col=2,lwd=2)
abline(a=0,b=1,lwd=2,lty=2,col="gray")

#predict class for test data set
pred<-predict(modelRF,testX, type="class")
table (predd,testXd$cluster)

varImpPlot(modelRFc, main=deparse(substitute("Determined modules and cytokines_7 clusters")))
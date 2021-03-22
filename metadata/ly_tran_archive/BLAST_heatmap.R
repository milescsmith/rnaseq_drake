library(readxl)
amodule7sc <- read_excel("//data/ADI/Informatics/Informatics_Personnel/Ly Tran/BLAST/data/BLAST_200312.xlsx", sheet = "SC")

admodule7sc<-amodule7sc[c("ID", "M1.1",	"M2.3",	"M3.1",	"M6.18",	"M1.2",	"M3.4",	"M5.12",	"M3.2",	"M4.2",	"M4.6",
                          "M4.13",	"M5.1",	"M5.7",	"M7.1",	"M3.6",	"M8.46",	"M4.1",	"M4.15",	"M4.10",	"M4.11",	
                          "M4.14",	"M5.15",	"ldg_1.1",	"ldg_2.1",	"M8.83",	"M6.6",	"M2.2",	"M3.3",	"M3.5",	"M4.7",
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

admodule7sc2 <- data.frame(admodule7sc, row.names=1)
admodulesct <- data.frame(t(admodule7sc2))

ano_col <- cbind(amodule7sc$id, amodule7sc$Cluster)
ano_col <- data.frame(ano_col, row.names=1)
colnames(ano_col)<-"Cluster"

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
         main = "Heatmap of determined module scores among BLAST patients_SC",
         breaks = breaks, show_colnames = T, angle_col=315,
         annotation_col=ano_col,
         color = pal,
         gaps_row= c(1,4,7,14,20,25),
         gaps_col=c(12))

amodule7sc <- read_excel("//data/ADI/Informatics/Informatics_Personnel/Ly Tran/BLAST/data/BLAST_200312.xlsx", sheet = "PN")
admodule7sc<-amodule7sc[c("ID", "M1.1",	"M2.3",	"M3.1",	"M6.18",	"M1.2",	"M3.4",	"M5.12",	"M3.2",	"M4.2",	"M4.6",
                          "M4.13",	"M5.1",	"M5.7",	"M7.1",	"M3.6",	"M8.46",	"M4.1",	"M4.15",	"M4.10",	"M4.11",	
                          "M4.14",	"M5.15",	"ldg_1.1",	"ldg_2.1",	"M8.83",	"M6.6",	"M2.2",	"M3.3",	"M3.5",	"M4.7",
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

admodule7sc2 <- data.frame(admodule7sc, row.names=1)
admodulesct <- data.frame(t(admodule7sc2))

pheatmap(admodulesct, 
         fontsize = 8,
         scale = "row",
         border_color = NA,
         cluster_rows = F,
         cluster_cols=F,annotation_col=ano_col,
         main = "Heatmap of determined module scores among BLAST patients_BL",
         breaks = breaks, 
         angle_col = 315, color = pal, 
         gaps_row= c(1,4,7,14, 20,25),show_colnames = T, 
         gaps_col=c(12))

amodule7sc <- read_excel("//data/ADI/Informatics/Informatics_Personnel/Ly Tran/BLAST/data/BLAST_200312.xlsx", sheet = "FV")
admodule7sc<-amodule7sc[c("ID", "M1.1",	"M2.3",	"M3.1",	"M6.18",	"M1.2",	"M3.4",	"M5.12",	"M3.2",	"M4.2",	"M4.6",
                          "M4.13",	"M5.1",	"M5.7",	"M7.1",	"M3.6",	"M8.46",	"M4.1",	"M4.15",	"M4.10",	"M4.11",	
                          "M4.14",	"M5.15",	"ldg_1.1",	"ldg_2.1",	"M8.83",	"M6.6",	"M2.2",	"M3.3",	"M3.5",	"M4.7",
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

admodule7sc2 <- data.frame(admodule7sc, row.names=1)
admodulesct <- data.frame(t(admodule7sc2))

pheatmap(admodulesct, 
         fontsize = 8,
         scale = "row",
         border_color = NA,
         cluster_rows = F,
         cluster_cols=F,annotation_col=ano_col,
         main = "Heatmap of determined module scores among BLAST patients_PN",
         breaks = breaks, 
         angle_col = 315, color = pal, 
         gaps_row= c(1,4,7,14, 20,25),show_colnames = T, 
         gaps_col=c(12))

amodule7sc <- read_excel("//data/ADI/Informatics/Informatics_Personnel/Ly Tran/BLAST/data/BLAST_200312.xlsx", sheet = "BL")
admodule7sc<-amodule7sc[c("ID", "M1.1",	"M2.3",	"M3.1",	"M6.18",	"M1.2",	"M3.4",	"M5.12",	"M3.2",	"M4.2",	"M4.6",
                          "M4.13",	"M5.1",	"M5.7",	"M7.1",	"M3.6",	"M8.46",	"M4.1",	"M4.15",	"M4.10",	"M4.11",	
                          "M4.14",	"M5.15",	"ldg_1.1",	"ldg_2.1",	"M8.83",	"M6.6",	"M2.2",	"M3.3",	"M3.5",	"M4.7",
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

admodule7sc2 <- data.frame(admodule7sc, row.names=1)
admodulesct <- data.frame(t(admodule7sc2))

pheatmap(admodulesct, 
         fontsize = 8,
         scale = "row",
         border_color = NA,
         cluster_rows = F,
         cluster_cols=F,annotation_col=ano_col,
         main = "Heatmap of determined module scores among BLAST patients_FV",
         breaks = breaks, 
         angle_col = 315, color = pal, 
         gaps_row= c(1,4,7,14, 20,25),show_colnames = T, 
         gaps_col=c(12))


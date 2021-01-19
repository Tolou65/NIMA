install.packages("devtools")
library(devtools)

devtools::install_github("DrugComb/TidyComb")


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

library("TidyComb")
library("dplyr")
library("xlsx")
setwd("~/drug_combo")
inhib_target_ann_tab <- read.xlsx("~/Desktop/drug_comb/Table_Inhibitor-target-profile-annotation.xlsx", 1,header=TRUE, colClasses=NA)
inhib_target_ann_tab <- read.xlsx("~/drug_combo/data2/Experiment_and_scores_list_summaryFIMM.xlsx", 1,header=TRUE, colClasses=NA)
CIDs <- inhib_target_ann_tab$PubChem.CID
CIDs <- CIDs[!is.na(CIDs)]
CIDs_unique <- as.character(unique(CIDs))
CIDs_unique_char <- replace(CIDs_unique, c(8), "1271002")
PubchemPro_table <- GetPubchemPro(CIDs_unique_char)
PubchemPro_table <- rbind(PubchemPro_table, GetPubchemPro(CIDs_unique_char[c(9,10,11,12)]))
PubchemPro_table <- cbind(PubchemPro_table, GetPubNames(CIDs_unique_char)[2])
PubchemPro_table <- cbind(PubchemPro_table, GetChembl(PubchemPro_table$inchikey)[c(2,3)])
PubchemPro_table <- cbind(PubchemPro_table, GetIds(PubchemPro_table$inchikey)[c(2,3)])
##For study table:
study_table <- GetPubmed(Drugs_list, tool = "TidyComb", email= "tolou.shadbaher@gmail.com")





##changing the column names accordingly
names(PubchemPro_table)[names(PubchemPro_table)=="uni_drugbank"] <- "drugbank_id"
names(PubchemPro_table)[names(PubchemPro_table)=="uni_kegg_c"] <- "kegg_id"
names(PubchemPro_table)[names(PubchemPro_table)=="chembl_phase"] <- "clinical_phase"
names(PubchemPro_table)[names(PubchemPro_table)=="name"] <- "dname"
###for changing the order of columns my_data2 <- my_data[, c(5, 4, 1, 2, 3)]




####new drugs
#PubChem
new_drug <- CheckDrug(PubchemPro_table$cid)$new
drugs_all <- CheckDrug(PubchemPro_table$cid)
drug_name <- GetPubNames(cids = new_drug)
drug_info <- GetPubchemPro(cids = new_drug)
drug <- full_join(drug_name, drug_info)
pub_phase <- GetPubPhase(drug$cid)
drug <- left_join(drug, pub_phase, by = "cid")
#ChEMBL
chembl <- GetChembl(drug$inchikey)
drug <- left_join(drug, chembl, by = "inchikey")
#UniChem
unichem <- GetIds(drug$inchikey)
drug <- left_join(drug, unichem, by = "inchikey") 
drug$chembl_id <- unlist(apply(drug[, c("chembl_id.x", "chembl_id.y")], 1,
                               function(x){
                                 if (length(na.omit(x)) == 0){
                                   return(NA)
                                 } else{
                                   return(paste(na.omit(x), collapse = "; "))
                                 }
                               }))
drug <- drug[, which(!colnames(drug) %in% c("chembl_id.x", "chembl_id.y"))]
##assigning the new id for the new drugs
drug$id <- seq(new_drug$n + 1, length.out = nrow(drug))


### drug bank file can be read by using read.table library as follow:
library("data.table")
system.time(db_mod <- fread("~/drug_combo/Data/drug_links.csv"))


### STITCH database reading
#library("ff")
#chemicall <- read.csv.ffdf(file="~/drug_combo/chemicals.v5.0.tsv", header=TRUE, sep="\t")
library(sqldf)
f <- file("~/drug_combo/chemicals.v5.0.tsv")
system.time(chemical <- sqldf("select * from f", dbname = tempfile(), file.format = list(header = T, row.names = F)))
## correction for the columns
drug <- select(drug, dname = "name", id, chembl_id, inchikey, smiles, cid, 
               molecular_formula, clinical_phase = "phase", #cid_m, cid_s, stitch_name, 
               drugbank_id = "uni_drugbank", kegg_id = "uni_kegg_c")



##cellosaurus dowloading


# Download "cellosaurus.html" file
download.file("ftp://ftp.expasy.org/databases/cellosaurus/cellosaurus.xml", "cellosaurus.xml")

cell_repl1 <- read.csv("~/drug_combo/Data/Primary_screen/Raw_data/NTNU_HTS_Repl1_raw.csv", header=TRUE, sep = "\t")
cell_repl2 <- read.csv("~/drug_combo/Data/Primary_screen/Raw_data/NTNU_HTS_Repl2_raw.csv", header=TRUE, sep="\t")
cell_repl3 <- read.csv("~/drug_combo/Data/Primary_screen/Raw_data/NTNU_HTS_Repl3_raw.csv", header=TRUE, sep="\t")

cells_1 <- as.character(unique(na.omit(cell_repl1$Cellline)))# The cell name vector must be unique and without NA
cells_2 <- as.character(unique(na.omit(cell_repl2$Cellline)))
cells_3 <- as.character(unique(na.omit(cell_repl3$Cellline)))
All_CellLines <- unique(c(cells_1,cells_2,cells_3))
cells.match <- MatchCellAcc(All_CellLines, file = file.path("~/Desktop/drug_comb/cellosaurus.xml"))
#All_cel <- c('SW872', '93T449')
#cells.match <- MatchCellAcc(All_cel, file = file.path("/u/63/shadbat1/unix/drug_combo/Data/cellosaurus.xml"))
print(cells.match)
cell.match.clean <- cells.match[-9, c("input_name", "cellosaurus_accession")]
#cell1.match <- MatchCellAcc(cells_1, file = system.file("~/drug_combo/cellosaurus.xml", package = "TidyComb")) """not working"""


## generating table according to cellosaurus accession table:

cell <- GenerateCell(cell.match.clean$cellosaurus_accession, file = file.path("~/Desktop/drug_comb/cellosaurus.xml"))
cell$cell_line
cell_index <- cell$cell_id %>% 
  left_join(cell.match.clean, by = "cellosaurus_accession")



## constructing tissue table
new_tisse <- CheckTissue(cell$cell_line$tissue)

############################################################################################
#########Calculation of synergy score from Drug Combo with TidyComb for new set:############
############################################################################################

ERI_GemCI <- read.csv("~/drug_combo/data2/ERI_GEMCI.csv", header = TRUE, row.names = 1)
names(ERI_GemCI) <- c('100','33.33', '11.11','3.70','1.23','0')
ERI_GEM <- ERI_GemCI[, c(6,5,4,3,2,1)]
m1<- 100-matrix(unlist(ERI_GEM), ncol=6, byrow=FALSE)
rownames(m1) <- c(row.names(ERI_GemCI))
colnames(m1) <- c(names(ERI_GEM))
ERI_GEMCI_Syng_93T449 <- CalculateMat(m1, summary.only = TRUE)


ERI_PALBO <- read.csv("~/drug_combo/data2/ERI_PALBO.csv", header = TRUE, row.names = 1)
names(ERI_PALBO) <- c('100','33.3333', '11.1111','3.7037','1.2345','0')
ERI_PAL <- ERI_PALBO[, c(6,5,4,3,2,1)]
m2<- 100-matrix(unlist(ERI_PAL), ncol=6, byrow=FALSE)
rownames(m2) <- c(row.names(ERI_PALBO))
colnames(m2) <- c(names(ERI_PAL))
ERI_PALBO_Syng_93T449 <- CalculateMat(m2, summary.only = TRUE)


ERI_PAZO <- read.csv("~/drug_combo/data2/ERI_PAZO.csv", header = TRUE, row.names = 1)
names(ERI_PAZO) <- c('100','33.3333', '11.1111','3.7037','1.2345','0')
ERI_PAZ <- ERI_PAZO[, c(6,5,4,3,2,1)]
m3<- 100-matrix(unlist(ERI_PAZ), ncol=6, byrow=FALSE)
rownames(m3) <- c(row.names(ERI_PAZO))
colnames(m3) <- c(names(ERI_PAZ))
ERI_PAZO_Syng_93T449 <- CalculateMat(m3, summary.only = TRUE)


ERI_DACAR <- read.csv("~/drug_combo/data2/ERI_DACAR.csv", header = TRUE, row.names = 1)
names(ERI_DACAR) <- c('10','3.3333', '1.1111','0.37','0.12','0')
ERI_DAC <- ERI_DACAR[, c(6,5,4,3,2,1)]
m4<- 100-matrix(unlist(ERI_DAC), ncol=6, byrow=FALSE)
rownames(m4) <- c(row.names(ERI_DACAR))
colnames(m4) <- c(names(ERI_DAC))
ERI_DACAR_Syng_sw872 <- CalculateMat(m4, summary.only = TRUE)

ERI_GEMCI_2 <- read.csv("~/drug_combo/data2/ERI_GEMCI_SW.csv", header = TRUE, row.names = 1)
names(ERI_GEMCI_2) <- c('10','5', '2.5','1.25','0.625','0')
ERI_GEMCI_2 <- ERI_GEMCI_2[, c(6,5,4,3,2,1)]
m5 <- 100-matrix(unlist(ERI_GEMCI_2), ncol=6, byrow=FALSE)
rownames(m5) <- c(row.names(ERI_GEMCI_2))
colnames(m5) <- c(names(ERI_GEMCI_2))
ERI_GEMCI_Syng_sw872 <- CalculateMat(m5, summary.only = TRUE)


ERI_GEMCI_2_2 <- read.csv("~/drug_combo/data2/ERI_GEMCI_SW_2.csv", header = TRUE, row.names = 1)
names(ERI_GEMCI_2_2) <- c('10','3.33', '1.11','0.37','0.12','0')
ERI_GEMCI_2_2 <- ERI_GEMCI_2_2[, c(6,5,4,3,2,1)]
m6 <- 100-matrix(unlist(ERI_GEMCI_2_2), ncol=6, byrow=FALSE)
rownames(m6) <- c(row.names(ERI_GEMCI_2_2))
colnames(m6) <- c(names(ERI_GEMCI_2_2))
ERI_GEMCI_Syng_sw872_190320 <- CalculateMat(m6, summary.only = TRUE)




ERI_IFOS <- read.csv("~/drug_combo/data2/ERI_IFOS.csv", header = TRUE, row.names = 1)
names(ERI_IFOS) <- c('10','3.33', '1.11','0.37','0.12','0')
ERI_IFOS <- ERI_IFOS[, c(6,5,4,3,2,1)]
m7<- 100-matrix(unlist(ERI_IFOS), ncol=6, byrow=FALSE)
rownames(m7) <- c(row.names(ERI_IFOS))
colnames(m7) <- c(names(ERI_IFOS))
ERI_IFOS_Syng_SW872 <- CalculateMat(m7, summary.only = TRUE)

ERI_PALBO_2 <- read.csv("~/drug_combo/data2/ERI_PALBO_SW.csv", header = TRUE, row.names = 1)
names(ERI_PALBO_2) <- c('10','3.3333', '1.1111','0.37','0.12','0')
ERI_PALBO_2 <- ERI_PALBO_2[, c(6,5,4,3,2,1)]
m8<- 100-matrix(unlist(ERI_PALBO_2), ncol=6, byrow=FALSE)
rownames(m8) <- c(row.names(ERI_PALBO_2))
colnames(m8) <- c(names(ERI_PALBO_2))
ERI_PALBO_Syng_sw872 <- CalculateMat(m8, summary.only = TRUE)


ERI_PAZO_2 <- read.csv("~/drug_combo/data2/ERI_PAZO_SW.csv", header = TRUE, row.names = 1)
names(ERI_PAZO_2) <- c('10','3.3333', '1.1111','0.37','0.12','0')
ERI_PAZO_2 <- ERI_PAZO_2[, c(6,5,4,3,2,1)]
m9<- 100-matrix(unlist(ERI_PAZO_2), ncol=6, byrow=FALSE)
rownames(m9) <- c(row.names(ERI_PAZO_2))
colnames(m9) <- c(names(ERI_PAZO_2))
ERI_PAZO_Syng_sw872 <- CalculateMat(m9, summary.only = TRUE)


ERI_PAZO_2_2 <- read.csv("~/drug_combo/data2/ER_PAZO_SW_2.csv", header = TRUE, row.names = 1)
names(ERI_PAZO_2_2) <- c('10','3.3333', '1.1111','0.37','0.12','0')
ERI_PAZO_2_2 <- ERI_PAZO_2[, c(6,5,4,3,2,1)]
m10<- 100-matrix(unlist(ERI_PAZO_2_2), ncol=6, byrow=FALSE)
rownames(m10) <- c(row.names(ERI_PAZO_2_2))
colnames(m10) <- c(names(ERI_PAZO_2_2))
ERI_PAZO_Syng_sw872_211019 <- CalculateMat(m9, summary.only = TRUE)



drugComb_result <- rbind(ERI_GEMCI_Syng_93T449, ERI_PALBO_Syng_93T449, ERI_PAZO_Syng_93T449, ERI_DACAR_Syng_sw872, ERI_GEMCI_Syng_sw872, ERI_GEMCI_Syng_sw872_190320, 
                         ERI_IFOS_Syng_SW872, ERI_PALBO_Syng_sw872, ERI_PAZO_Syng_sw872, ERI_PAZO_Syng_sw872_211019)

rownames(drugComb_result) <- c('ERI_GEMCI_93T449', 'ERI_PALBO_93T449', 'ERI_PAZO_93T449', 'ERI_DACAR_sw872', 'ERI_GEMCI_sw872', 'ERI_GEMCI_sw872_190320', 'ERI_IFOS_SW872'
                               ,'ERI_PALBO_sw872', 'ERI_PAZO_sw872', 'ERI_PAZO_sw872_211019')


write.csv(drugComb_result, "/u/63/shadbat1/unix/drug_combo/data2/DrugComb_Syng_results.csv")

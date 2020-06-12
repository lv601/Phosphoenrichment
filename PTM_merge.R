library(readr)
library(tidyverse)
library(Hmisc)
library(data.table)
library(dplyr)
library( plyr )
#install.packages("doParallel")
library(doParallel)
library(foreach)
library(readr)
library(dpe)
library(Hmisc)
library(fs)
options(digits = 6)
require(data.table)
setDTthreads(threads = 0)


#6 Replikas reinladen:
#DP_C_CYT_1_1 <- read_csv("//dc/Users/private_mg/Documents/Projekt_PTM_merge/results_PTM_melanoma/results/DP-C_CYT_1_1.mgf.ionbot.csv")
#DP_C_CYT_1_2 <- read_csv("//dc/Users/private_mg/Documents/Projekt_PTM_merge/results_PTM_melanoma/results/DP-C_CYT_1_2.mgf.ionbot.csv")
#DP_C_CYT_2_1 <- read_csv("//dc/Users/private_mg/Documents/Projekt_PTM_merge/results_PTM_melanoma/results/DP-C_CYT_2_1.mgf.ionbot.csv")
#DP_C_CYT_2_2 <- read_csv("//dc/Users/private_mg/Documents/Projekt_PTM_merge/results_PTM_melanoma/results/DP-C_CYT_2_2.mgf.ionbot.csv")
#DP_C_CYT_3_1 <- read_csv("//dc/Users/private_mg/Documents/Projekt_PTM_merge/results_PTM_melanoma/results/DP-C_CYT_3_1.mgf.ionbot.csv")
#DP_C_CYT_3_2 <- read_csv("//dc/Users/private_mg/Documents/Projekt_PTM_merge/results_PTM_melanoma/results/DP-C_CYT_3_2.mgf.ionbot.csv")
###Fragestellung: Welche Unterschiede habe ich Zwischen Fraktionen innerhalb einer Zelllinie: DP-C-Cyt und DP-C-NE
###               Welche Unterschiede habe ich zwischen Zuständen: DP-C-Cyt und DP-CB-Cyt; DP-C-NE und DP-CB-NE
###Finde vorher high-confidence peptides: q_value <= 0.05; DB != D; eventuell weitere Kriterien;
###Finde vorher high-confidence peptides: matched_peptide in allen 6 Replikaten (3 Biologische jeweils 2x technisch) AND (precursor_mass + z) ident
###Dann Datenpunkt behalten
###Unterschiede zwischen Zelllinien, zwischen C und CB (also metastasierungszuständen) und Fraktionen (NE,Cyt,SN)
##Loop over csv; parse info from filename to column
#variable names:
#[1] "scan_id"                 "charge"                  "precursor_mass"          "matched_peptide"         "modifications"          
#[6] "ionbot_psm_score"        "DB"                      "unexpected_modification" "ms2pip_pearsonr"         "proteins"               
#[11] "num_unique_pep_ids"      "percolator_psm_score"    "q_value"                 "PEP"                     "title"

#data_path = "//DC/Users/private_mg/Documents/Projekt_PTM_merge/results_PTM_melanoma/results/"
data_path = "C:/Users/601/Documents/Projekt_PTM_merge/results_PTM_melanoma/results"
setwd(data_path)
#create a list of the files from your target directory
file_list <- list.files(path=data_path)
#initiate a blank data frame, each iteration of the loop will append the data from the given file to this variable
dataset <- data.frame()

###Loop for separate dataframes
#data_path = "//DC/Users/private_mg/Documents/Projekt_PTM_merge/results_PTM_melanoma/results/"
file_list <- list.files(path=data_path, pattern = "*.csv", include.dirs = FALSE)
#file_list = gsub(".mgf.ionbot.csv","",file_list)

table_list = gsub(".mgf.ionbot.csv","",file_list); table_list
state_list = gsub("_[0-9]_[0-9]","",table_list);state_list=unique(state_list);state_list
list_of_tables = data.frame()
a = 0
for(i in 1:length(file_list)) {
  temp_data <- fread(file_list[i], stringsAsFactors = F, nThread = 12, header = TRUE, sep = ",",sep2 =";"
                     ,data.table = TRUE) #read in files using the fread function from the data.table package
  temp_data = temp_data[which(temp_data$q_value<=0.05), ]
  temp_data = temp_data[which(temp_data$DB != "D"), ]
  temp_data$title = sub("File PROT-DATAMS StoragePRIDEPRIDE-Melanoma_WitzRAW"
                        , path_sanitize(temp_data$title), replacement = "",fixed=TRUE)
  fraction_state = strsplit(file_list[i],"_[0-9]")[[1]][1]
  print(fraction_state)
  
  ##Wenn neuer zustand beginnt
  if(a == 0 || fraction_state != state_list[a]) {
    print("nah")
    a = a+1
    
  }
  else {
    assign(paste(state_list[a],sep=""),setDT(rbindlist(list(dataset, temp_data), use.names = T, fill=TRUE)))
    #assign(paste("DT_",state_list[a],sep=""),setDT(rbindlist(list(paste("DT_",state_list[a],sep=""), temp_data), use.names = T, fill=TRUE)))
    #DTT = lapply(paste(state_list[a],sep=""), data.table)
  }
    
  dataset <- setDT(rbindlist(list(dataset, temp_data), use.names = T, fill=TRUE)); setkey(dataset, "matched_peptide") #for each iteration, bind the new data to the building dataset
  assign(paste(gsub(".mgf.ionbot.csv","",file_list[i]),sep="_"),rbind(temp_data, use.names = T, fill=TRUE))
  list_of_tables = c(list_of_tables,gsub(".mgf.ionbot.csv","",file_list[i]))
  #temp_data = ""
  #paste(gsub(".mgf.ionbot.csv","",file_list[i])) = as.data.frame(paste(gsub(".mgf.ionbot.csv","",file_list[i])))
}






#Create list of matched_peptides occuring in all 6 Replikas

##Test full out join multiple key columns
# make main table with key col(s)
kc = c("modifications","matched_peptide","precursor_mass")
DT = setkey(unique(rbindlist(lapply(DTs, `[`, j = ..kc))))
# get non-key cols
#kc = "modifications"
DTs = lapply(list(`DP-C_CYT_1_1`,`DP-C_CYT_1_2`,`DP-C_CYT_2_1`,`DP-C_CYT_2_2`
                  ,`DP-C_CYT_3_1`,`DP-C_CYT_3_2`), data.table, key=c("modifications","matched_peptide","precursor_mass"))
DTs = lapply()
DT = unique(rbindlist(lapply(DTs, `[`, j = ..kc)))
DT = setkey(unique(rbindlist(lapply(DTs, `[`, j = ..kc))))
for (d in DTs){
  cols = setdiff(names(d), kc)
  DT[d, (cols) := mget(sprintf("i.%s", cols)),on=c("modifications","matched_peptide","precursor_mass") ][]
}
##apply rounding on precursor mass and use column too merge
##what about charge
DT$precursor_mass = round(DT$precursor_mass,digits=1)
d$precursor_mass = round(d$precursor_mass,digits=1)
setkey(d,modifications,matched_peptide,precursor_mass,charge,unexpected_modification)
d = unique(d,by=key(d))

####Old code
DTs = lapply(list(`DP-C_CYT_1_1`,`DP-C_CYT_1_2`,`DP-C_CYT_2_1`,`DP-C_CYT_2_2`
                  ,`DP-C_CYT_3_1`,`DP-C_CYT_3_2`), data.table)
testapply = lapply(DTs, `[`, j = c("scan_id","modifications"))

Id_keep = rbindlist(lapply(DTs, `[`, j = "matched_peptide"))[, .N, by=matched_peptide][N >= 6L, matched_peptide];Id_keep
#hier ein j=lapply matched_pep+mod+prec
#Id_keep = rbindlist(lapply(DTs, `[`, j = "matched_peptide"))[, .N, by=matched_peptide][N >= 6L, matched_peptide];Id_keep
#multiple columns
"%IN%" <- function(x, y) interaction(x) %in% interaction(y)
`DP-C_CYT_1_1` %IN% `DP-C_CYT_1_2`


DT_id_keep = Reduce(funion, DTs)[matched_peptide %in% Id_keep];DT_id_keep

#Create modifications_keep; merge id_keep and mod_keep
mod_keep = rbindlist(lapply(DTs, `[`, j = "modifications"))[, .N, by=modifications][N >= 6L, modifications];mod_keep
DT_mod_keep = Reduce(funion, DTs)[modifications %in% mod_keep];DT_mod_keep

merge_keep = merge(DT_id_keep,DT_mod_keep,by=c("matched_peptide","modifications","unexpected_modification","proteins","charge","scan_id"), allow.cartesian = T)
merge_distinct = merge_keep[!duplicated(merge_keep$scan_id), ]

mergeres = merge(`DP-C_CYT_1_1`,`DP-C_CYT_1_2`, by = c("matched_peptide","modifications","proteins","charge"), suffixes = c(table_list[1],table_list[2]), allow.cartesian = TRUE)
merger2 = merge(mergeres, `DP-C_CYT_2_1`,by = c("matched_peptide","modifications","proteins","charge"), suffixes = c(table_list[2],table_list[3]), allow.cartesian = TRUE)
merger3 = merge(merger2,`DP-C_CYT_2_2`,by = c("matched_peptide","modifications","proteins","charge"), suffixes = c(table_list[3],table_list[4]),allow.cartesian = TRUE)
merger4 = merge(merger3,`DP-C_CYT_3_1`,by = c("matched_peptide","modifications","proteins","charge"), suffixes = c(table_list[4],table_list[5]),allow.cartesian = TRUE)
merger5 = merge(merger4,`DP-C_CYT_3_2`,by = c("matched_peptide","modifications","proteins","charge"), suffixes = c(table_list[5],table_list[6]), allow.cartesian = TRUE)



Merge(DT_keep)




########################################################################
#my_data %>% distinct(Sepal.Length, Petal.Width, .keep_all = TRUE)
###Merging
a = setDT(`d_DP-C_CYT_1_1.mgf.ionbot.csv`,key="matched_peptide")
b = setDT(`d_DP-C_CYT_1_2.mgf.ionbot.csv`,key="matched_peptide")
c = setDT(`d_DP-C_CYT_2_1.mgf.ionbot.csv`); setkey(c, "matched_peptide")
d = setDT(`d_DP-C_CYT_2_2.mgf.ionbot.csv`); setkey(d, "matched_peptide")
setkey(a, "Modified sequence")
setkey(b, "Modified sequence")
setkey(c, "Modified sequence")
setkey(d, "Modified sequence")

#character data table
anew = a[, lapply(.SD, as.character)]
bnew = b[, lapply(.SD, as.character)]

bnew[bnew %in% anew]

merger = merge(a,b, by=c("Sequence","Modifications","Modified Sequence","Mass"))
merger = merger[, lapply(.SD, as.character)]
library(reshape2)

Merge(a,b, id = "title", by=c("matched_peptide","precursor_mass","modifications"))

Merge(a,b, all=FALSE, id=~matched_peptide)

#Merge Master dataset with single tables
Merge(dataset,`d_DP-C_CYT_1_1.mgf.ionbot.csv`,`d_DP-C_CYT_1_2.mgf.ionbot.csv`,`d_DP-C_CYT_2_1.mgf.ionbot.csv`,`d_DP-C_CYT_2_2.mgf.ionbot.csv`,all=FALSE)

######
#had to specify columns to get rid of the total column
for (i in 1:length(file_list)){
  temp_data <- fread(file_list[i], stringsAsFactors = F) #read in files using the fread function from the data.table package
  #temp_data2 = temp_data[,2:16]
  dataset <- setDT(rbindlist(list(dataset, temp_data), use.names = T, fill=TRUE)); setkey(dataset, "matched_peptide") #for each iteration, bind the new data to the building dataset
}

DP_C_CYT_1_1 <- read_csv("//dc/Users/private_mg/Documents/Projekt_PTM_merge/results_PTM_melanoma/results/DP-C_CYT_1_1.mgf.ionbot.csv", col_types = cols(title = col_skip()), locale = locale())
df = select(data, scan_id:unexpected_modification)
x = df
x[!duplicated(x$matched_peptide), ] #negate duplicates - remove duplicates by column matched_peptide
df %>% distinct(df$matched_peptide, .keep_all = TRUE)
data_filtered = data[which(data$q_value<=0.05), ]

data_filtered2 = data_filtered[which(data_filtered$DB != "D"), ]
colnms=c("charge","precursor_mass")
data_filtered2$new_col<-rowSums(mtcars[,colnms])
h_mass = 1.00784

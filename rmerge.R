DP_C_Cyt_SN_1_1 <- read_excel("~/test_data/DP_C_Cyt_SN_1_1.xlsx")
DP_C_Cyt_SN_1_2 <- read_excel("~/test_data/DP_C_Cyt_SN_1_2.xlsx")
DP_C_Cyt_SN_1_3 <- read_excel("~/test_data/DP_C_Cyt_SN_1_3.xlsx")
DP_C_Cyt_SN_2_1 <- read_excel("~/test_data/DP_C_Cyt_SN_2_1.xlsx")
DP_C_Cyt_SN_2_2 <- read_excel("~/test_data/DP_C_Cyt_SN_2_2.xlsx")
DP_C_Cyt_SN_2_3 <- read_excel("~/test_data/DP_C_Cyt_SN_2_3.xlsx")


df = merge(DP_C_Cyt_SN_1_1, DP_C_Cyt_SN_1_2, by.x = c("Modified sequence"), by.y = c("Modified sequence"), all=TRUE)
write.table(df1, file="C:/Users/601/Documents/test_data/a.csv")
df_match_1_2 = match(DP_C_Cyt_SN_1_1,DP_C_Cyt_SN_1_2)

install.packages("Hmisc")
library(Hmisc)
library(data.table)

a <- data.frame(sid=1:3, age=c(20,30,40))
b <- data.frame(sid=c(1,2,2), bp=c(120,130,140))
d <- data.frame(sid=c(1,3,4), wt=c(170,180,190))
all <- Merge(a, b, d, id = ~ sid)
all

all = Merge(DP_C_Cyt_SN_1_1,DP_C_Cyt_SN_1_2, id = ~ "Mass")

a <- data.table(a); setkey(a, sid)
# data.table also does not allow duplicates without allow.cartesian=TRUE
b <- data.table(sid=1:2, bp=c(120,130)); setkey(b, sid)
d <- data.table(d); setkey(d, sid)
all <- Merge(a, b, d)

a = setDT(DP_C_Cyt_SN_1_1)
b = setDT(DP_C_Cyt_SN_1_2)
c = setDT(DP_C_Cyt_SN_1_3)
setkey(a, "Modified sequence")
setkey(b, "Modified sequence")
setkey(c, "Modified sequence")

#character data table
anew = a[, lapply(.SD, as.character)]
bnew = b[, lapply(.SD, as.character)]

bnew[bnew %in% anew]

merger = merge(a,b, by=c("Sequence","Modifications","Modified Sequence","Mass"))
merger = merger[, lapply(.SD, as.character)]
library(reshape2)
#testen

df1 <- data.frame(Id=c(1,2,3,4),a=c(0,1,0,2),b=c(1,0,1,0),c=c(0,0,4,0),d=c(0,0,4,0)) 
df2 <- data.frame(Id=c(7,2,5,9),a=c(4,1,9,2),b=c(1,0,1,5),c=c(3,0,7,0),d=c(0,0,4,0))
df3 <- data.frame(Id=c(5,3,2,6),a=c(9,0,1,5),b=c(1,1,0,0),c=c(7,4,0,0),d=c(0,0,4,0))
##do.call(rbind, mget(ls(pattern = 'df[0-9]+'))) to make the initial data frame 
##If you are concerned about dataframe binding efficiency,
##rbindlist function in data.table package does the dataframes binding very efficiently.
##Just try combined <- rbindlist(list(df1, df2, df3))

library(data.table)

DTs = lapply(list(`DP-C_CYT_1_1`,`DP-C_CYT_1_2`,`DP-C_CYT_2_1`,`DP-C_CYT_2_2`
                  ,`DP-C_CYT_3_1`,`DP-C_CYT_3_2`), data.table)

Id_keep = rbindlist(lapply(DTs, `[`, j = "matched_peptide"))[, .N, by=matched_peptide][N >= 6L, matched_peptide];Id_keep
##Wert auf 3 wenn in allen 3 gefunden werden soll, Wert muss in mind. 3 Tabellen enthalten sein
Id_keep = rbindlist(lapply(DTs, `[`, j = c("Id","d")))[, .N, by=c("Id","d")][N >= 3L, Id];Id_keep 


DT_keep = Reduce(funion, DTs)[matched_peptide %in% Id_keep];DT_keep


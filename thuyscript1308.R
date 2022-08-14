
library(dplyr)
com <- readRDS("../Inp/Compound_properties.RDS")
gene <- readRDS("../Inp/GDSC_cancerdrivers_mut.RDS")
ic50 <- readRDS("../Inp/IC50dt_lungcancer_MTT_inctime2.RDS")

inctime <- ic50 %>% filter(inc_time == "68-72 hrs")

nic <- inctime %>% select (assay_cell_type, molregno, standard_value, inc_time)

matches <- inner_join(nic, gen, c("assay_cell_type" = "model_name")) # 5781 obs
# Model name and assay_cell_type are not in the sainme syntax, 
# so we curate to increase the number of observation
nic$assay_cell_type[which(nic$assay_cell_type == 'A-427')] = 'A427' # ChemBL says A-427,GDSC says A427
nic$assay_cell_type[which(nic$assay_cell_type == 'HCC827')] = 'HCC-827'# ChemBL says HCC827,GDSC says HCC-827
nic$assay_cell_type[which(nic$assay_cell_type == 'HCC78')] = 'HCC-78'
nic$assay_cell_type = toupper(nic$assay_cell_type) # for e.g CALU-3 and Calu-3
gen$model_name =  toupper(gen$model_name) # so we change them to all uppercase

matches <- inner_join(nic, gen, c("assay_cell_type" = "model_name")) # 5877 obs
match <- inner_join(nic, gen, c("assay_cell_type" = "model_name")) #5931 obs

tab1 <- table (match$molregno, match$assay_cell_type)
a1 <- as.data.frame(tab1)
testone <- a1 %>% filter(Freq == 1) #use this to test solution2
one1 <- a1 %>% filter(Freq == 1)
#-----------------------------------------#
library(tidyr)
unite_pair1<- unite(one1, col= "pair_1", c("Var1", "Var2"), sep=".")
unite_pair2 <- unite(match, col="pair_2", c("molregno", "assay_cell_type"), sep= ".")
pair.1 <- merge(unite_pair1, unite_pair2, by.x='pair_1', by.y='pair_2', all.x = TRUE)
#-----------------------------------------#

#now I want to split molregno and assay_cell_type in pair.1 but I met error using strsplit() 
#I can now only separate() because I know it better than strsplit()
#But separate() separates strings based on any non-alphanumeric value, 
#so if I have the name '234.NCI-H123' then it will only appear '234' and 'NCI'.
# to solve this, I found 2 solutions: 
#1: I can install package 'stringr' and use str_split_fixed
#2: I will change all assay_cell_type name so that the only non-alphanumeric appear in the name is '.', the use separate()
#I will only try solution2

# >table(unique(one1$Var2))
# there are 14 assay_cell_type names only
#[1] "A549"      "NCI-H446"  "NCI-H460"  "CALU-3"    "SK-LU-1"  
#[6] "SK-MES-1"  "NCI-H226"  "HCC-827"   "HOP-92"    "NCI-H358" 
#[11] "HOP-62"    "NCI-H1299" "HCC-78"    "NCIH446"

two$Var1[which(two$Var1 == 'NCI-H446')] = 'NCIH446'
two$Var1[which(two$Var1 == 'NCI-H460')] = 'NCIH460'
two$Var1[which(two$Var1 == 'NCI-H226')] = 'NCIH226'
two$Var1[which(two$Var1 == 'NCI-H358')] = 'NCIH358'
two$Var1[which(two$Var1 == 'NCI-H1299')] = 'NCIH1299'
two$Var1[which(two$Var1 == 'CALU-3')] = 'CALU3'
two$Var1[which(two$Var1 == 'SK-LU-1')] = 'SKLU1'
two$Var1[which(two$Var1 == 'SK-MES-1')] = 'SKMES1'
two$Var1[which(two$Var1 == 'HCC-827')] = 'HCC827'
two$Var1[which(two$Var1 == 'HOP-92')] = 'HOP92'
two$Var1[which(two$Var1 == 'HOP-62')] = 'HOP62'
two$Var1[which(two$Var1 == 'HCC-78')] = 'HCC78'

#using add row I accidentally add 2 more row in 'two' table so now i have to remove them from 'two'
# with the 2 extra rows, I notice they have N/A in Var2 column >> try to remove row with N/A
#there are 3 ways to remove a row with N/A
#1:two1 <- two %>% drop_na(Var2)
#2: two[!is.na(two$Var2,)]
#3: subset(two, !is.na(Var2))

#after rename 'two', do the same with 'match'

match$assay_cell_type[which(match$assay_cell_type == 'SK-MES-1')] = 'SKMES1'
match$assay_cell_type[which(match$assay_cell_type == 'NCI-H460')] = 'NCIH460'
match$assay_cell_type[which(match$assay_cell_type == 'NCI-H358')] = 'NCIH358'
match$assay_cell_type[which(match$assay_cell_type == 'HCC-827')] = 'HCC827'
match$assay_cell_type[which(match$assay_cell_type == 'NCI-H446')] = 'NCIH446'
match$assay_cell_type[which(match$assay_cell_type == 'NCI-H226')] = 'NCIH226'
match$assay_cell_type[which(match$assay_cell_type == 'CALU-3')] = 'CALU3'
match$assay_cell_type[which(match$assay_cell_type == 'HOP-62')] = 'HOP62'
match$assay_cell_type[which(match$assay_cell_type == 'HOP-92')] = 'HOP92'
match$assay_cell_type[which(match$assay_cell_type == 'NCI-H522')] = 'NCIH522'
match$assay_cell_type[which(match$assay_cell_type == 'NCI-H1299')] = 'NCIH1299'
match$assay_cell_type[which(match$assay_cell_type == 'CALU-6')] = 'CALU6'
match$assay_cell_type[which(match$assay_cell_type == 'HCC-78')] = 'HCC78'
match$assay_cell_type[which(match$assay_cell_type == 'SK-LU-1')] = 'SKLU1'

unite_pair1 <- unite_pair1 %>% separate(pair_1, c('assay_cell_type', 'molregno'))
unite_pair1 <- unite_pair1 %>% unite(pair_1, molregno, assay_cell_type, sep = '.')
unite_pair2 <- unite(match, col="pair_2", c("molregno", "assay_cell_type"), sep= ".")
pair.2 <- merge(unite_pair1, unite_pair2, by.x='pair_1', by.y='pair_2', all.x = TRUE)
pair.2 <- pair.2 %>% separate(pair_1, c('molregno', 'assay_cell_type'))
finaltable <- merge(pair.2, com, by.x="molregno", by.y="molregno", all.x=TRUE)
View(finaltable)

#save my table
setwd("/Users/siri/Desktop/intern/thuy_dissertation_2022.Rproj")
saveRDS(finaltable, file = 'finaltable.RData')
f <- readRDS('finaltable.RData')
View(f)

#####################################################

#try linear regression model
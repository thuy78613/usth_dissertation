
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

tab <- table (matches$molregno, matches$assay_cell_type)
a <- as.data.frame(tab)
one <- a %>% filter(Freq == 1)
#cách 1 là e paste Var1 và Var2 thành 1 cột mới (pairs), xong paste molregno và assay-cell_type thành 1 cột mới (pairs), xong e join 2 bàng với nhau theo cột pairs
#cách 2 là e làm một vòng for, for (i in 1:nrow(one)){ }
#sau đó đối với từng row, thì lọc row ở bảng matches có molregno == one$var1[i] & assay_cell_type == one$var2[i]
#rồi rbind() các dòng thành 1 bảng matches mới

#1: 
library(tidyr)
unite_pair_1<- unite(one, col= "pair_1", c("Var1", "Var2"), sep="-")
unite_pair_22 <- unite(matches, col="pair_2", c("molregno", "assay_cell_type"), sep= "-")

#1:
?strsplit()


unlist(strsplit(pair$pair_1[1],'-'))[1]
paste(unlist(strsplit(pair$pair_1[1],'-'))[-1],collapse = '-')

pair$molregno =unlist(lapply(strsplit(pair$pair_1,'-'),'[',1))
pair$cellname = unlist(lapply(strsplit(pair$pair_1,'-'),
              function(x){
                paste(unlist(x)[-1],collapse = '-')
              }))

table(pair$cellname)
length(unique(pair$cellname))
#2:
dotunite_pair_1 <- unite(one, col= "pair_1", c("Var1", "Var2"), sep=".")
dotunite_pair_2 <- unite(matches, col="pair_2", c("molregno", "assay_cell_type"), sep= ".")
fpair <- merge(dotunite_pair_1, dotunite_pair_2, by.x="pair_1", by.y="pair_2", all.x=TRUE)
head(strsplit(fpair$pair_1, split = '.')) #have unwanted output
counttable <- separate(fpair, col=pair_1, into= c('molregno', 'assay_cell_type'), sep='.') #still have error splitting

#I want to try paste() but got an error
colnames(one)
# >"Var1" "Var2" "Freq"
head(paste(Var1, Var2, sep = '.'))
# >Error in paste(Var1, Var2, sep = ".") : object 'Var1' not found




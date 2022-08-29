library(rpart)
library(rpart.plot)
library(caTools)
data <- read.csv("/Users/siri/Desktop/intern/thuy_dissertation_2022.Rproj/finaltable_script/FinalTableIncludedYear.csv")
ft <- subset(data, select=-c(Freq, inc_time, full_molformula, cx_most_bpka, cx_most_apka, assay_cell_type,molregno, year))
ft <- na.omit(ft)
set.seed(1)
sample <- sample.split(ft$standard_value, SplitRatio = 0.8)
train <- subset(ft, sample==TRUE)
train$Response=ifelse(train$standard_value >= median(train$standard_value),'NonR','Resp')
test <- subset(ft, sample==FALSE)
test$Response=ifelse(test$standard_value >= median(train$standard_value),'NonR','Resp')

fml = paste('Response~', paste(colnames(train)[-c(1,635)], collapse = '+'))
tableTree <- rpart(data = train,fml, method="class")

fml = paste('Response~',paste(colnames(train)[-c(1,632)], collapse = '+'))
tableTree <- rpart(data = train,fml, method="class")

tableTree$frame$var


#tree.pred <- predict(tableTree, test)
#head(tree.pred)
#NonR      Resp
#8  0.1475862 0.8524138
#9  0.2776911 0.7223089
#14 0.6849315 0.3150685
#61 0.6849315 0.3150685
#62 0.6849315 0.3150685
#63 0.6849315 0.3150685

#tree <- rpart(standard_value ~ ., method="class", data=ft)
tree.pred <- predict(tableTree, test,type = 'class')
head(tree.pred)
#8    9   14   61   62   63 
#Resp Resp NonR NonR NonR NonR 
#Levels: NonR Resp
test$Response
rownames(test)

mcc = function(pred,obs){
  TP = length(which(pred == "Resp" & obs == 'Resp'))    # cells predicted sensitive which are actually sensitive (cell sensitive is +ve)
  TN = length(which(pred == "NonR" & obs == 'NonR'))  # cells predicted resistant which are actually resistant (cell resistant is -ve)
  FP = length(which(pred == "Resp" & obs == 'NonR')) # cells predicted sensitive which are actually resistant (wrong prediction)
  FN = length(which(pred == "NonR" & obs == 'Resp')) # cells predicted resistant which are actually sensitive (wrong prediction)
  dum1 = ((TP*TN)-(FP*FN))
  dum2 = sqrt(1.0*(TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)) # 1.0 makes the multiplication float, thus avoiding integer overflow
  if(dum2 == 0) dum2 = 1
  mcc = dum1/dum2 
  return(mcc)
}
mcc(tree.pred,test$Response)
#[1] 0.5900372

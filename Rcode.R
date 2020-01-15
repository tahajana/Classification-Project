#rm(list=ls())


# All the packages used
library(MASS)
library(ggplot2)
library("ggthemes")
library("GGally")
library("extracat")
library(hdrcde)
library(KernSmooth)
library("ggplot2")
library("gridExtra")
library("vcd")
library(class)
library(rpart) 
library(rattle)
library(tree) 
library(randomForest)
library(gbm)
library(e1071)
library(mclust)
library(xtable)

# reading in the data
data1 <- read.csv("pulsar_stars.csv")
# renaming the data to smaller names:
names(data1) <- c("meanIP","sdIP","excessIP","skewnessIP","meanDM","sdDM","excessDM","skewnessDM","target")
summary(data1)

ggpairs(data1,aes(color = V9, alpha = 0.4))


# Scaling the data:
data2 <- scale(data1[,-9])
data2 <- cbind(data2,data1[,9])
data3 <- as.data.frame(data2)
data3$V9 = as.factor(data3$V9)
summary(data3)

# Making up our TRAIN and TEST sample using stratified sampling (set.seed(1))
p.tr <- 0.7
y.1 <- data3$V9 == 1
y.0 <- data3$V9 == 0
set.seed(1)
train.1 <- sample(row.names(data3[y.1,]), floor(nrow(data3[y.1,]) * p.tr) )
train.0 <- sample(row.names(data3[y.0,]), floor(nrow(data3[y.0,]) * p.tr) )
train <- as.numeric( c(train.0, train.1) )
d.train <- data3[train, ]
d.test <- data3[-train, ]



#LDA 
lda_data<-lda(x=data3[,-9],grouping=data3[,9],subset=train)
pred_class<-predict(lda_data,data3[-train,-9])$class
class.table <- table(data3[-train,9],pred_class)
xtable(class.table)
ari.LDA <- classAgreement(class.table)$crand
misc.LDA <- 1-classAgreement(class.table)$diag

# QDA 
set.seed(1)
qda_data<-qda(x=data3[,-9],grouping=data3[,9],subset=train)
pred_class<-predict(qda_data,data3[-train,-9])$class
class.table <- table(data3[-train,9],pred_class)
xtable(class.table)
ari.QDA <- classAgreement(class.table)$crand
misc.QDA <- 1-classAgreement(class.table)$diag


labels <-data3[train,9]
kn <- seq(1:30)
error = kn
 for(k in kn){
   cat(paste("There are",k,"nearest neighbours"),"\n")
   set.seed(1)
   output <- knn(train=d.train[,-9], test=d.test[,-9], cl=labels, k=k, prob=TRUE)
   class.table <- table(d.test[,9],output)
   error[k] <- classAgreement(class.table)$crand
}
plot(kn,error,typ="l",col=2,xlab= "k", xlim=c(0,30),ylab="ARI",ylim=c(0.8,0.86))
# We get our best ARI when k = 3
set.seed(1)
output <- knn(train=d.train[,-9], test=d.test[,-9], cl=labels, k=3, prob=TRUE)
class.table <- table(d.test[,9],output)
xtable(class.table)
ari.knn <- classAgreement(class.table)$crand
misc.knn <- 1-classAgreement(class.table)$diag


# tree
data3_tree <- rpart(V9~., data=data3, method="class")
# Plot decision tree for whole data
fancyRpartPlot(data3_tree, main="Classification Tree for the Pulsar Data")
# Classification application
set.seed(1)
data3_tree2 <- rpart(V9 ~., data=data3, subset = train,method="class",minsplit=25)
fancyRpartPlot(data3_tree2, main="Classification Tree for the train set")

# boxplot
boxplot(excessIP~V9, data=d.train, main="Boxplot of the training set", xlab = "group")
abline(v = 0,h =0.64,typ="l",col=2)

class.table <- table(data3[-train, "V9"],predict(data3_tree2, data3[-train,-9], type = "class"))
xtable(class.table)
ari.dtree <- classAgreement(class.table)$crand
misc.dtree <- 1-classAgreement(class.table)$diag


# Bagging
# I commented parameter tuning code because they take a lot of time to run.

# lets check for m=1,2,3,4
# M<-(1:20)*50
# error<-M
#  for(m in M){
#    set.seed(1)
#    bag.data3=randomForest(V9~.,data=data3,subset=train,ntree = m,mtry=8,type="class")
#    data3.pred=predict(bag.data3,data3[-train,],type="class")
#    tab<-table(data3.pred,d.test[,9])
#    error[m/50]<-classAgreement(tab)$crand
# }
#  bag<-error
# # # in bagging, we got our highest ARI for number of trees = 200: 0.8611219
# for(m in M){
#    set.seed(1)
#    bag.data3=randomForest(V9~.,data=data3,subset=train,ntree = m,mtry=3,importance=TRUE,type="class")
#    data3.pred=predict(bag.data3,data3[-train,],type="class")
#    tab<-table(data3.pred,d.test[,9])
#    error[m/50]<-classAgreement(tab)$crand
#  }
#  rf3<-error
#  for(m in M){
#    set.seed(1)
#    bag.data3=randomForest(V9~.,data=data3,subset=train,ntree = m,mtry=2,importance=TRUE,type="class")
#    data3.pred=predict(bag.data3,data3[-train,],type="class")
#    tab<-table(data3.pred,d.test[,9])
#    error[m/50]<-classAgreement(tab)$crand
#  }
#  rf2<-error
#  for(m in M){
#    set.seed(1)
#    bag.data3=randomForest(V9~.,data=data3,subset=train,ntree = m,mtry=1,importance=TRUE,type="class")
#    data3.pred=predict(bag.data3,data3[-train,],type="class")
#    tab<-table(data3.pred,d.test[,9])
#    error[m/50]<-classAgreement(tab)$crand
#  }
#  rf1<-error
#  for(m in M){
#    set.seed(1)
#    bag.data3=randomForest(V9~.,data=data3,subset=train,ntree = m,mtry=4,importance=TRUE,type="class")
#    data3.pred=predict(bag.data3,data3[-train,],type="class")
#    tab<-table(data3.pred,d.test[,9])
#    error[m/50]<-classAgreement(tab)$crand
#  }
#  rf4<-error
# 
# # in random forest, we got our highest ARI for m = 4 and number of trees = 500: 0.8625785
# plot(M,bag,typ="l",col=2,xlab= "number of trees", xlim=c(0,1000),ylab="ARI",ylim=c(0.8,0.88))
# lines(M,rf3,typ="l",col=1)
# lines(M,rf2,typ="l",col=4)
# lines(M,rf1,typ="l",col=5)
# lines(M,rf4,typ="l",col=6)
# legend(x=1,y=0.844,c("m=19 (bagging)","m=3","m=2","m=1","m=4"),col=c(2,1,4,5,6),pch=c(2,3,4,5,6))

# Bagging using ntree = 200 ( chosen from the plot)
set.seed(1)
bag.data3=randomForest(V9~.,data=data3,subset=train,ntree = 200,mtry=8,importance=TRUE,type="class")
data3.pred=predict(bag.data3,data3[-train,],type="class")
tab<-table(d.test[,9],data3.pred)
xtable(tab)
ari.bagging<-classAgreement(tab)$crand
misc.bagging <- 1 - classAgreement(tab)$diag

# rf using ntree = 500 and m = 4 ( chosen from the plot)
set.seed(1)
rf.data3=randomForest(V9~.,data=data3,subset=train,ntree = 500,mtry=4,importance=TRUE,type="class")
data3.pred=predict(rf.data3,data3[-train,],type="class")
tab<-table(d.test[,9],data3.pred)
xtable(tab)
ari.rf<-classAgreement(tab)$crand
misc.rf <- 1 - classAgreement(tab)$diag

# variable importance
varImpPlot(rf.data3, main = "Variable importance Plot for Rf with m = 4 and ntree = 500")

# Boosting:
# I commented out parameter tuning for boostig, because they require a lot of time to run 

# lamda1 <- c(0.01,0.001)
# d <- c(1,2,3,4,5,6)
# n1 = (1:5)* 1000 
# a = data.frame(numberofTrees = NULL,d. = NULL, Lamda = NULL, testError = NULL)
#  for(m in c(4000,5000)) {
#    for (j in d) {
#      for (k in lamda1) {
#        set.seed(1)
#        boost.data3=gbm(V9~.,data= d.train,distribution="multinomial",n.trees=m,interaction.depth=j,shrinkage=k)
#        yhat.data3=predict(boost.data3,newdata=d.test,n.trees=m,distribution="multinomial",type="response")
#        class.pred<-rep(0,5370)
#        for(i in 1:5370){
#          which(yhat.data3[i,,1]==max(yhat.data3[i,,1]))->class.pred[i]
#        }
#        tab<-table(class.pred,d.test[,9])
#        ari = classAgreement(tab)$crand
#        a = rbind(a,data.frame(numberofTrees = m,d. = j, Lamda = k, testError = ari))
#        print(h)
#      }
#    }
# }

# write.csv(a,"boostParam3.csv")
# a <- read.csv(file = "boostingParameters.csv")


# The highest ARI was for the parameters: d = 6 and lamda = 0.01,and number of trees = 1000 
set.seed(1)
boost.data3=gbm(V9~.,data= d.train,distribution="multinomial",n.trees=1000,interaction.depth=6,shrinkage=0.01)
# Now I can use gbm.perf to get the optimal number of trees for prediction.
#ntree = gbm.perf(boost.data3)
yhat.data3=predict(boost.data3,newdata=d.test,n.trees=295,distribution="multinomial",type="response")
class.pred<-rep(0,5370)
for(i in 1:5370){
  which(yhat.data3[i,,1]==max(yhat.data3[i,,1]))->class.pred[i]
}
tab<-table(d.test[,9],class.pred)
xtable(tab)
ari.boost = classAgreement(tab)$crand
misc.boost <- 1 - classAgreement(tab)$diag


# The plots

# a <- read.csv(file = "boostingParameters3.csv")
# a = a[,-1]
# plot(a$numberofTrees[(a$d. == 1) & (a$Lamda == 0.010)],a$testError[(a$d. == 1) & (a$Lamda == 0.010)],typ="l",col=2,xlab= "number of trees",
#      xlim=c(1000,6500),ylab="ARI",ylim=c(0.79,0.87), main="Lambda = 0.01")
# lines(a$numberofTrees[(a$d. == 2) & (a$Lamda == 0.010)],a$testError[(a$d. == 2) & (a$Lamda == 0.010)],typ="l",col=1)
# lines(a$numberofTrees[(a$d. == 3) & (a$Lamda == 0.010)],a$testError[(a$d. == 3) & (a$Lamda == 0.010)],typ="l",col=4)
# lines(a$numberofTrees[(a$d. == 4) & (a$Lamda == 0.010)],a$testError[(a$d. == 4) & (a$Lamda == 0.010)],typ="l",col=5)
# lines(a$numberofTrees[(a$d. == 5) & (a$Lamda == 0.010)],a$testError[(a$d. == 5) & (a$Lamda == 0.010)],typ="l",col=6)
# lines(a$numberofTrees[(a$d. == 6) & (a$Lamda == 0.010)],a$testError[(a$d. == 6) & (a$Lamda == 0.010)],typ="l",col=3)
# legend("bottomright",c("d = 1","d = 2","d = 3","d = 4","d = 5", "d = 6"),col=c(2,1,4,5,6,3),pch=c(2,3,4,5,6,7))
# 
# plot(a$numberofTrees[(a$d. == 1) & (a$Lamda == 0.001)],a$testError[(a$d. == 1) & (a$Lamda == 0.001)],typ="l",col=2,xlab= "number of trees",
#      xlim=c(1000,6500),ylab="ARI",ylim=c(0.79,0.87), main="Lambda = 0.001")
# lines(a$numberofTrees[(a$d. == 2) & (a$Lamda == 0.001)],a$testError[(a$d. == 2) & (a$Lamda == 0.001)],typ="l",col=1)
# lines(a$numberofTrees[(a$d. == 3) & (a$Lamda == 0.001)],a$testError[(a$d. == 3) & (a$Lamda == 0.001)],typ="l",col=4)
# lines(a$numberofTrees[(a$d. == 4) & (a$Lamda == 0.001)],a$testError[(a$d. == 4) & (a$Lamda == 0.001)],typ="l",col=5)
# lines(a$numberofTrees[(a$d. == 5) & (a$Lamda == 0.001)],a$testError[(a$d. == 5) & (a$Lamda == 0.001)],typ="l",col=6)
# lines(a$numberofTrees[(a$d. == 6) & (a$Lamda == 0.001)],a$testError[(a$d. == 6) & (a$Lamda == 0.001)],typ="l",col=3)
# legend("bottomright",c("d = 1","d = 2","d = 3","d = 4","d = 5", "d = 6"),col=c(2,1,4,5,6,3),pch=c(2,3,4,5,6,7))



# MixtureDA
set.seed(1)
data3MclustDA <- MclustDA(data3[train,-9], data3[train,9])
summary(data3MclustDA, parameters = TRUE)
summ = summary(data3MclustDA, newdata = data3[-train,-9], newclass = data3[-train,9])
print(summ)
class.tab <- summ$tab.newdata
xtable(class.tab)
mixtureDA.ari <- classAgreement(class.tab)$crand
mixtureDA.misc <-  1 - classAgreement(class.tab)$diag

results <- data.frame(set_seed= 1,LDA = ari.LDA, ari.QDA, knn = ari.knn, tree = ari.dtree, bagging = ari.bagging, rf = ari.rf, boosting = ari.boost, mixtureDA = mixtureDA.ari)
results2 <- data.frame(set_seed= 1,LDA = misc.LDA, misc.QDA, knn = misc.knn, tree = misc.dtree, bagging = misc.bagging, rf = misc.rf, boosting = misc.boost, mixtureDA = mixtureDA.misc)


# Getting the errors for 10 random traing sets
for (j in 2:10) {
  set.seed(j)
  p.tr <- 0.7
  y.1 <- data3$V9 == 1
  y.0 <- data3$V9 == 0
  train.1 <- sample(row.names(data3[y.1,]), floor(nrow(data3[y.1,]) * p.tr) )
  train.0 <- sample(row.names(data3[y.0,]), floor(nrow(data3[y.0,]) * p.tr) )
  train <- as.numeric( c(train.0, train.1) )
  d.train <- data3[train, ]
  d.test <- data3[-train, ]
  #LDA 
  set.seed(j)
  lda_data<-lda(x=data3[,-9],grouping=data3[,9],subset=train)
  pred_class<-predict(lda_data,data3[-train,-9])$class
  class.table <- table(data3[-train,9],pred_class)
  ari.LDA <- classAgreement(class.table)$crand
  misc.LDA <- 1-classAgreement(class.table)$diag
  # QDA 
  set.seed(j)
  qda_data<-qda(x=data3[,-9],grouping=data3[,9],subset=train)
  pred_class<-predict(qda_data,data3[-train,-9])$class
  class.table <- table(data3[-train,9],pred_class)
  ari.QDA <- classAgreement(class.table)$crand
  misc.QDA <- 1-classAgreement(class.table)$diag
  #knn
  set.seed(j)
  labels <-data3[train,9]
  output <- knn(train=d.train[,-9], test=d.test[,-9], cl=labels, k=3, prob=TRUE)
  class.table <- table(d.test[,9],output)
  ari.knn <- classAgreement(class.table)$crand
  misc.knn <- 1-classAgreement(class.table)$diag
  #trees
  set.seed(j)
  data3_tree2 <- rpart(V9 ~., data=data3, subset = train,method="class",minsplit=25)
  class.table <- table(data3[-train, "V9"],predict(data3_tree2, data3[-train,-9], type = "class"))
  ari.dtree <- classAgreement(class.table)$crand
  misc.dtree <- 1-classAgreement(class.table)$diag
  #rf
  set.seed(j)
  rf.data3=randomForest(V9~.,data=data3,subset=train,ntree = 500,mtry=4,importance=TRUE,type="class")
  data3.pred=predict(rf.data3,data3[-train,],type="class")
  tab<-table(d.test[,9],data3.pred)
  ari.rf<-classAgreement(tab)$crand
  misc.rf <- 1-classAgreement(tab)$diag
  #bagging
  set.seed(j)
  bag.data3=randomForest(V9~.,data=data3,subset=train,ntree = 200,mtry=8,importance=TRUE,type="class")
  data3.pred=predict(bag.data3,data3[-train,],type="class")
  tab<-table(d.test[,9],data3.pred)
  ari.bagging<-classAgreement(tab)$crand
  misc.bagging <- 1-classAgreement(tab)$diag
  #boosting
  set.seed(j)
  boost.data3=gbm(V9~.,data= d.train,distribution="multinomial",n.trees=1000,interaction.depth=6,shrinkage=0.01)
  yhat.data3=predict(boost.data3,newdata=d.test,n.trees=295,distribution="multinomial",type="response")
  class.pred<-rep(0,5370)
  for(i in 1:5370){
    which(yhat.data3[i,,1]==max(yhat.data3[i,,1]))->class.pred[i]
  }
  tab<-table(d.test[,9],class.pred)
  ari.boost = classAgreement(tab)$crand
  misc.boost <- 1-classAgreement(tab)$diag
  # MixtureDA
  set.seed(j)
  data3MclustDA <- MclustDA(data3[train,-9], data3[train,9])
  summ = summary(data3MclustDA, newdata = data3[-train,-9], newclass = data3[-train,9])
  class.tab <- summ$tab.newdata
  mixtureDA.ari <- classAgreement(class.tab)$crand
  mixtureDA.misc <- 1-classAgreement(class.tab)$diag
  results <- rbind(results,data.frame(set_seed= j,LDA = ari.LDA, ari.QDA, knn = ari.knn, tree = ari.dtree, bagging = ari.bagging, rf = ari.rf, boosting = ari.boost, mixtureDA = mixtureDA.ari))
  results2 <- rbind(results2,data.frame(set_seed= j,LDA = misc.LDA, misc.QDA, knn = misc.knn, tree = misc.dtree, bagging = misc.bagging, rf = misc.rf, boosting = misc.boost, mixtureDA = mixtureDA.misc))
}

results <- read.csv("results.csv")
data5 = data5[,-1]
results2 <- read.csv("results2.csv")
data6 = data6[,-1]
xtable(data5, digits = 4, caption = "ARI for different training sets")
xtable(data6, digits = 4,caption = "Missclassification error for different training sets")

# Calculating the average ARI and missclassification rate of the 10 random traing sets
LDAARI = sum(results$LDA)/10
QDAARI = sum(results$ari.QDA)/10
knnARI = sum(results$knn)/10
treeARI = sum(results$tree)/10
bagARI = sum(results$bagging)/10
rfARI = sum(results$rf)/10
boostARI = sum(results$boosting)/10
mixtARI = sum(results$mixtureDA)/10

LDAmisc = sum(results2$LDA)/10
QDAmisc = sum(results2$misc.QDA)/10
knnmisc = sum(results2$knn)/10
treemisc = sum(results2$tree)/10
bagmisc = sum(results2$bagging)/10
rfmisc = sum(results2$rf)/10
boostmisc = sum(results2$boosting)/10
mixtmisc = sum(results2$mixtureDA)/10

# The average ARI and missclassification rate for the 10 different training sets.
resultsARI <- data.frame(LDA = LDAARI,QDA = QDAARI, knn = knnARI, tree = treeARI, bagging = bagARI, rf = rfARI, boosting = boostARI, mixtureDA = mixtARI)
resultsMisc <- data.frame(LDA = LDAmisc, QDA = QDAmisc, knn = knnmisc, tree = treemisc, bagging = bagmisc, rf = rfmisc, boosting = boostmisc, mixtureDA = mixtmisc)




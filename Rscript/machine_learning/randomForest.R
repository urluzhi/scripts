library('e1071')
require(randomForest)
input=read.csv("bins.training-5classes.sampled.csv")
output='5classes.rf.out'
roc_file="5classes.rf.roc.pdf"
importance_file="5classes.rf.importance.pdf"
model_file="5classes.rf.model"
dataall=input[,c(8:12,15:16,18:19)]
classesall=subset(input,select=X1.Annotation)

#Generate training and testing sets
nall=nrow(input)
ntrain=2*floor(nall/3)
datatrain <- dataall[1:ntrain,]
classestrain <- classesall[1:ntrain,]
ntrain=ntrain+1
datatest <- dataall[ntrain:nall,]
classestest <- classesall[ntrain:nall,]



#10 fold cross-validation on training set
obj=tune(randomForest,train.x=datatrain,train.y=classestrain,ntree=1000)
write(paste("Error estimation of randomForest using 10-fold cross validation:",obj$performance$error,spe="\t"),file=output)



#test on testing set
#predict on testing test
rf=randomForest(datatrain,classestrain,xtest=datatest,ytest=classestest,ntree=1000,importance = TRUE,proximity = TRUE)
save(rf, file=model_file)
prob=rf$test$votes
#AUC of ROC
library('ROCR')
write("Performance on testing set",file=output,append=T)
for(i in c('UTR','exon_CDS','ncRNA_sampled','unexpressed_intergenic','pgene_exon')) {
	test=as.vector(classestest)
	test[test==i]=1
	test[test!=1]=0
	predout=prediction(prob[,i],test)
	auc=performance(predout, "auc")
	write(paste(i," AUC:\t",auc@y.values[[1]],spe=''),file=output,append=T)
}
#Draw ROC plot
pref=performance(predout,"tpr","fpr")
pdf(roc_file)
plot(pref, col='red',main=i )
dev.off()
#Draw importance plot
write("Importance of each feature",file=output,append=T)
write.table(rf$importance,sep="\t", file=output,append=T,row.name=T,col.name=T)
train=as.vector(classestrain)
train[train=='pgene_exon']=1
train[train!=1]=0
train=as.factor(train)
rf=randomForest(datatrain,train,ntree=1000,importance = TRUE, proximity = TRUE)
pdf(importance_file,width=15,height=5)
varImpPlot(rf)
dev.off()

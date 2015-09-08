library('e1071')
input=read.csv("bins.training-5classes.sampled.csv")
output='5classes.svm.out'
model_file="5classes.svm.model"
pdf("5classes.svm.roc.pdf")
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
#model_train<-svm(datatrain,classestrain,type='C',kernel='radial',cross=10,probability=T)
#write(paste(" 10 fold cross-validation Accuracy:",model_train$tot.accuracy,sep='\t'),file=output)
#obj=tune(svm, train.x=datatrain, train.y=classestrain, validation.x=datatest, validation.y=classestest, ranges = list(gamma = 2^(-1:1), cost = 2^(2:4)), control = tune.control(sampling = "fix"))
obj=tune(svm,train.x=datatrain,train.y=classestrain)
write(paste("Error estimation of SVM using 10-fold cross validation:",obj$performance$error,spe="\t"),file=output)


#test on testing set
#predict on testing test
model<-svm(datatrain, classestrain,type='C',kernel='radial',probability=T)
save(model, file=model_file)
pred<-predict(model,datatest,probability=T)
prob=attr(pred, "probabilities")
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
plot(pref, col='red',main=i )
dev.off()

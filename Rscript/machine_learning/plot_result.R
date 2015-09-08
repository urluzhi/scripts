library('e1071')
require(randomForest)
require(RColorBrewer)
input=read.csv("bins.training-5classes.sampled.csv")
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

#load model
load(model_file)
pdf("mdsplot2.pdf")
MDSplot(rf,classestrain,k=2)
dev.off()
	
pdf("mdsplot3.pdf")
MDSplot(rf,classestrain,k=3)
dev.off()
	
pdf("mdsplot9.pdf")
MDSplot(rf,classestrain,k=9)
dev.off()

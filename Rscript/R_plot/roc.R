library(ROCR)
	data=read.table("roc.example.in")
	pdf("roc.pdf")
	type=data[,1]
	pred <- prediction(data[,3], type)
	perf <- performance(pred, measure = "tpr", x.measure = "fpr") 
	plot(perf, col=2)
	dev.off()


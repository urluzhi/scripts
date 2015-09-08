##load the needed package
library("class")

data=read.table("roc.example.in");

##plot the ppv curve
pdf("ppv.pdf")
tp <- 0
fp <- 0
tn <- 0
fn <- 0
ppv <- c()
sen <- c()
positives <- c()


cutoff <- seq(from=0,to=1,by=0.01)
x<- 1
for (i in cutoff) {
	tp <- 0
	fp <- 0
	tn <- 0
	fn <- 0
	for (index in (1:length(data[,1]))) {
		if(data[index,3] >= i && data[index,1] == 1) tp <- tp + 1
		if(data[index,3] >= i && data[index,1] == 0) fp <- fp + 1
		if(data[index,3] < i && data[index,1] == 0) tn <- tn + 1
		if(data[index,3] < i && data[index,1] == 1) fn <- fn + 1
	}
	ppv[x] <- tp/(tp+fp)
	sen[x] <- tp/(tp+fn)
	positives[x]=(tp+fp)/length(data[,1])
	x <- x + 1
}
plot(sen,ppv,type="l",main="PPV",col=2,xlab="X",ylab="PPV")
lines(positives,ppv,type="l",col=3)
legend("topright",legend=c("X=Sensitivity","X=Rate of Positives"),text.col=c(2,3))
dev.off()

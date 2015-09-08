#data=read.csv("bins.chr1.known_types.sample2000.csv")
data=read.csv("bins.chr1.known_types.csv")

#normalize the data
data$X8.Identity=data$X8.Identity/100
data$X10.SCI=data$X10.SCI/data$X8.Identity
for (i in 10) { #normalize SCI from 0 to 1
		maxValue=max(data[,i])
				minValue=min(data[,i])
					range=maxValue-minValue
						data[,i]=(data[,i]-minValue)/range
}
for (i in 12:17){ # log expression values
		data[,i]=log(data[,i]+1)
}

#Merge sub-types of different elements
data[,2]=sub ("ncRNA_3019","ncRNA_New", data[,2])
data[,2]=sub ("ncRNA_selected","ncRNA", data[,2])
data[,2]=sub ("exon_CCDS","CDS", data[,2])
data[,2]=sub ("five_prime_UTR","UTR", data[,2])
data[,2]=sub ("three_prime_UTR","UTR", data[,2])
data[,2]=sub ("ancestral_repeat","AP", data[,2])

#plot
png("scatterplot.png",width=1024,height=1024)
#par(mfrow=c(5,5))
for (j in 17) { 
	for (k in 9) {
		for (i in c("AP","CDS","UTR","ncRNA_New","ncRNA") ){
				x=data[data[,2]==i,j]
						y=data[data[,2]==i,k]
							if (i=="AP") {
										plot(x,y, col=rgb(0,255,0,100,maxColorValue=255),pch=16,ylim=c(-7,5),xlab=colnames(data)[j],ylab=colnames(data)[k])
												}
					else if (i=="CDS") {
								points(x,y, col=rgb(255,0,0,100,maxColorValue=255), pch=16)
										}
					else if (i=="UTR") {
								points(x,y, col=rgb(255,255,0,100,maxColorValue=255), pch=16)
										}
					else if(i=="ncRNA_New") { 		
								points(x,y, col=rgb(173,216,230,150,maxColorValue=255), pch=16)
										}
					else if(i=="ncRNA") { 		
								points(x,y, col=rgb(0,0,255,180,maxColorValue=255), pch=16)
										}

		}
	}
}
dev.off()


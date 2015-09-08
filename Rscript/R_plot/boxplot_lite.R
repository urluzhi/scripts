#load library
library("lattice")
#read the data table
data=read.table("boxplot_lite_example.tsv",header=T)
#name the output file format
png("boxplot_lite.png",width=800,height=600)
#plot
attach(data)
boxplot(GC~Annotation,xlab="Gene class",ylab="GC%",ylim=c(0,1),col=c("red","yellow","green","blue","light pink","light blue","grey"))
#Write the file
dev.off()

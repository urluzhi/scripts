
# need gplots to use heatmap.2
require("gplots")
#for ( i in c("PHA4_PQM")) {
for ( i in c("HOX")) {
for ( j in c(".novelcro")) {
mysig=paste(i,j,".bins.signal",sep="")
myinput=paste(i,j,".bins.input",sep="")
mypdf=paste(i,j,".pdf",sep="")

# open an image file
pdf(mypdf,width=10,height=10)
sig = read.table(mysig,header=TRUE )
input = read.table(myinput,header=TRUE )
data=sig/input
noise1=apply(sig,1,min)
noise2=apply(input,1,min)
mdata=data[noise1>20&noise2>20,]
m=cor(log(mdata))


mn <- -0.8
mx <- 0.8
n <- 98		# the color steps, 98 color levels.
density.info <- "none"	# whether to display density in the key
bias <- 1	# parameter for colorRampPalette. 1 is what I used
mc <- matrix(as.character(round(m, 2)), ncol=dim(m)[2])     # the note on the matrix
breaks <- seq(mn, mx, (mx-mn)/(n))
cr <- colorRampPalette(colors = c("#2927FF","#FFFFFF","#DF5C5C"), bias=bias)
#heatmap.2(m, col = cr(n), breaks=breaks, trace="none", keysize=0.80000, density.info=density.info, margins=c(12.000000,12.000000), cexRow=1.200000, cexCol=1.200000)
#heatmap.2(m, col = cr(n),cellnote=mc,notecol="black",notecex=1.8,dendrogram = "none",Rowv = NA, Colv = NA, breaks=breaks, trace="none", keysize=0.80000, density.info=density.info, margins=c(12.000000,12.000000), cexRow=1.200000, cexCol=1.200000)
heatmap.2(m, col = cr(n),breaks=breaks, trace="none", cellnote=mc,notecol="black",keysize=0.80000, notecex=1.8,density.info=density.info, margins=c(12.000000,12.000000), cexRow=1.200000, cexCol=1.200000)
dev.off()
}
}

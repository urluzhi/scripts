
# need gplots to use heatmap.2
require("gplots")

# open an image file
png("sample_heatmap.R.png",width=2048,height=2048)
m <- matrix(runif(n=30*30, min=-1.0, max=1.0), ncol=30)
#m <- read.csv("file.csv")  #as.matrix
mn <- -1
mx <- 1
n <- 98		# the color steps, 98 color levels.
density.info <- "none"	# whether to display density in the key
bias <- 1	# parameter for colorRampPalette. 1 is what I used
mc <- matrix(as.character(round(m, 2)), ncol=dim(m)[2])     # the note on the matrix
breaks <- seq(mn, mx, (mx-mn)/(n))
cr <- colorRampPalette(colors = c("#2927FF","#FFFFFF","#DF5C5C"), bias=bias)
#heatmap.2(m, col = cr(n), breaks=breaks, trace="none",  notecol="black", notecex=1.800000, keysize=0.500000, density.info=density.info, margins=c(27.000000,27.000000), cexRow=2.200000, cexCol=2.200000)
heatmap.2(m, col = cr(n), breaks=breaks, trace="none", cellnote=mc, notecol="black", notecex=1.800000, keysize=0.500000, density.info=density.info, margins=c(27.000000,27.000000), cexRow=2.200000, cexCol=2.200000)
dev.off()

data=read.table("output/ALL-UNION.minrun50_maxgap5_threshold1.csv", header = TRUE, sep = ",")
pdf("srna_TAR_expression.hist.pdf")
#small RNAseq max
plot (density(log2(data[data$X27.Covered_Annotation...50..=="exon_CDS","X20.sRNAseq_max_all"]) ), col=2, main="sRNAseq (Max)", xlab="log2(DCPM)")
lines (density(log2(data[data$X27.Covered_Annotation...50..=="ncRNA_sampled","X20.sRNAseq_max_all"]) ), col=3)
lines (density(log2(data[data$X27.Covered_Annotation...50..=="intronic","X20.sRNAseq_max_all"]) ), col=4)
#lines (density(log2(data[data$X27.Covered_Annotation...50..=="intergenic","X20.sRNAseq_max_all"]) ), col=5)
lines (density(log2(data[data$X27.Covered_Annotation...50..=="three_prime_UTR","X20.sRNAseq_max_all"]) ), col=5)
lines (density(log2(data[data$X27.Covered_Annotation...50..=="five_prime_UTR","X20.sRNAseq_max_all"]) ), col=6)
legend("topright",legend=c("Confirmed CDS","Sampled ncRNA","Intronic Region","Three Prime UTR","Five Prime UTR"),text.col=c(2,3,4,5,6))


#poly-A RNAseq
#data$X19.RNAseq_max_all	
plot (density(log2(data[data$X27.Covered_Annotation...50..=="exon_CDS","X19.RNAseq_max_all"]) ), col=2, main="Poly-A RNAseq (Max)", xlab="log2(DCPM)")
lines (density(log2(data[data$X27.Covered_Annotation...50..=="ncRNA_sampled","X19.RNAseq_max_all"]) ), col=3)
lines (density(log2(data[data$X27.Covered_Annotation...50..=="intronic","X19.RNAseq_max_all"]) ), col=4)
#lines (density(log2(data[data$X27.Covered_Annotation...50..=="intergenic","X19.RNAseq_max_all"]) ), col=5)
lines (density(log2(data[data$X27.Covered_Annotation...50..=="three_prime_UTR","X19.RNAseq_max_all"]) ), col=5)
lines (density(log2(data[data$X27.Covered_Annotation...50..=="five_prime_UTR","X19.RNAseq_max_all"]) ), col=6)
legend("topright",legend=c("Confirmed CDS","Sampled ncRNA","Intronic Region","Three Prime UTR","Five Prime UTR"),text.col=c(2,3,4,5,6))

#Tiling array Poly-A
#data$X22.Array_L2_polyA
plot (density(data[data$X27.Covered_Annotation...50..=="exon_CDS","X22.Array_L2_polyA"]) , col=2, main="Poly-A RNA Tiling Array (Max)", xlab="log2(Intensity)")
lines (density(data[data$X27.Covered_Annotation...50..=="ncRNA_sampled","X22.Array_L2_polyA"]) , col=3)
lines (density(data[data$X27.Covered_Annotation...50..=="intronic","X22.Array_L2_polyA"]) , col=4)
#lines (density(data[data$X27.Covered_Annotation...50..=="intergenic","X22.Array_L2_polyA"]), col=5)
lines (density(data[data$X27.Covered_Annotation...50..=="three_prime_UTR","X22.Array_L2_polyA"]), col=5)
lines (density(data[data$X27.Covered_Annotation...50..=="five_prime_UTR","X22.Array_L2_polyA"]), col=6)
legend("topright",legend=c("Confirmed CDS","Sampled ncRNA","Intronic Region","Three Prime UTR","Five Prime UTR"),text.col=c(2,3,4,5,6))

#Tiling array total
#data$X21.Array_max_no_polyA
plot (density(data[data$X27.Covered_Annotation...50..=="exon_CDS","X21.Array_max_no_polyA"]) , col=2, main="Total RNA Tiling Array (Max)", xlab="log2(Intensity)")
lines (density(data[data$X27.Covered_Annotation...50..=="ncRNA_sampled","X21.Array_max_no_polyA"]) , col=3)
lines (density(data[data$X27.Covered_Annotation...50..=="intronic","X21.Array_max_no_polyA"]) , col=4)
#lines (density(data[data$X27.Covered_Annotation...50..=="intergenic","X21.Array_max_no_polyA"]), col=5)
lines (density(data[data$X27.Covered_Annotation...50..=="three_prime_UTR","X21.Array_max_no_polyA"]), col=5)
lines (density(data[data$X27.Covered_Annotation...50..=="five_prime_UTR","X21.Array_max_no_polyA"]), col=6)
legend("topright",legend=c("Confirmed CDS","Sampled ncRNA","Intronic Region","Three Prime UTR","Five Prime UTR"),text.col=c(2,3,4,5,6))
dev.off()



data=read.table("features/intergenic_gold_bins.expression.max", sep = ",")

cat("RNAseq_max_all\n")
i=2
mean(log2(data[data[,i]>0,i]))
median(log2(data[data[,i]>0,i]))
sd(log2(data[data[,i]>0,i]))

cat("sRNAseq_max_all")
i=3
mean(log2(data[data[,i]>0,i]))
median(log2(data[data[,i]>0,i]))
sd(log2(data[data[,i]>0,i]))

cat("Array_max_polyA")
i=4
mean(data[data[,i]>1,i])
median(data[data[,i]>1,i])
sd(data[data[,i]>1,i])

cat("Array_max_no_polyA")
i=5
mean(data[data[,i]>1,i])
median(data[data[,i]>1,i])
sd(data[data[,i]>1,i])

data=read.table("features/exon_CDS_bins.expression.max", sep = ",")

cat("RNAseq_max_all\n")
i=2
mean(log2(data[data[,i]>0,i]))
median(log2(data[data[,i]>0,i]))
sd(log2(data[data[,i]>0,i]))

cat("sRNAseq_max_all")
i=3
mean(log2(data[data[,i]>0,i]))
median(log2(data[data[,i]>0,i]))
sd(log2(data[data[,i]>0,i]))

cat("Array_max_polyA")
i=4
mean(data[data[,i]>1,i])
median(data[data[,i]>1,i])
sd(data[data[,i]>1,i])

cat("Array_max_no_polyA")
i=5
mean(data[data[,i]>1,i])
median(data[data[,i]>1,i])
sd(data[data[,i]>1,i])

data=read.table("features/exon_CDS_and_ncRNA_sample_bins.expression.max", sep = ",")

cat("RNAseq_max_all\n")
i=2
mean(log2(data[data[,i]>0,i]))
median(log2(data[data[,i]>0,i]))
sd(log2(data[data[,i]>0,i]))

cat("sRNAseq_max_all")
i=3
mean(log2(data[data[,i]>0,i]))
median(log2(data[data[,i]>0,i]))
sd(log2(data[data[,i]>0,i]))

cat("Array_max_polyA")
i=4
mean(data[data[,i]>1,i])
median(data[data[,i]>1,i])
sd(data[data[,i]>1,i])

cat("Array_max_no_polyA")
i=5
mean(data[data[,i]>1,i])
median(data[data[,i]>1,i])
sd(data[data[,i]>1,i])

data=read.table("features/intergenic_gold_bins.expression.max", sep = ",")

cat("RNAseq_max_all\n")
i=2
mean(log2(data[data[,i]>0,i]))
median(log2(data[data[,i]>0,i]))
sd(log2(data[data[,i]>0,i]))

cat("sRNAseq_max_all")
i=3
mean(log2(data[data[,i]>0,i]))
median(log2(data[data[,i]>0,i]))
sd(log2(data[data[,i]>0,i]))

cat("Array_max_polyA")
i=4
mean(data[data[,i]>1,i])
median(data[data[,i]>1,i])
sd(data[data[,i]>1,i])

cat("Array_max_no_polyA")
i=5
mean(data[data[,i]>1,i])
median(data[data[,i]>1,i])
sd(data[data[,i]>1,i])

cat("RNAseq_max_all\n")
i=2
mean(log2(data[,i]+1))
median(log2(data[,i]+1))
sd(log2(data[,i]+1))

cat("sRNAseq_max_all")
i=3
mean(log2(data[,i]+1))
median(log2(data[,i]+1))
sd(log2(data[,i]+1))


cat("Array_max_polyA")
i=4
mean(data[,i])
median(data[,i])
sd(data[,i])

cat("Array_max_no_polyA")
i=5
mean(data[,i])
median(data[,i])
sd(data[,i])

mean(log2(data[data[,i]>0,i]+1))
median(log2(data[data[,i]>0,i]+1))
sd(log2(data[data[,i]>0,i]+1))

cat("sRNAseq_max_all")
i=3
mean(log2(data[data[,i]>0,i]+1))
median(log2(data[data[,i]>0,i]+1))
sd(log2(data[data[,i]>0,i]+1))



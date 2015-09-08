#total mapped reads
reads=read.table("../../DATA/total_mapped_reads",header=T)

#corrlations
out=matrix(ncol=7,nrow=29)
lists=read.table("../../list_final")
a=as.vector(lists[,1])
for (name in c("high_conf","low_conf")) {
n=0
#for (i in c(a,"uniq","all27"))  {
for (i in c("uniq","all27"))  {
	n=n+1
	data=read.table(paste("../known_genes/",i,".overlapped.ws170.targets.",name,".expression",sep=""))
	data=data[data[,11]!='na',]
	if(i=="all27" ) {	data=data[sample(1:nrow(data),1000),] }  # sample all data with same size
	if(i=="all27" || i=="uniq") {
		data[,1]=as.character(data[,1])
		for(j in 1:nrow(data)) {
			title=unlist(strsplit(data[j,1],':'))[1]
			data[j,2]=data[j,2]/reads[reads$TF==title,2]
			data[j,3]=data[j,3]/reads[reads$TF==title,3]
		}	
	}
	ratio=data[,2]/data[,3]
	expression=as.numeric(data[,11])
	pval=data[,4]
	corr1=cor.test(log(ratio),log(expression))
	corr2=cor.test(ratio,expression,method='spearman')
	corr3=cor.test(log(pval),log(expression))
	corr4=cor.test(pval,expression,method='spearman')
	corr5=cor.test(log(data[,2]),log(expression))
	corr6=cor.test(data[,2],expression,method='spearman')
	colnames(out)=c('TF','ratio-pearson','ratio-spearman','pval-pearson','pval-spearman','sig-pearson','sig-spearman')
	out[n,1]=i
	out[n,2]=corr1$estimate	
	out[n,3]=corr2$estimate	
	out[n,4]=corr3$estimate	
	out[n,5]=corr4$estimate	
	out[n,6]=corr5$estimate	
	out[n,7]=corr6$estimate	
}
	write.table(out, file =paste("corr.",name,".out",sep=""), quote=F, sep="\t",col.names=T,row.names=F)
}

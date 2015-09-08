for (l in c("L3","PHA4","HOX","all27")) {
	out_file=paste("shared_targets",l,sep='.')
	write("Shared Targers by differen number of Factors",out_file)
	list=paste("../../../list_",l,sep="")
	list_file=read.table(list)
	prefix=as.vector(list_file[,1])
	for (i in c("No_filter","velcro15","velcro20")) {
	for (j in c("coding","miRNA","nonmiRNA")) {
	for (k in c("high_conf","low_conf","both_conf")){
		out=list();
		share=vector(length=length(prefix))
		for (m in c(prefix,i)) {
			target=paste("../../known_genes/",m,".overlapped.ws170.targets-",j,".",k,".gene_id",sep="")
			if(m != "No_filter") {
			if(file.info(target)$size != 0 ) { 
				data=read.table(target)
				for(gene in data[,1]) {
					out[[gene]]=c(out[[gene]],m)
				}
			}
			}
		}
		# count how many gene in how many factors
		for (gene in out){
			n=length(gene)
			has_velcro=length(gene[gene==i])
			if (has_velcro==0) {
				share[n]=share[n]+1 
			}
		}
		write(paste(l,"-",i,".",j,".",k,sep=""),file=out_file, append=T)	
		write.table(share,file=out_file,append=T)	
	}
	}
	}
}


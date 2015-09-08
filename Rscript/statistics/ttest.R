out=matrix(ncol=2,nrow=32)

for ( j in c("lin35")) {
for ( k in c("uniq")) {
for ( i in c("germline")) {
	inputfile1=paste(i,j,k,sep="_")
	for ( h in c("soma")) {
		inputfile2=paste(i,j,k,sep="_")
		# read file
		input1 = read.table(inputfie1,header=TRUE )
		input2 = read.table(inputfie2,header=TRUE )
		# only high-signal peaks are kept (q val.< 1e-5)	
		qinput1=input1[input1[6]<0.00001,]
		qinput2=input2[input2[6]<0.00001,]
		for (n in c(7:38)) {
			ttest=t.test(qinput1[n],qinput2[n],alternative="two.sided",paired=F)
			out[n,1]=ttest$estimate[1]-ttest$estimate[2]
			out[n,2]=ttest$p.value
		}
		#format output
		colnames(out)=c("diff","p_val")
		rownames(out)=colnames(input1)[7:38]	
		outputfile=paste(i,"-",h,"_",j,"_",k,".ttest",sep="")
		write.table(out, file =outputfile, quote=F, sep="\t",col.names=T,row.names=T)	
	}
}
}
}

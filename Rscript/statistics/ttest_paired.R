coil_prob=0
H_prob=0
E_prob=0
name='NA'
output=matrix(ncol=6,nrow=1)
for (i in dir("prob") ) {
if (substr(i,nchar(i)-4,nchar(i)) == '.prob') {
	rna=read.table(paste("prob/",i,sep=""))
	pro_name=paste (substr(i,1,nchar(i)-4),"pss",sep="")
	pro_name=paste("protein/",pro_name,sep="")
	protein=read.table(pro_name)
	if(nrow(protein)==nrow(rna) & nrow(rna) > 300) {
		data=cbind(protein,rna)
#		for(i in 91:(nrow(rna)-90) ) { data[i,5]=mean(data[(i-25):(i+25),5]) }
		data=data[91:(nrow(rna)-90),] 
		data=data[data[,2]>5,]
		coil_prob=c(coil_prob,mean(data[data[,3]=='C',5]))
		E_prob=c(E_prob,mean(data[data[,3]=='E',5]))
		H_prob=c(H_prob,mean(data[data[,3]=='H',5]))
		name=c(name,i)
	}
}
}
CH=t.test(coil_prob,H_prob,paired=TRUE,alternative="greater")
CE=t.test(coil_prob,E_prob,paired=TRUE,alternative="greater")
EH=t.test(E_prob,H_prob,paired=TRUE,alternative="greater")
j=1
output[j,1]=CH$estimate
output[j,2]=CH$p.value
output[j,3]=CE$estimate
output[j,4]=CE$p.value
output[j,5]=EH$estimate
output[j,6]=EH$p.value
write.table(output, file ="out", quote=F, sep="\t",row.names=F,col.names=F)

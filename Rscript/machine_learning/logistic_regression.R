
Matrix1=read.table("yourfile",head=T,sep="\t");#this is your prepared matrix; should be scaled to 0~1 or -1~1;
classesID=as.vector(unique(Matrix1[,1]));#Matrix1 contains an column (here it is the first column) recording each instances' class label, for example, classA, classB, classC; and classesID gets all these classes label categories, which in this case is ("classA", "classB", "classC");  
Matrix_mean=matrix(NaN,nrow=(ncol(Matrix1)-1),ncol=length(classesID));#this is the logistic regression weights matrix
colnames(Matrix_mean)=classesID;
rownames(Matrix_mean)=colnames(Matrix1[,-1]);
# since logistic regression can only do binary classification, we here use a loop to do "classA(\B\C) to others"
for (i in classesID){
        label=as.vector(Matrix1[,1]);
        label[label==i]=1;
        label[label!=1]=0;
        Matrix2=cbind(label,Matrix1[,-1]);
#logistic regression need a formula in the format of Y~X1+X2+X3... 
        glm.sol<-glm(label ~ GC+blastn+blastx+randfold+sci_max+RNAcode+RNAseq_polyA+RNAseq_total+smRNAseq+tiling_total+X01.H3K18Ac+X02.H3K27me1+X03.H3K27me3+X04.H3K36me2+X05.H3K36me3+X06.H3K4me2+X07.H3K4me3+X08.H3K9Ac+X09.H3K9me2+X10.H3+metCG+metCHG+metCHH,family=binomial,data=Matrix2);
        tmp=glm.sol$coefficients;
        for (j in colnames(Matrix1[,-1])){
                Matrix_mean[j,i]=tmp[j];
        }
}
write.table(Matrix_mean,"logistic_weights",sep="\t",col.names=T,row.names=T,quote=F);



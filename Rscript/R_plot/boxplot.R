#Read the data table
data=read.csv("boxplot_example.csv")
###################
#I.Prepare the data
#1.Normalize the data, etc
for (i in 12:17){
	data[,i]=log(data[,i]+1e-3) # log some expression values
}
for (i in 9:17) {
	maxValue=max(data[,i])  #scale the data into 0-1
	minValue=min(data[,i])
	range=maxValue-minValue
	data[,i]=(data[,i]-minValue)/range
}
data$X8.Identity=data$X8.Identity/100

#2.Make the new matrix for boxplot: cleaning the data table
library("reshape")
m=melt(data[,c(2,7:12,14:17)], id=1)# remove some columns not to show and reshape the matrix into 3 columns for boxplot drawing in bwplot
colnames(m)=c("Type","Feature","Normalized_Value")  #define the new column names

#3.Clean the names of each type and each feature
#Merge sub-types of different elements
m[,1]=sub ("ncRNA_selected","RNAI", m[,1])
m[,1]=sub ("ncRNA_3019","RNAII", m[,1])
m[,1]=sub ("exon_CCDS","CDS", m[,1])
m[,1]=sub ("five_prime_UTR","UTR", m[,1])
m[,1]=sub ("three_prime_UTR","UTR", m[,1])
m[,1]=sub ("ancestral_repeat","AP", m[,1])
#Rename the feature
m[,2]=sub('X7.GC','01.GC Content',m[,2])
m[,2]=sub('X8.Identity','02.DNA Conservation',m[,2])
m[,2]=sub('X9.z_score','03.RNA Struc. Free Energy',m[,2])
m[,2]=sub('X10.SCI','04.RNA Struc. Cons.',m[,2])
m[,2]=sub('X11.tblastx_score','05.Protein Conservation',m[,2])
m[,2]=sub('X12.polyA_RNAseq_MAX','06.PolyA+ RNA-seq',m[,2])
m[,2]=sub('X14.small_RNAseq_MAX','07.Small RNA-seq',m[,2])
m[,2]=sub('X15.Array_totalRNA_MAX','08.Total RNA Array',m[,2])
m[,2]=sub('X16.Array_polyA_MAX','09.PolyA+ RNA Array',m[,2])
m[,2]=sub('X17.Array_nonpolyA_MAX','10.PolyA- RNA Array',m[,2])

###########################
#Making Boxplot
library("lattice")
png("boxplot.png",width=1500,height=500) # pdf is recommended for most cases, or png for figure with huge amount of data points
#pdf("boxplot.pdf") 
attach(m)
bwplot(Normalized_Value ~ Type|Feature,fill=c("green","red","yellow","blue","light blue"),layout=c(10,1))
dev.off()



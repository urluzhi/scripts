

#include "stdafx.h"

//#if !defined(STRUCTURE_CPP)
//#define STRUCTURE_CPP

//#include <stdlib.h>
#include <fstream>
#include <cmath>
#include "structure.h"
#include "platform.h"

using namespace std;





structure::structure(int structures)
{
	allocatedstructures = structures;
	int i;

	
	allocatestructure();

	for (i=1;i<allocatedstructures;i++) {
		energy[i]=0;
	}
	allocated = false;
	
	nnopair=0;
	nopairmax=maxforce-1;
	nopair= new short int[maxforce];//allocate nopair 
	nforbid=0;
	npair=0;
	ndbl=0;
	intermolecular = false;
	ngu = 0;
	templated = false;
	nmod = 0;
	min_gu = 0;
	min_g_or_u = 0;
	nneighbors = 0;
	nregion = 0;
	nmicroarray=0;
	for (i=0;i<maxregions;i++) rnneighbors[i]=0;
	//for (i = 0; i < 100;i++) {
	//	for ( j = 0; j < 25; j++) {
	//		neighbors[i][j]=0;
	//	}
	//}

	stacking = false;
	limitdistance=false;//toogle to true to limit base pairing distance
	maxdistance=600;//default maximum distance between paired nucs
	shaped = false;//be default, a structure does not have SHAPE data associated with it

}

void structure::allocatestructure() {
	int i;
	energy = new int [allocatedstructures];
	ctlabel = new char *[allocatedstructures];
	for (i=0;i<allocatedstructures;i++) ctlabel[i]=new char [ctheaderlength];

}

void structure::deletestructure() {
	int i;

	if (allocated) {
		for (i=0;i<allocatedstructures;i++) {
    		delete[] basepr[i];
   		}

		delete[] basepr;	
	}

	delete[] energy;
	for (i=0;i<allocatedstructures;i++) delete[] ctlabel[i];
	delete[] ctlabel;
}

structure::~structure()
{
	int i;
	delete[] nopair;
	
	deletestructure();
	if (allocated) {
		delete[] numseq;
   		
		delete[] hnumber;
		delete[] nucs;
	}
	if (templated) {
		for (i=0;i<=numofbases;i++) {
    		delete[] tem[i];
   		}

   		delete[] tem;
	}
	if (shaped) delete[] SHAPE;
}



void structure::allocate(int size)

{
   int i;
   numseq = new short int [2*size+1];
   hnumber = new short int [size+1];
   nucs = new char [size+2];
   basepr = new short int *[allocatedstructures];
   for (i=0;i<allocatedstructures;i++) {
    	if (stacking) basepr[i] = new short int [2*size+1];
		else basepr[i] = new short int [size+1];
   }
   allocated = true;

   

}


//this allocates space in an array that is used for folding with phylogenetic data
//	tem == template for allowed and disallowed pairs
void structure::allocatetem()

{
	int i,j;
	//Size = size;//save the size of the array so that the destructor can
   				//deallocate the space

   tem = new bool *[numofbases+1];
   for (i=0;i<=numofbases;i++) {
    	tem[i] = new bool [i+1];
   }
   templated = true;

   //initialize all positions to true:
	for (i=0;i<=numofbases;i++) {
		for (j=i+1;j<=numofbases;j++) {
    		tem[j][i] = true;
		}
   }

}


void structure::checknopair() {
	short int *temp;
	int i;
	if (nnopair >= nopairmax) {
	nopairmax= nopairmax*2;
	temp = new short int [nopairmax];
	for (i=1;i<=nnopair;i++) {
		temp[i]=nopair[i];
	}
	delete [] nopair;
	nopair=temp;
	}
}

void structure::checknumberofstructures() {
	structure *ct2;
	int i,j;

	if (numofstructures>=allocatedstructures-1) {
		//there will be no room for the next structure, double the allocated space for structures
		ct2=new structure (allocatedstructures);
		ct2->stacking = stacking;
		ct2->allocate(numofbases);
		
		for (i=1;i<allocatedstructures;i++) {
			strcpy(ct2->ctlabel[i],ctlabel[i]);
			ct2->energy[i]=energy[i];
			
			for (j=1;j<=numofbases;j++) {
				ct2->basepr[i][j] = basepr[i][j];
				if (stacking) ct2->basepr[i][j+numofbases]=basepr[i][j+numofbases];
			}

		}


		deletestructure();

		allocatedstructures=2*allocatedstructures;

		allocatestructure();
		basepr = new short int *[allocatedstructures];
		for (i=0;i<allocatedstructures;i++) {
    		if (stacking) basepr[i] = new short int [2*numofbases+1];
			else basepr[i] = new short int [numofbases+1];
		}
		

		for (i=1;i<=numofstructures;i++) {
			strcpy(ctlabel[i],ct2->ctlabel[i]);
			energy[i]=ct2->energy[i];
			 
			for (j=1;j<=numofbases;j++) {
				basepr[i][j] = ct2->basepr[i][j];
				if (stacking) basepr[i][j+numofbases] =ct2->basepr[i][j+numofbases];
			}
	
			
		}

		delete ct2;
	}

}

//sort the structures from lowest to highest free energy
void structure::sort() {
	int i;

	for (i=1;i<=numofstructures;i++) basepr[i][0]=(short int) energy[i];

	//take advantage of qsort--
	//store the energies in basepr[i][0] so that they are associated with that array
	
	basepr++;

	qsort((void*) ((basepr)),(size_t) numofstructures,/*(numofbases+1)**/(sizeof(short int *)),ecompare);
	
	basepr--;
	for (i=1;i<=numofstructures;i++) energy[i]=(int)basepr[i][0];

}


//This function reads a SHAPE reactivity datafile and saves the data for a linear penalty.
//calculate (default true) indicate whether these data are being read for folding.  (false means
	//the raw values need to be stored.)
void structure::ReadSHAPE(char *filename, bool calculate) {
	ifstream in;
	int position;
	float data;

	SHAPE = new double [2*numofbases+1];

	in.open(filename);
	shaped = true;

	for (position=0; position <= 2*numofbases; position++) SHAPE[position]=-999.0;

	in >> position;
	in >> data;

	while (!in.eof()) {
		//read and parse all data
			//required format is rows with sequence position followed by reactivity

		if (position<=numofbases) {
			SHAPE[position]=data;
			//SHAPE[position+numofbases]=data;
		}

		in >> position;
		in >> data;
	}
	in.close();

	if (calculate) {
		for (position=0;position<=numofbases;position++) {
			if (SHAPE[position]>0) {
				SHAPE[position]= (log(SHAPE[position]+1.0)*SHAPEslope+SHAPEintercept);
			
			}
			else if (SHAPE[position]>-500) {
				SHAPE[position]= SHAPEintercept;
			
			}
			else {
				SHAPE[position]=0.0;
			
			}
			SHAPE[position+numofbases]=SHAPE[position];
		}
	}
}




//This function reads a SHAPE reactivity datafile and parse the data into single-stranded amd chemical modification constraints.
void structure::ReadSHAPE(char *filename, float SingleStrandThreshold, float ModificationThreshold) {
	ifstream in;
	int position;
	float data;

	in.open(filename);

	in >> position;
	in >> data;

	while (!in.eof()) {
		//read and parse all data
			//required format is rows with sequence position followed by reactivity
	 if (position<=numofbases) {

		if (data>=SingleStrandThreshold) {
			nnopair++;
			nopair[nnopair]=position;
		}
		else if (data>=ModificationThreshold) {
			nmod++;
			mod[nmod]=position;
		}
	  }
		in >> position;
		in >> data;
	}
	in.close();

}


//comparison function used by qsort in structure::sort
int ecompare(const void *i, const void *j) {

	return (**((short int**)i)-**((short int**)j));
	
}


//read a ct file with sequence and structural information
#define linelength 20



void openct(structure *ct,char *ctfile) {
	int count,i,j;
	char base[2],header[ctheaderlength],temp[1000];
	ifstream in;
	in.open(ctfile);
	in >> count;
	j = 0;


	if (count == -100) { //this is a CCT formatted file:
		in >> ct->numofbases;
   
		in >> count;
	

		in.getline(header,ctheaderlength-1);


		ct->allocate(ct->numofbases);
   
		
		for (i=1;i<=ct->numofbases;i++) {
   			in >> ct->numseq[i];
			ct->nucs[i]=*tobase(ct->numseq[i]);
			ct->hnumber[i] = i;
		}
		for (ct->numofstructures=1;ct->numofstructures<=count;ct->numofstructures++) {
			ct->checknumberofstructures();
			strcpy(ct->ctlabel[ct->numofstructures],header);
			in >> ct->energy[ct->numofstructures];
		
			ct->checknumberofstructures();
    		for (i=1;i<=ct->numofbases;i++) {
       			in >> ct->basepr[ct->numofstructures][i];
			}
		}
		ct->numofstructures=count;
		return;

	}
	else {//this is a ct file:
		//first decide if it contains nucleotide stacking data at the far right
		//determine this based on the length of the line:
		in.getline(temp,1000);
		in.getline(temp,1000);
		i = strlen(temp);
		if (i==35) {
			ct->stacking = true;
		}
		//by default, ct->stacking is false

		ct->allocate(count);
		in.close();
		in.open(ctfile);
		in >> ct->numofbases;
		for (ct->numofstructures = 1;!in.eof();(ct->numofstructures)++)	{
			ct->checknumberofstructures();
			strcpy (header,"");
			

			in.getline(header,ctheaderlength-1);


			strcpy((ct->ctlabel[ct->numofstructures]),header);
			for (count=1;count<=((ct->numofbases));count++)	{
		
		
				in >> temp;//ignore base number in ctfile
				in >> base;//read the base
				
				ct->nucs[count]=base[0];
				if (base[0]=='A'||base[0]=='a') ct->numseq[count]=1;
				else if (base[0]=='C'||base[0]=='c') ct->numseq[count]=2;
				else if (base[0]=='G'||base[0]=='g') ct->numseq[count]=3;
				else if (base[0]=='U'||base[0]=='u'||base[0]=='T'||base[0]=='t') ct->numseq[count]=4;
				else if (base[0]=='I') ct->numseq[count]=5;
				else ct->numseq[count]=0;

				if (ct->numseq[count]==5&&ct->numofstructures==1) {
      				ct->intermolecular = true;
       				ct->inter[j] = count;
					j++;
				}
				in >> temp;//ignore numbering
				in >> temp;//ignore numbering
				in >> ct->basepr[ct->numofstructures][count];//read base pairing info
				in >> ct->hnumber[count];//read historical numbering
				if (ct->stacking) in >> ct->basepr[ct->numofstructures][count+ct->numofbases];
			}
			in >> ct->numofbases; //start on next structure and see whether the end of file is reached
		}
	}
	(ct->numofstructures)--;
	return;
}


//takes a nucleotide input and stores a numeric representation
void tonum(char *base,structure *ct,int count)	{
if (!strcmp(base,"A")) (ct->numseq[count] = 1);
else if(!strcmp(base,"B")) {
	(ct->numseq[count] = 1);

}
else if(!strcmp(base,"a")) {
	ct->numseq[count]=1;
	ct->nnopair++;
	ct->checknopair();
	ct->nopair[ct->nnopair] = count;
}
else if(!strcmp(base,"C")) (ct->numseq[count] = 2);
else if(!strcmp(base,"Z")) {
	(ct->numseq[count] = 2);

}
else if(!strcmp(base,"c")) {
	ct->numseq[count] = 2;
	ct->nnopair++;
	ct->checknopair();
	ct->nopair[ct->nnopair] = count;
}
else if(!strcmp(base,"G")) (ct->numseq[count] = 3);
else if(!strcmp(base,"H")) {
	(ct->numseq[count] = 3);

}
else if(!strcmp(base,"g")) {
	ct->numseq[count] = 3;
	ct->nnopair++;
	ct->checknopair();
	ct->nopair[ct->nnopair] = count;
}

else if(!strcmp(base,"U")||!strcmp(base,"T")) (ct->numseq[count] = 4);
else if(!strcmp(base,"V")||!strcmp(base,"W")) {
	(ct->numseq[count] = 4);

}
else if(!strcmp(base,"u")||!strcmp(base,"t")) {
	ct->numseq[count] = 4;
	ct->nnopair++;
	ct->checknopair();
	ct->nopair[ct->nnopair] = count;
}

else if(!strcmp(base,"I")) {
	ct->numseq[count] = 5;
	ct->intermolecular= true;
}

else (ct->numseq[count]=0);  //this is for others, like X
return;
}

short int tonumi(char *base)	{
short int	a;
if (!strcmp(base,"A")||!strcmp(base,"B")) (a = 1);
else if(!strcmp(base,"C")||!strcmp(base,"Z")) (a = 2);
else if(!strcmp(base,"G")||!strcmp(base,"H")) (a = 3);
else if(!strcmp(base,"U")||!strcmp(base,"V")) (a = 4);
else if(!strcmp(base,"T")||!strcmp(base,"W")) (a = 4);
else (a=0);  //this is for others, like X
return a;
}


#define seqline 300 //maximum line length in .seq file
//returns 0 on error
int openseq (structure *ct,char *seqfile) {
char temp[seqline],seq[seqline],base[seqline],test[seqline];
int i,j,length,nucs;

ct->nnopair = 0;

//nucs = 0;

FILE *se;


//read the sequence file to get the number of nucleotides
se=fopen(seqfile,"r");

do {
	fgets(temp,seqline,se);
	strncpy(test,temp,1);
	strcpy(test+1,"\0");
} while (!strcmp(test,";"));

//fgets(ct->ctlabel[1],seqline,se);
strcpy(ct->ctlabel[1],temp);

nucs = 1;

while (1) {
	fgets(seq,seqline,se);
	length = strlen (seq);
	for (j=0;j<length;j++) {//up to length-1 bcs of /n character
		strncpy (base,seq+j,1);
		strcpy (base+1,"\0");
		if (!strcmp(base,"1")) break;
      if (!(!strcmp(base,"A")||!strcmp(base,"a")||!strcmp(base,"C")||!strcmp(base,"c")
      	||!strcmp(base,"G")||!strcmp(base,"g")||!strcmp(base,"T")
         ||!strcmp(base,"t")||!strcmp(base,"U")||!strcmp(base,"u")
         ||!strcmp(base,"X")||!strcmp(base,"x")||!strcmp(base," ")||
         !strcmp(base,"\n")||!strcmp(base,"N"))) return 0;
	  
	 
//		tonum(base,ct,(i));
      if (strcmp(base," ")&&strcmp(base,"\n")) nucs++;
	}
	if (!strcmp(base,"1")) break;
}
ct->numofbases = nucs - 1;
nucs--;
fclose (se);

if (nucs==0) return 0;

ct->allocate(nucs);



//now read the file
se=fopen(seqfile,"r");

do {
	fgets(temp,seqline,se);
	strncpy(test,temp,1);
	strcpy(test+1,"\0");
} while (!strcmp(test,";"));

//fgets(ct->ctlabel[1],seqline,se);
strcpy(ct->ctlabel[1],temp);

i = 1;
while (1) {
	fgets(seq,seqline,se);
	length = strlen (seq);
	for (j=0;j<length;j++) {//up to length-1 bcs of /n character
		strncpy (base,seq+j,1);
		strcpy (base+1,"\0");
		if (!strcmp(base,"1")) break;
      if (!(!strcmp(base,"A")||!strcmp(base,"a")||!strcmp(base,"C")||!strcmp(base,"c")
      	||!strcmp(base,"G")||!strcmp(base,"g")||!strcmp(base,"T")
         ||!strcmp(base,"t")||!strcmp(base,"U")||!strcmp(base,"u")
         ||!strcmp(base,"X")||!strcmp(base,"x")||!strcmp(base," ")||
         !strcmp(base,"\n")||!strcmp(base,"N"))) return 0;
		
	  if ((!strcmp(base,"A")||!strcmp(base,"a")||!strcmp(base,"C")||!strcmp(base,"c")
      	||!strcmp(base,"G")||!strcmp(base,"g")||!strcmp(base,"T")
         ||!strcmp(base,"t")||!strcmp(base,"U")||!strcmp(base,"u")
         ||!strcmp(base,"X")||!strcmp(base,"x")||!strcmp(base,"N"))) {
	  
	  
			tonum(base,ct,(i));
			ct->nucs[i]=base[0];
			ct->hnumber[i] = i;
			i++;
	  }
	}
	if (!strcmp(base,"1")) break;
}
//ct->numofbases = i - 1;

fclose (se);




return 1;
}

//outputs a ct file (connection table)
void ctout (structure *ct,char *ctoutfile) {
	int count,i;//length
	char line[2*ctheaderlength],number[2*numlen];//base[2]

	FILE *ctfile;
	ctfile=fopen(ctoutfile,"w");

	for (count=1;count<=(ct->numofstructures);count++) {

		strcpy(line,"");
		sprintf(line,"%5i",ct->numofbases);


	
		sgifix  //this corrects a difference on the sgi computers
		if (ct->energy[count]!=0) {
   			strcat(line,"  ENERGY = ");

			//gcvt((float (ct->energy[count]))/conversionfactor,6,number);
			if (conversionfactor==10)
				sprintf(number,"%.1f",(float (ct->energy[count]))/conversionfactor);
			else if (conversionfactor==100)
				sprintf(number,"%.2f",(float (ct->energy[count]))/conversionfactor);
			else sprintf(number,"%f",(float (ct->energy[count]))/conversionfactor);
	
   			strcat(line,number);
   			strcat(line,"  ");
		}
		else strcat(line,"  ");
		strcat(line,ct->ctlabel[count]);
		fputs (line,ctfile);
		for (i=1;i<ct->numofbases;i++) {
			if (ct->stacking) {
				sprintf(line,"%5i%2c%8i%5i%5i%5i%5i\n",
					i,ct->nucs[i],(i-1),(i+1),ct->basepr[count][i],ct->hnumber[i],ct->basepr[count][i+ct->numofbases]);
			}
			else {
				sprintf(line,"%5i%2c%8i%5i%5i%5i\n",
					i,ct->nucs[i],(i-1),(i+1),ct->basepr[count][i],ct->hnumber[i]);
			}
			fputs(line,ctfile);
		}
		
   
		//last nucleotide not connected--
		i = ct->numofbases;
		if (ct->stacking) {
			sprintf(line,"%5i %1c%8i%5i%5i%5i%5i\n",
				i,ct->nucs[i],(i-1),0,ct->basepr[count][i],ct->hnumber[i],ct->basepr[count][i+ct->numofbases]);
		}
		else {
			sprintf(line,"%5i %1c%8i%5i%5i%5i\n",
				i,ct->nucs[i],(i-1),0,ct->basepr[count][i],ct->hnumber[i]);
		}
		fputs(line,ctfile);

	}

	fclose (ctfile);
	return;
}


char *tobase (int i)

{  //function is a look up table to convert the base
	// 	integer represention to a familiar character

	if (i==1) return "A";
	 else if (i==2) return "C";
	 else if (i==3) return "G";
	 else if (i==4) return "U";
	 else if (i==0) return "X";
    else if (i==5) return "I";
	 else return "?";   //to track down mistakes

}

void sortstructures (structure *ct) {//this function reorders the structures by
													//the efn energies
register int c;
int cur,i;
char ctheader[ctheaderlength];

for (c = 2; c<=(ct->numofstructures);c++){

	cur = c;

	while (cur>1) {
		if ((ct->energy[cur])<(ct->energy[cur-1])) {
      	swap(&ct->energy[cur],&ct->energy[cur-1]);
         //also swap the ct labels:
         strcpy(ctheader, ct->ctlabel[cur]);
         strcpy(ct->ctlabel[cur],ct->ctlabel[cur-1]);
         strcpy(ct->ctlabel[cur-1],ctheader);
         for (i=1;i<=(ct->numofbases);i++) {
         	swap(&ct->basepr[cur][i],&ct->basepr[cur-1][i]);
         }
         cur--;
   	}
		else {
		break;
		}
	}
}

}

void swap(int *a,int *b) {
int temp;

	temp = *a;
   *a = *b;
   *b = temp;

return;

}

void swap (short int *a,short int *b) {
short int temp;

	temp = *a;
   *a = *b;
   *b = temp;

return;

}

void swap(float *a,float *b) {
float temp;

	temp = *a;
   *a = *b;
   *b = temp;

return;

}

//#endif

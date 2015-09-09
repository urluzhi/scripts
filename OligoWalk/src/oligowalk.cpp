/*====================================================================================================

oligowalk.cpp : This is the file including main function
oligowalk_test: This is the program to explore the database finding correlation between inhibition efficacy and energy

OligoWalk calculate the binding affinity of oligomer with strucutured RNA
They are revised based on Mathews' code from RNAStructure

oligomer (such as siRNA) is L in lengther
optimal and sub-optimal structures of target can be calculated
partition calculation and stochastic smapling method can also be 
used to predict the breaking energy of target and oligo structure

intermolecular.cpp: handle the options of folding and binding
pclass.cpp: includes the important classed(Pclass for normal partion function, OligoPclass for refilling 
of constrained sequence and scan folding at different site				
Created:Nov. 2005						Modified: 
zhi_lu@urmc.rochester.edu
=======================================================================================================*/
#include <iostream>
#include "stdafx.h"
#include "globals.h"
#include "intermolecular.h"


#define seqline 300
//maximum line length in .seq file
//returns 0 on error
int openseq_frac (structure *ct,char *seqfile, int start, int end); //open one fraction of a sequence

void getdh(char *loop, char *stackf, char *tstackh, char *tstacki,char *tloop, char *miscloop, char *danglef, char *int22,char *int21,char *coax, char *tstackcoax,char *coaxstack, char *tstack, char *tstackm,char *triloop, char *int11,char *hexaloop,char *tstacki23, char *tstacki1n ,char *datapath);


int main(int argc, char *argv[]) {
	
	char stackf[maxfil], tstackh[maxfil], tstacki[maxfil],
 		 tloop[maxfil], miscloop[maxfil], danglef[maxfil], int22[maxfil],
		 int21[maxfil],coax[maxfil], tstackcoax[maxfil],coaxstack[maxfil],
		 tstack[maxfil], tstackm[maxfil],triloop[maxfil],int11[maxfil],
		 loop2[maxfil],hexaloop[maxfil],
		 tstacki23[maxfil], tstacki1n[maxfil];
	char seqfilename[100], reportfilename[100];
	int i,j,M,N,start,end,foldsize,unit;
	char shapefile[maxfil];
	shapefile[0]='\0';	//default: shape is not used
	int distance=0;
	int mode,suboptimal,length,conc, useprefilter;
	int **table,**numofsubstructures;
	int *TEST,TESTnum=-1;
	int *TESTon;
	double c;
	int isdna=0;
	bool scoreit=false;
	bool WRITE=false;
	structure *ct;
	datatable data,dhdata;
	datatable *ddata;
	rddata *hybriddata;
	thermo *helixstack;
	//TProgressDialog PD;
	siPREFILTER *prefilter;
	//define the datapath from shell
	char *datapath;
	datapath=getenv("DATAPATH");
	if (datapath == NULL) {
		datapath="../../NN_data/";
	}
	//default oligo tested number
	TESTon=new int[5];
	
/*--------------------------------------------------------------------------------------------------
isdna 0 - rna(oligo)-rna hybridyzation
isdna 1 - dna(oligo)-rna hybridyzation
isdna 2 - dna(oligo)-dna hybridyzation

mode 1 - break local target structure to bind oligo 
mode 2 - refold target RNA after oligo binding
mode 3 - no target structure considered

suboptimal 0 - only consider optimal structure
suboptimal 1 - like choice 3,using suboptimal structures,but the whole set from alltarce() function prediction
suboptimal 2 - using partition funcion considering every possible structure  of target
			   suboptimal 2 can only used with mold 2
suboptimal 3 - using suboptimal structures (heuristic method) for both oligo-free and oligo-bound target RNA
suboptimal 4 - using stochasit sampling method to sample 1000 structures

useprefileter 1 - using creteria to prefill functional siRNA; (-- you may not want to type -test )
foldsize >0  - only folding a fragment with size=foldsize+binding length, 
			   which is centered on the siRNA binding region
			   when foldsize>1, only option 2 plus Usesub 2 is the availabe option
test 0         testing all sites without prefilter score
test 3 2 5 8   testing 3 sites 2,5 and 8 without prefilter score -- do not type -test if you want prefiter score
----------------------------------------------------------------------------------------------------*/

	if (argc > 1) {
	

//get the input information	
		i=1;
		while (i<argc) {
			if (!strcmp(argv[i],"-seq")) {//input sequnce file
				strcpy(seqfilename,argv[i+1]);
				i=i+2;

			}
			else if (!strcmp(argv[i],"-o")) {//output file stored the result in debug mode
				strcpy(reportfilename,argv[i+1]);
				i=i+2;

			}
			else if (!strcmp(argv[i],"-m")) {//choice of mode, see detatil at the beginning of main function
				mode = atoi(argv[i+1]);
				i=i+2;

			}
			else if (!strcmp(argv[i],"-st")){//start position of folding region of target
				start = atoi(argv[i+1]);
				i+=2;
			
			}
			else if (!strcmp(argv[i],"-en")){
				end = atoi(argv[i+1]);
				i+=2;
			
			}
			else if (!strcmp(argv[i],"-M")){//start posion of scanning on folded region of target
				M = atoi(argv[i+1]);
				i+=2;
			
			}
			else if (!strcmp(argv[i],"-N")){
				N = atoi(argv[i+1]);
				i+=2;
			
			}
			else if (!strcmp(argv[i],"-l")) {//length of hybridiztion
				length = atoi(argv[i+1]);
				i=i+2;

			}
			else if (!strcmp(argv[i],"-type")) {//if it is a DNA-RNA or RNA-RNA(default) or DNA-DNA
				if(!strcmp(argv[i+1],"d"))		{	isdna=1;	} //DNAoligo-RNA
				else if (!strcmp(argv[i+1],"dd")) {	isdna=2;	} //DNAoligo-DNA
				i+=2;

			}
			else if (!strcmp(argv[i],"-s")) {//suboptimal structure choices
				suboptimal = atoi(argv[i+1]);
				i+=2;

			}
			else if (!strcmp(argv[i],"-co")) {//concentration of oligo
				conc = atoi(argv[i+1]);
				i=i+2;
			}
			else if (!strcmp(argv[i],"-unit")) {//concentration of oligo
				unit = atoi(argv[i+1]);
				i=i+2;
			}
			else if (!strcmp(argv[i],"-fi")) {//prefilter of siRNA
				useprefilter = atoi(argv[i+1]);
				i=i+2;
			}
			else if (!strcmp(argv[i],"-score")) {//prefilter of siRNA
				scoreit = true;
				i=i+1;
			}
			else if (!strcmp(argv[i],"-fold")) {//fold substructure size centered on binding
				foldsize = atoi(argv[i+1]);
				i=i+2;
			}
			else if (!strcmp(argv[i],"-dist")) {//fold substructure size centered on binding
				distance = atoi(argv[i+1]);
				i=i+2;
			}
			else if (!strcmp(argv[i],"-shape")) {//fold substructure size centered on binding
				strcpy(shapefile,argv[i+1]);
				i=i+2;
			}
			else if (!strcmp(argv[i],"-test")) {//fold substructure size centered on binding
				TESTnum = atoi(argv[i+1]);
				delete [] TESTon;
				TESTon=new int[TESTnum+1]; 
				for (j=1;j<=TESTnum;j++)		TESTon[j]=atoi(argv[i+1+j]);
				i=i+2+TESTnum;
			}
			else if (!strcmp(argv[i],"-write")) {//write sav files to save time in test mode
				WRITE= true;
			}
			else {
				cout<<"unrecognized modifier:!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
				cout<< argv[i]<<"\n";
				cout<<"Please type command only to see an example of input.\n";
				return 1;
				break;
			}
		}

		c = (double)conc * pow(10,((double) unit));
		if (suboptimal==2 && mode != 2) {
			cout << "Partition function (-s 2) can only used with refolding target (-m 2).\n";
			return 1;
		}
	
		//open the sequence into ct 
		ct=new structure;	
		openseq_frac(ct,seqfilename, start, end);
		
		if (TESTnum==-1)	scoreit=true;
		else				scoreit=false;
		TEST = new int[ct->numofbases+1];
		for (i=1;i<=ct->numofbases;i++)	{
			if (TESTnum==0)		TEST[i]=1;
			else				TEST[i]=0;
			for (j=1;j<=TESTnum;j++) {
				if (i==TESTon[j])	TEST[i]=1;
			}
		}
		//allocate memory of energy tables
		table = new int*[ct->numofbases - length + 2];
		for (i = 0; i < ct->numofbases - length + 2; i++) {
   			table[i] = new int[6];
		}
		//allocate memory of number of suboptimal structures
		numofsubstructures= new int*[ct->numofbases - length +2];
		for (i = 0; i < ct->numofbases - length + 2; i++)	{
			numofsubstructures[i]= new int [2];
			numofsubstructures[i][0]=0;
			numofsubstructures[i][1]=0;
		}
		
		//show the inputed information
//		cout <<"Content-Type: text/html\n\n";
//		cout << "<html><head><title>OligoWalk Calculation Result</title></head>\n";
//		cout << "<body bgcolor=\"#AFCAE0\">\n";	
		cout << "<h1>RNAstructure OligoWalk calculation:</h1>\n\n";
		cout<<"<HR>\n";
	
		cout << "Target sequence: "<<ct->ctlabel[1] <<"<br>";
		cout << "Total size of the target -> "<< ct->numofbases<<" nucleotides <br>\n";
		if(N==0) cout<<"Scanned position on target: "<<M<<" (nt) -- end <br>\n";
		else	 cout << "Scanned position on target: "<<M<<" (nt) -- "<<N<<" (nt) <br>\n";
	
		cout <<"\n<br>Oligonucleotides:<br>\n";
		if (isdna == 0) cout << "Type -> RNAoligo-RNA ";
		else if (isdna == 1)	   cout << "Type -> DNAoligo-RNA ";
		else if (isdna ==2)	cout <<"Type -> DNAoligo-DNA ";
		cout << "	length -> "<<length;
		cout << "	concentration: -> "<<c<<"  M<br><br>\n\n";
		
		cout << "Method options: <br>\n";
		if (mode ==1) cout << "Break local structure of the target to bind oligo<br>\n";
		else if (mode == 2) cout << "Refold target structure after oligo binding<br>\n";
		else if (mode == 3) cout<< "Not considering the structure of target<br>\n";
		if (mode ==2 || mode==1) {	
			if (foldsize!=0)	cout << "Folding region size -> "<<foldsize+length<<" nt <br>\n";
			else 			cout << "Folding region size -> Global target region<br>\n";
		}
		if (suboptimal ==1) {
			cout << "All Suboptimal structures within a energy difference window were considered .<br>\n";
			cout << ".<br>\n";
		}
		else if (suboptimal == 0) {
			cout << "Only one optimal structure was considered.<br>\n";
		}
		else if (suboptimal == 2){
			cout<< "All possible structures were considered using Partition Function.<br>\n";
		}
		else if (suboptimal == 3){
			cout<< "Heuristic suboptimal structures were considered for both oligo-free and olig-bound target.<br>\n";
			cout<< "There are at most 1000 (free energy is less than 10% off lowest one) suboptimal structures.<br>\n";
		}
		else if (suboptimal == 4) {
			cout<< "1000 stochastic sampled structures were considered for both oligo-free and olig-bound target.<br>\n";
		}
		else	cout<< " Unrecognized modifier for -s (Suboptimal structures).<br>\n";
		if (distance >0)	cout <<"The base pairs with distance larger than "<<distance<<" nt are not allowed.<br>\n";
		if (useprefilter != 0)  cout<<"\n<br>Prefilter was used.<br>\n";	
		if (shapefile[0] != '\0')  cout<<"\n<br>The shape information was used: "<<shapefile<<" <br>\n\n";
		cout<<"<HR>\n\n";
		
		//read the data talbes
		//Prepare for DNA data
		if (isdna != 0) 	ddata = new datatable();
		//Prepare for DNA-RNA hybrid data
		if (isdna ==1 ) hybriddata = new rddata;
		//Prepare helixstack for all (default: RNA)		
		helixstack = new thermo(datapath);
		//Read RNA data
		if (isdna <= 1) {
			getdat (loop2,stackf,tstackh,tstacki,tloop,miscloop,danglef,
   					int22,int21,coax,tstackcoax,coaxstack,tstack,tstackm,triloop,int11,
					hexaloop,tstacki23, tstacki1n,datapath,true);
			if (opendat (loop2,stackf,tstackh,tstacki,tloop,miscloop,danglef,int22,int21,
   						coax,tstackcoax,coaxstack,tstack, tstackm, triloop, int11,hexaloop,
						tstacki23, tstacki1n, &data)==0) {
				cout << "A Thermodynamic Data File was NOT Found.  Abort.<br>\n";	
				return 1;
			}
			//open the enthalpy files
			getdh (loop2,stackf,tstackh,tstacki,tloop,miscloop,danglef,
   					int22,int21,coax,tstackcoax,coaxstack,tstack,tstackm,triloop,int11,
					hexaloop,tstacki23, tstacki1n,datapath);
			if (opendat (loop2,stackf,tstackh,tstacki,tloop,miscloop,danglef,int22,int21,
   						coax,tstackcoax,coaxstack,tstack, tstackm, triloop, int11,hexaloop,
						tstacki23, tstacki1n, &dhdata)==0) {
				cout << "A Thermodynamic Data File (.dh) was NOT Found.  Abort.<br>\n";	
				return 1;
			}
		}
		//Read DNA data
		if (isdna != 0 ) {		
			getdat (loop2,stackf,tstackh,tstacki,tloop,miscloop,danglef,
   	      		int22,int21,coax,tstackcoax,coaxstack,tstack,tstackm,triloop,int11,
				hexaloop,tstacki23, tstacki1n,datapath,false);
			if (opendat (loop2,stackf,tstackh,tstacki,tloop,miscloop,danglef,int22,int21,
   	      		coax,tstackcoax,coaxstack,tstack,tstackm,triloop,int11,hexaloop,
				tstacki23, tstacki1n,ddata)==0) {

				cout << "A Thermodynamic Data File was NOT Found.  Abort.<br>\n";	
				return 1;
			}
		}
		//Read DNA-RNA hybrid data
		if (isdna == 1 ) {
			strcpy(stackf,datapath);
			strcat (stackf,"stackdr.dat");
			readrd (hybriddata,stackf);

			strcpy(helixstack->DH,datapath);
			strcat(helixstack->DH,"stackdr.dh");

			strcpy(helixstack->DS,datapath);
			strcat(helixstack->DS,"stackdr.ds");

			strcpy(helixstack->HELIX,datapath);
			strcat(helixstack->HELIX,"helixdr.dat");
		}
		//Read DNA-DNA data (HELIX,DH and DS)
		else if (isdna == 2) {
			//no dh and ds data for DNA-DNA yet; so the default (RNA-RNA) is used in class helixstack 
		}	
		//Read helixstack for all
		helixstack->read();

//define prefilter
		prefilter = new siPREFILTER(data,dhdata,useprefilter,scoreit,ct->numofbases - length + 2,isdna);

//oligowalk
		if (N==0) N=ct->numofbases-length+1;
		else N=N-length+1; 
		if ( isdna <=1) {
			olig(isdna, mode, ct, length, 
				c, table,numofsubstructures, data, *ddata, 
				hybriddata,suboptimal,NULL,helixstack,
				M,N,prefilter,foldsize,distance,shapefile,TEST,WRITE);
		}
		else if (isdna ==2) { //DNA-DNA replace data with *ddata
			olig(isdna, mode, ct, length, 
				c, table,numofsubstructures, *ddata, *ddata, 
				hybriddata,suboptimal,NULL,helixstack,
				M,N,prefilter,foldsize,distance,shapefile,TEST,WRITE);
		}

		report(reportfilename,ct,table,numofsubstructures, 
			length,isdna, c,suboptimal,M,N,prefilter,foldsize);

//delete the dynamic parameters
		for (i = 0; i < ct->numofbases - length + 2; i++) {
   			delete[] table[i];
			delete[] numofsubstructures[i];
		}
		delete[] table;
		delete[] numofsubstructures;
		delete[] TEST;	
		delete[] TESTon;	
		
		if (isdna != 0) delete ddata;
		if (isdna ==1)	delete hybriddata;
				
		delete helixstack;
		delete ct;
		delete prefilter;
		
//		cout << "RNAstructure Calculation Complete\n";

	}

	else {
		//The parameters in command line were wrong, explain that to user
		cout <<"Usage oligowalk [-modifiers]\nSee Online Help for more Information.";
		cout <<"Input format:\n";
		cout <<"oligowalk -type r -seq foo.seq -o foo.out -m 1 -l 19 -s 1 -co -3\n";

	}
	

	return 0;
}





//Open a fraction of a sequence

int openseq_frac (structure *ct,char *seqfile,int start,int end) {
char temp[seqline],seq[seqline],base[seqline],test[seqline];
int i,j,length,nucs;
char *base_str;
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

if (end > ct->numofbases){
		cout<<"Sequence is shorter than "<<end<<"\n";
		return 0;
		}
if (end==0)	 end = ct->numofbases;
ct->numofbases = end-start+1;
	


if (nucs==0) return 0;

ct->allocate(ct->numofbases);


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
base_str=new char[nucs+2];
base_str[0]='X';
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
	  
			base_str[i]=base[0];	  
			i++;
		 }
	}
	if (!strcmp(base,"1")) break;
}
//ct->numofbases = i - 1;
	base_str[i]='\0';
	i=1;
	for (j=start;j<=end;j++) {

			base[0]=base_str[j];
			strcpy (base+1,"\0");
			tonum(base,ct,(i));
			ct->nucs[i]=base[0];
			ct->hnumber[i] = i;
			i++;
	}			
delete base_str;
fclose (se);




return 1;
}

/*      Function getdh      
Function gets the names of dhdata files to open */

void getdh(char *loop, char *stackf, char *tstackh, char *tstacki,char *tloop, char *miscloop, char *danglef, char *int22,char *int21,char *coax, char *tstackcoax,char *coaxstack, char *tstack, char *tstackm,char *triloop, char *int11,char *hexaloop,char *tstacki23, char *tstacki1n,char *datapath )
{
strcpy (loop,datapath);

strcpy (stackf,datapath);

strcpy (tstackh,datapath);

strcpy (tstacki,datapath);

strcpy (tloop,datapath);

strcpy (miscloop,datapath);

strcpy (danglef,datapath);

strcpy (int22,datapath);

strcpy (int21,datapath);

strcpy (triloop,datapath);

strcpy (coax,datapath);

strcpy (tstackcoax,datapath);

strcpy (coaxstack,datapath);

strcpy (tstack,datapath);

strcpy (tstackm,datapath);

strcpy (int11,datapath);

strcpy (hexaloop,datapath);

strcpy (tstacki23,datapath);

strcpy (tstacki1n,datapath);

strcat (loop,"loop.dh");
strcat (stackf,"stack.dh");
strcat (tstackh,"tstackh.dh");
strcat (tstacki,"tstacki.dh");
strcat (tloop,"tloop.dh");
strcat (miscloop,"miscloop.dh");
strcat (danglef,"dangle.dh");
strcat (int22,"int22.dh");
strcat (int21,"int21.dh");
strcat (coax,"coaxial.dh");
strcat (tstackcoax,"tstackcoax.dh");
strcat (coaxstack,"coaxstack.dh");
strcat (tstack,"tstack.dh");
strcat (tstackm,"tstackm.dh");
strcat (triloop,"triloop.dh");
strcat (int11,"int11.dh");
strcat (hexaloop,"hexaloop.dh");
strcat (tstacki23,"tstacki23.dh");
strcat (tstacki1n,"tstacki1n.dh");
}


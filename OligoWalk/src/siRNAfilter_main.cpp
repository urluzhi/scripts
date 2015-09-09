/*====================================================================================================

siRNAfilter_main.cpp : This is the file including main function
filter non-functional siRNA with rules list in Reynolds, Nautre Biotech. 2004

input file fomrat: ( antisense siRNA sequence)

AAGGUUU 
GGCCUUU
...
end


Created:July 14. 2006						Modified: July. 14, 2006
zhi_lu@urmc.rochester.edu
=======================================================================================================*/
#include "stdafx.h"
#include "structure.cpp"
#include "algorithm.cpp"
#include "rna_library.cpp"
#include "globals.cpp"
#include "thermo.cpp"
#include "siRNAfilter.cpp"
#include <math.h>
#include <string.h>
#define datapath "../../NN_data/"

struct rddata {

   short int stack[5][5][5][5];
   short int init;
};
int readrd (rddata* data,char* dnarna) ;
void getdh(char *loop, char *stackf, char *tstackh, char *tstacki,char *tloop, char *miscloop, char *danglef, char *int22,char *int21,char *coax, char *tstackcoax,char *coaxstack, char *tstack, char *tstackm,char *triloop, char *int11,char *hexaloop,char *tstacki23, char *tstacki1n );


int main(int argc, char *argv[]) {

	bool isdna=false;
	int i;
	char *seqfile=argv[1]; //store the sequence file name
	char seq[500];  //store the sequence 
	char base[2];
	int numofbases;
	structure *ct;
	siPREFILTER *prefilter;
	datatable data,dhdata;
	datatable *ddata;
	rddata *hybriddata;
	thermo *helixstack;
	
	//begin to read the data tables
	////////////////////////////////////////////////
	char stackf[maxfil], tstackh[maxfil], tstacki[maxfil],
 		 tloop[maxfil], miscloop[maxfil], danglef[maxfil], int22[maxfil],
		 int21[maxfil],coax[maxfil], tstackcoax[maxfil],coaxstack[maxfil],
		 tstack[maxfil], tstackm[maxfil],triloop[maxfil],int11[maxfil],
		 loop[maxfil],hexaloop[maxfil],
		 tstacki23[maxfil], tstacki1n[maxfil];
	helixstack = new thermo(datapath);

	//open free energy parameters files
	getdat (loop,stackf,tstackh,tstacki,tloop,miscloop,danglef,
   			int22,int21,coax,tstackcoax,coaxstack,tstack,tstackm,triloop,int11,
			hexaloop,tstacki23, tstacki1n,datapath,true);
	if (opendat (loop,stackf,tstackh,tstacki,tloop,miscloop,danglef,int22,int21,
   			coax,tstackcoax,coaxstack,tstack, tstackm, triloop, int11,hexaloop,
			tstacki23, tstacki1n, &data)==0) {
		cout << "A Thermodynamic Data File was NOT Found.  Abort.<br>\n";	
		return 1;
	}
	//open the enthalpy files
	getdh (loop,stackf,tstackh,tstacki,tloop,miscloop,danglef,
           int22,int21,coax,tstackcoax,coaxstack,tstack,tstackm,triloop,int11,hexaloop,tstacki23, tstacki1n);
	if (opendat (loop,stackf,tstackh,tstacki,tloop,miscloop,danglef,int22,int21,
        coax,tstackcoax,coaxstack,tstack,tstackm,triloop,int11,hexaloop,tstacki23, tstacki1n,&dhdata)==0) {
			cout << "A dhdata file was lost";
			cin >> i;
	}



	//read dna parameters

	if (isdna) {
		ddata = new datatable();
		hybriddata = new rddata;
		getdat (loop,stackf,tstackh,tstacki,tloop,miscloop,danglef,
   	      		int22,int21,coax,tstackcoax,coaxstack,tstack,tstackm,triloop,int11,
				hexaloop,tstacki23, tstacki1n,datapath,false);
		if (opendat (loop,stackf,tstackh,tstacki,tloop,miscloop,danglef,int22,int21,
   	      		coax,tstackcoax,coaxstack,tstack,tstackm,triloop,int11,hexaloop,
				tstacki23, tstacki1n,ddata)==0) {
			cout << "A Thermodynamic Data File  was NOT Found.  Abort.<br>\n";	
			return 1;
		}
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
	helixstack->read();
	//////////////////////////////////////////////////
	//finish reading datatables
	
	//read sequence from seqfile
	ifstream in(seqfile);
	in >> seq;
	ofstream out(argv[2]);
	prefilter = new siPREFILTER(data,dhdata,1,1,1,isdna); //1(using filter) 1(use score) 1(one oligo) 
//	out << "score \t melting \t end_diff \n";
	while (strcmp(seq,"end") !=0 ) {
		//define the structure of the oligo
		numofbases=strlen(seq);
		ct=new structure;
		ct->allocate(numofbases+1);
		strcpy(ct->ctlabel[1],"oligo");
		ct->numofbases = numofbases;
		ct->intermolecular = false;
		for (i=1;i<=numofbases;i++) {
			base[0]=seq[i-1];
			base[1]='\0';
			tonum(base,ct,i);
		}
		//run prefilter to get the score
		prefilter->count(ct,1,0);//1 (first oligo), 0 (not the using test option in oligo())
		//output the final score
		//cout << prefilter->enddiff[1]<<"\n";
		out << prefilter->score[1]<<"\n";
//		out << prefilter->melt[1]<<"\t";
//		out << prefilter->enddiff[1]<<"\n";
		//cout << prefilter->melt[1]<<"\n";
		delete ct;
		in>>seq;
	}
	in.close();
	delete prefilter;
	out.close();

	if (isdna) {
		delete ddata;
		delete hybriddata;
	}
	delete helixstack;
	return 0;
}



//=======================================================================
int readrd (rddata* data,char* dnarna) {
	
	int count,i,k,j,l;
	ifstream dr;
	char lineoftext[100];

	//make sure the file exists
	FILE *check;
	if ((check = fopen(dnarna, "r"))== NULL) {
	return 0;
	}

	fclose(check);
	dr.open(dnarna);
	/* Read info from stackdr */
	//add to the stack table the case where X (represented as 0) is looked up:
	for (count=1;count<=2;count++) dr >> lineoftext;//get past text in file
	dr >> lineoftext;
	data->init =(int)floor(conversionfactor*(atof(lineoftext)));

	for (count=1;count<=42;count++) dr >> lineoftext;//get past text in file
	for (i=0;i<=4;i++) {
		if (i!=0) for (count=1;count<=60;count++) dr >> lineoftext;
		for (k=0;k<=4;k++) {
			for (j=0;j<=4;j++) {
				for (l=0;l<=4;l++) {
					if ((i==0)||(j==0)||(k==0)||(l==0)) {
						data->stack[i][j][k][l]=0;
					}
					else {
						dr >> lineoftext;
						if (strcmp(lineoftext,".")){
							data->stack[i][j][k][l] =(int)floor(conversionfactor*(atof(lineoftext))+.5);
						}
						else data->stack[i][j][k][l] = infinity;
					}
					//cout <<"stack "<<i<<" "<<j<<" "<<k<<" "<<l<<"  "<<data->stack[i][j][k][l]<<"\n";
				}
				//cin >> m;
			}
		}
	}

	return 1;
}

/*      Function getdh      
Function gets the names of dhdata files to open */

void getdh(char *loop, char *stackf, char *tstackh, char *tstacki,char *tloop, char *miscloop, char *danglef, char *int22,char *int21,char *coax, char *tstackcoax,char *coaxstack, char *tstack, char *tstackm,char *triloop, char *int11,char *hexaloop,char *tstacki23, char *tstacki1n )
{
strcpy (loop,"../../NN_data/loop.dh");
strcpy (stackf,"../../NN_data/stack.dh");
strcpy (tstackh,"../../NN_data/tstackh.dh");
strcpy (tstacki,"../../NN_data/tstacki.dh");
strcpy (tloop,"../../NN_data/tloop.dh");
strcpy (miscloop,"../../NN_data/miscloop.dh");
strcpy (danglef,"../../NN_data/dangle.dh");
strcpy (int22,"../../NN_data/int22.dh");
strcpy (int21,"../../NN_data/int21.dh");
strcpy (coax,"../../NN_data/coaxial.dh");
strcpy (tstackcoax,"../../NN_data/tstackcoax.dh");
strcpy (coaxstack,"../../NN_data/coaxstack.dh");
strcpy (tstack,"../../NN_data/tstack.dh");
strcpy (tstackm,"../../NN_data/tstackm.dh");
strcpy (triloop,"../../NN_data/triloop.dh");
strcpy (int11,"../../NN_data/int11.dh");
strcpy (hexaloop,"../../NN_data/hexaloop.dh");
strcpy (tstacki23,"../../NN_data/tstacki23.dh");
strcpy (tstacki1n,"../../NN_data/tstacki1n.dh");
}

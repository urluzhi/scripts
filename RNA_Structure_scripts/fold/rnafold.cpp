#include <iostream>
#include "stdafx.h"
#include "TProgressDialog.cpp"
#include "arrayclass.cpp"
#include "dotarray.cpp"
#include "stackclass.cpp"
#include "stackstruct.cpp"
#include "forceclass.cpp"
#include "rna_library.cpp"
#include "structure.cpp"
#include "algorithm.cpp"

#define maxsequences 2000 //maximum number of sequences that can be folded
#define cntrl8 20
#define cntrl6 750
#define cntrl9 0
using namespace std;

void getdat(char *loop, char *stackf, char *tstackh, char *tstacki,
	char *tloop, char *miscloop, char *danglef, char *int22,
	char *int21,char *coax, char *tstackcoax,
    char *coaxstack, char *tstack, char *tstackm, char *triloop,
    char *int11, char *hexaloop, char *tstacki23, char *tstacki1n,
	char *datapath, bool isRNA);

int main(int argc, char *argv[]) {


	char loop[maxfil],stackf[maxfil],tstackh[maxfil],tstacki[maxfil],
		tloop[maxfil],miscloop[maxfil],danglef[maxfil],int22[maxfil],
      int21[maxfil],coax[maxfil],tstackcoax[maxfil],
      coaxstack[maxfil],tstack[maxfil],tstackm[maxfil],triloop[maxfil],int11[maxfil],hexaloop[maxfil],
	  tstacki23[maxfil], tstacki1n[maxfil];
   char list[maxfil],outfile[maxfil],label[maxfil],prefix[maxsequences][maxfil];
   char sequence[maxfil],ctfile[maxfil];
   char ctpath[maxfil],seqpath[maxfil],*datapath;
   
   register int c,i,j,ir;
   //int cntrl6,cntrl8,cntrl9,ir;
   datatable *data;
   structure *ct;
   //dynamic alloc datatalbes
   data=new datatable();
  //define the datapath from shell
	datapath=getenv("DATAPATH");
	if (datapath == NULL) {
		datapath="";
	}
	
	
   //Get required information:
   if (argc!=4) {
	cout << "rnafold.o rna.lis seq_dir/ ct_dir/\n";
	return 1;
	}
	else {
   	  strcpy(list,argv[1]);
   	  strcpy(seqpath,argv[2]);
   	  strcpy(ctpath,argv[3]);
   }



	ifstream inf;
	inf.open(list);


   //read the list of structures to be folded from a list file:
	for (ir=1;ir<=maxsequences;ir++) {
		
		inf >> prefix[ir];
		if(inf.eof()) break;
		cout << "ir = "<<ir<<" ct = "<<prefix[ir]<<"\n"<<flush;

}

	ir--;


	//open the data files -- must reside with executable
   getdat (loop,stackf,tstackh,tstacki,tloop,miscloop,danglef,
   	int22,int21,coax,tstackcoax,coaxstack,tstack,tstackm,triloop,int11,hexaloop,tstacki23, tstacki1n,datapath,true);
	if (opendat (loop,stackf,tstackh,tstacki,tloop,miscloop,danglef,int22,int21,
   	coax,tstackcoax,coaxstack,tstack,tstackm,triloop,int11,hexaloop,tstacki23, tstacki1n,data)==0) {
      	cout << "A data file was lost";
         cin >> i;
   }
  			    
   for (i=1;i<=ir;i++) {//fold and score every structure
	  ct = new structure;
      strcpy(sequence,seqpath);
      strcat(sequence,prefix[i]);
      strcat(sequence,".seq");
  	  cout << i << "  "<<sequence<<"\n"<<flush;
      strcpy(ctfile,ctpath);
      strcpy(ctfile,prefix[i]);
      strcat(ctfile,".ct");

  	  //open the current sequence
      openseq (ct,sequence);
  
	  //predict the secondary structure
	  dynamic (ct,data,cntrl6,cntrl8,cntrl9);

	  ctout(ct,ctfile);
	  delete ct;

	}
    delete data;
   return 0;

}



void errmsg(int err,int erri) {

if (err==30) {
	cout << "End Reached at traceback #"<<erri<<"\n";
   exit(1);
}
if (err==100) {
	cout << "error # "<<erri;
   exit(1);
}
switch (err) {
	case 1:
   	cout << "Could not allocate enough memory";
      break;
   case 2:
   	cout << "Too many possible base pairs";
      break;
   case 3:
   	cout << "Too many helixes in multibranch loop";
   case 4:
   	cout << "Too many structures in CT file";
   default:
   	cout << "Unknown error";
}
cin >> err;
exit(1);
return;

}

/*	Function getdat

	Function gets the names of data files to open

*/

void getdat(char *loop, char *stackf, char *tstackh, char *tstacki,
	char *tloop, char *miscloop, char *danglef, char *int22,
	char *int21,char *coax, char *tstackcoax,
    char *coaxstack, char *tstack, char *tstackm, char *triloop,
    char *int11, char *hexaloop, char *tstacki23, char *tstacki1n,
	char *datapath, bool isRNA)

{

 
  if( !isRNA) strcat( datapath,"dna");
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
 
  strcat (loop,"loop.dat");
  strcat (stackf,"stack.dat");
  strcat (tstackh,"tstackh.dat");
  strcat (tstacki,"tstacki.dat");
  strcat (tloop,"tloop.dat");
  strcat (miscloop,"miscloop.dat");
  strcat (danglef,"dangle.dat");
  strcat (int22,"int22.dat");
  strcat (int21,"int21.dat");
  strcat (triloop,"triloop.dat");
  strcat (coax,"coaxial.dat");
  strcat (tstackcoax,"tstackcoax.dat");
  strcat (coaxstack,"coaxstack.dat");
  strcat (tstack,"tstack.dat");
  strcat (tstackm,"tstackm.dat");
  strcat (int11,"int11.dat");
  strcat (hexaloop,"hexaloop.dat");
  strcat (tstacki23,"tstacki23.dat");
  strcat (tstacki1n,"tstacki1n.dat");
}


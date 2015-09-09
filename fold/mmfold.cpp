// fold a single RNA sequence.  Programmed by DHM, 8/17/05
//

#include "stdafx.h"
//#include "mseq.h"
#include "platform.h"
#include "rna_library.h"
#include "structure.h"
#include "algorithm.h"
#include <iostream>
#include <fstream>
#include <string>
using namespace std;

int openstring (structure *ct, string  *idtr, string *seqstr) ;

void getdat(char *loop, char *stackf, char *tstackh, char *tstacki,
	char *tloop, char *miscloop, char *danglef, char *int22,
	char *int21,char *coax, char *tstackcoax,
    char *coaxstack, char *tstack, char *tstackm, char *triloop,
    char *int11, char *hexaloop, char *tstacki23, char *tstacki1n,
	char *datapath, bool isRNA);

int main(int argc, char* argv[])
{
	char loop[maxfil],stackf[maxfil],tstackh[maxfil],tstacki[maxfil],
		tloop[maxfil],miscloop[maxfil],danglef[maxfil],int22[maxfil],
		int21[maxfil],coax[maxfil],tstackcoax[maxfil],
		coaxstack[maxfil],tstack[maxfil],tstackm[maxfil],triloop[maxfil],int11[maxfil],hexaloop[maxfil],
		tstacki23[maxfil], tstacki1n[maxfil],datapath[maxfil],*pointer;

	int maxtracebacks,percent,windowsize;
	structure *ct1,*ct2;
	datatable data;
	int ir,i;
	if(argc!=2) {
		cout<<"Command line only: ./mmfold seq.fa\n";	
		return 1;
	}
	pointer = getenv("DATAPATH");
	if (pointer!=NULL) strcpy(datapath,pointer);
	else strcpy(datapath,"");
	//open the data files -- must reside with executable
	//open the thermodynamic data tables
	getdat (loop, stackf, tstackh, tstacki,tloop, miscloop, danglef, int22,
          int21,coax, tstackcoax,coaxstack, tstack, tstackm, triloop,
          int11, hexaloop, tstacki23, tstacki1n, datapath, true);//the true indicates RNA parameters
	if (opendat (loop,stackf,tstackh,tstacki,tloop,miscloop,danglef,int22,int21,
   		coax,tstackcoax,coaxstack,tstack,tstackm,triloop,int11,hexaloop,tstacki23, tstacki1n,&data)==0) {
      	
		cout << "A data file was lost";
		cin >> i;
	}

	//open the sequences
	ifstream OpenFile(argv[1]);
	string id,seq1,seq2;
	while(!OpenFile.eof())	{
		getline(OpenFile,id);
 		if (id.compare(0,1,">") != 0) continue;
		getline(OpenFile,seq1);
		getline(OpenFile,seq2);
		ct1=new structure(1+1);
		ct2=new structure(1+1);
		openstring(ct1,&id,&seq1);
		openstring(ct2,&id,&seq2);
		//now do the structure prediction:
		dynamic(ct1,&data,0,0,0,0,1); 
		dynamic(ct2,&data,0,0,0,0,1);
	   	cout<<id<<"\t"<<ct1->energy[1]/10.0<<"\t"<<ct2->energy[1]/10.0<<"\n";
		//cout<<seq2<<endl;
		//cout <<ct2->nucs<<"seq2"<<endl;
		delete ct1;
		delete ct2;
	}
	OpenFile.close();

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
//open sequence string
int openstring (structure *ct, string  *idstr, string *seqstr) {
	int length,j,i;
	char base[2];
	ct->nnopair = 0;
	strcpy(ct->ctlabel[1],idstr->c_str());
	length = seqstr->length();
	ct->numofbases = length;
	ct->allocate(1+length);
	i=1;
	//ct->nucs[0]='N';
	//ct->nucs[length+1]='\0';
	for (j=0;j<length;j++) {//up to length-1 bcs of /n character
		seqstr->copy(base,1,j);
		base[1]='\0';
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
}

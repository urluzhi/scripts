//This file contains the global interface functions 

#include "stdafx.h"
#include <iostream>
using namespace std;

void errmsg(int err,int erri) 
{
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

/////////////////////////////////////////////////////////////////////////////
// CRNAstructureApp message handlers




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

/*

strcat (loop,"\\");
strcat (stackf,"\\");
strcat (tstackh,"\\");
strcat (tstacki,"\\");
strcat (tloop,"\\");
strcat (miscloop,"\\");
strcat (danglef,"\\");
strcat (int22,"\\");
strcat (int21,"\\");
strcat (triloop,"\\");
strcat (coax,"\\");
strcat (tstackcoax,"\\");
strcat (coaxstack,"\\");
strcat (tstack,"\\");
strcat (tstackm,"\\");
strcat (int11,"\\");
strcat (hexaloop,"\\");
strcat (tstacki23,"\\");
strcat (tstacki1n,"\\");
*/


if (!isRNA) {
	//load DNA parameters


	strcat (loop,"dna");
	strcat (stackf,"dna");
	strcat (tstackh,"dna");
	strcat (tstacki,"dna");
	strcat (tloop,"dna");
	strcat (miscloop,"dna");
	strcat (danglef,"dna");
	strcat (int22,"dna");
	strcat (int21,"dna");
	strcat (triloop,"dna");
	strcat (coax,"dna");
	strcat (tstackcoax,"dna");
	strcat (coaxstack,"dna");
	strcat (tstack,"dna");
	strcat (tstackm,"dna");
	strcat (int11,"dna");
	strcat (hexaloop,"dna");
	strcat (tstacki23,"dna");
	strcat (tstacki1n,"dna");

}

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

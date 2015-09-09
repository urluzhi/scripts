#if !defined(INTERMOLECULAR_H)
#define INTERMOLECULAR_H


/*=======================================================================
intermolecular.h and intermolecular.cpp inculde funcions calculating and report
different free energy for binding in OligoWalk.
They are revised based on Mathews' code from RNAStructure.
olig() generate the energy data;
report save the generated data;
siprefileter and filterbysirna() were used to prefilter and postfileter the 
functional siRNA respectively, using criterias other than thermodynamics


															----Nov., 2005
															John(Zhi Lu)

=======================================================================*/

#include <math.h>
#include "stdafx.h"
#include "rna_library.h"
#include "thermo.h"
#include "pclass.h"
#include "algorithm.h"
#include "siRNAfilter.h"
#include "alltrace.h"
#include "stochastic.h"


#define FILTER_PASS 6	//The criteria for prefiltering of functional siRNA in OligoWalk
#define maxfile 75		//maximum length of file names

//=======================================================================
//this structure rddata contains all the thermodynamic parameters needed
//for determining the stability of an RNA-DNA duplex
struct rddata {

   short int stack[5][5][5][5];
   short int init;
};


//=======================================================================
//olig is the backend function for oligo walk
void olig(int isdna, int option, structure *ct, int length,double c, int **table,int **numofsubstructures,datatable& data,datatable& ddata, rddata *hybriddata, int Usesub,TProgressDialog *update,
		  thermo* helixstack,int start, int stop, siPREFILTER *prefilter,int foldsize,int distance,char *shapefile,int *TEST,bool WRITE=false); 

//=======================================================================
//this function reads data into structure rddata
int readrd (rddata* data,char* dnarna);

//=======================================================================
//this function writes a tab delimited file with the oligowalk data
void report(char* filename, structure *ct, int **table,int **numofsubstructures, int length, int isdna,
			double conc, int Usesub,int start,int stop,siPREFILTER *prefilter,int scansize,
			bool *mask=NULL, double asuf=0, double tofe=0, double fnnfe=0);

//=======================================================================
//returns the base complementary to base i
//with 1 == A, 2 == C, 3 == G, 4 == U
int complement(int i, structure *ct);

//=======================================================================
//post-filter of functional siRNA written by Dave.
void filterbysirna(structure *ct, int **table, int length, datatable *data, 
				   bool *mask, double asuf, double tofe, double fnnfe);

//=======================================================================
//copy the arrays to be reused with different index when the folded region move one nucleotide to the right
inline void scancopy(OligoPclass *region, OligoPclass *copyregion);
//copy every arrays to be reused when the folded region did not move right yet 
inline void scancopyend(OligoPclass *region, OligoPclass *copyregion) ;

#endif


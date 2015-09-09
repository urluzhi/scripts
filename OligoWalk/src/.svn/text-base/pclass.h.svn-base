#if !defined(PCLASS_H)
#define PCLASS_H


/*======================================================================
	Two files (pclass.h and pclass.cpp) are made based on pfunction.h,
pfunction.cpp from the old version of RNAStructure.

	A new class ,Pclass, would be used instead of pfunction(...).
An inherited class from Pclass, OligoPclass, was also added to reused the 
partition function array (v,w,w5,...) for OligoWalk program.

														----Nov., 2005
															John(Zhi Lu)
  =======================================================================*/
//#pragma once 

#include <math.h>
#include "rna_library.h"
#include "algorithm.h"


//=======================================================================
/*
pfunctionclass encapsulates the large 2-d arrays of w and v, used by the 
partition function
*/
class pfunctionclass {
	
	int Size;

public:
	int k;
	PFPRECISION **dg;
	PFPRECISION infinite;
	
	pfunctionclass(int size);//the constructor allocates the space needed by the arrays
	~pfunctionclass();//the destructor deallocates the space used

     
   	inline PFPRECISION &f(int i, int j);
						//f is an integer function that references the correct
						//element of the array
};

//f is an integer function that references the correct element of the array
inline PFPRECISION &pfunctionclass::f(int i, int j) {

   


   if (i>j) {
        return infinite;
    }
   else if (i>Size) return f(i-Size,j-Size);
   else return dg[i][j-i];
         
}


//=======================================================================
//this class contains all the info read from thermodynamic data files
class pfdatatable {
 
public:
	int numoftriloops,numoftloops,numofhexaloops;
	int itloop[maxtloop+1],itriloop[maxtloop+1],ihexaloop[maxtloop+1];
	int maxintloopsize;

	PFPRECISION poppen [5],maxpen,eparam[11],dangle[6][6][6][3],inter[31],bulge[31],
				hairpin[31],stack[6][6][6][6],tstkh[6][6][6][6],tstki[6][6][6][6],
				tloop[maxtloop+1],iloop22[6][6][6][6][6][6][6][6],
				iloop21[6][6][6][6][6][6][6],iloop11[6][6][6][6][6][6],
				coax[6][6][6][6],tstackcoax[6][6][6][6],coaxstack[6][6][6][6],
				tstack[6][6][6][6],tstkm[6][6][6][6],auend,gubonus,cint,cslope,c3,
				efn2a,efn2b,efn2c,triloop[maxtloop+1],init,mlasym,strain,
				singlecbulge,tstki23[6][6][6][6],tstki1n[6][6][6][6];
	PFPRECISION hexaloop[maxtloop+1];
	PFPRECISION scaling;
	PFPRECISION prelog;
	PFPRECISION temp;
	
	pfdatatable(datatable *indata, PFPRECISION Scaling, PFPRECISION Temp=310.15);
	pfdatatable();
	void rescaledatatable(PFPRECISION rescalefactor); //rescale the entries in datatable
};



//=======================================================================
//some functions used for calculation and storage of partition function

PFPRECISION pfchecknp(bool lfce1,bool lfce2);
PFPRECISION erg1(int i,int j,int ip,int jp,structure *ct,pfdatatable *data);
		//calculates quilibrium constant of stacked base pairs
PFPRECISION erg2(int i,int j,int ip,int jp,structure *ct,pfdatatable *data,char a,
	char b);
		//calculates quilibrium constant of a bulge/internal loop
PFPRECISION erg3(int i,int j,structure *ct,pfdatatable *data,char dbl);
		//calculates quilibrium constant of a hairpin loop
PFPRECISION erg4(int i,int j,int ip,int jp,structure *ct,pfdatatable *data,
	bool lfce);
		//calculates quilibrium constant of a dangling base
PFPRECISION erg2in(int i,int j,int ip,int jp,structure *ct, pfdatatable *data);
		//calculate the quilibrium constant of the intorior part of a internal loop
		//where i is paired to j; ip is paired to jp; ip > i; j > jp
PFPRECISION erg2ex(int i,int j,int size, structure *ct, pfdatatable *data);
		//calculate the quilibrium constant of the exterior part of internal loop
		//where i is paired to j; ip is paired to jp; ip > i; j > jp
PFPRECISION penalty(int i,int j,structure *ct,pfdatatable *data);
		//calculates end of helix penalty
PFPRECISION penalty2(int i, int j, pfdatatable *data);
PFPRECISION ergcoaxflushbases(int i, int j, int ip, int jp, structure *ct, pfdatatable *data);
PFPRECISION ergcoaxinterbases1(int i, int j, int ip, int jp, structure *ct, pfdatatable *data);
PFPRECISION ergcoaxinterbases2(int i, int j, int ip, int jp, structure *ct, pfdatatable *data);


//functions for the storage of arrays of partition function

void write(ofstream *out,double *i);
void read(ifstream *out,double *i) ;

void readpfsave(char *filename, structure *ct,PFPRECISION *w5, PFPRECISION *w3, 
				pfunctionclass *v, pfunctionclass *w, pfunctionclass *wmb, 
				pfunctionclass *wl, pfunctionclass *wmbl, pfunctionclass *wcoax,
				forceclass *fce, PFPRECISION *scaling, bool *mod, bool *lfce, 
				pfdatatable *data);
PFPRECISION calculateprobability(int i, int j, pfunctionclass *v, 
								 PFPRECISION *w5, structure *ct, pfdatatable *data,
								 bool *lfce, bool *mod, PFPRECISION scaling,
								 forceclass *fce);
inline void rescale(int i, int j,structure *ct, pfdatatable *data, pfunctionclass *v, pfunctionclass *w, pfunctionclass *wl, pfunctionclass *wcoax,
			 pfunctionclass *wmb,pfunctionclass *wmbl, PFPRECISION *w5, PFPRECISION *w3, PFPRECISION **wca, PFPRECISION **curE, PFPRECISION **prevE, PFPRECISION rescalefactor); 
		//function to rescale all arrays when partition function calculation is headed out of bounds

int pfshape(structure *ct, PFPRECISION  temp);//scale the shape energy into partitionfunction

//=======================================================================
/*
Pclass would be used instead of pfunction(). 
This will fill the array (v,w,...) of partition function calculation with
O(N3) in calculation time, with limited or unlimited internal loop sizes.
The old algorithm of O(N4) is also avaiable by calling oldfill,oldpartition()
*/
class Pclass {

protected:	
	int i,j,h,d,maxj,lowi;
	int dp,ip,jp,ii,jj,jpf,jf,bl,ll;
	int k,l,m,n,o,p;
	int before,after;
	int inc[6][6];
	
	bool *lfce,*mod;//[maxbases+1][maxbases+1];
	PFPRECISION e;
	PFPRECISION twoscaling,rarray;
	
	forceclass *fce;	
	
	

public:
	int number;	
	pfdatatable* data;
	structure* ct;
	PFPRECISION *w5,*w3,**wca,**curE,**prevE,**tempE;
	pfunctionclass *w,*v,*wmb,*wl,*wmbl,*wcoax;
	
	Pclass(structure *CT,pfdatatable *DATA);
	~Pclass();
	
	inline void limitdist();//This function handles the case where base pairs are not
					//not allowed to form between nucs more distant
					//than ct->maxdistance
	void store(char *save);//This function will save the arrays in a binary file

	//-------------------------------------------------------------------------------------------------
	//old fill routine with O(N4) for unlimited internal loops
	void fillw3();//This is a inline function to calc w3 in the old fill routine (O(N4) in time)
	void oldfill();//This is a old fill routine,O(N4) for unlimited internall loop
	void oldpartition(bool quickQ,PFPRECISION *Q=NULL,TProgressDialog* update=NULL,char *save=NULL);
					//This is partition function using oldfill()
	
	//--------------------------------------------------------------------------------------------------
	//improved fill routine with O(N3) for unlimited internal loops
	void fillw5();//This is a inline function to calc w5 in the improved fill routine (O(N3) in time)
	void interprefill();//This function prefill curE and prevE for internal loops's energy for the improved algorithm
	void fill();//This function fill the array for a partition function of the sequence
	void partition( bool quickQ,PFPRECISION *Q=NULL,TProgressDialog* update=NULL,char *save=NULL);
					//If quickQ == true, return the partition function value in Q
					//(Q is scaled in case overflowed, not a real value)
					//otherwise, save the partial partition functions
					//in the datafile named save
					//If updates on progress are unwanted, set update=NULL
};



//=======================================================================
/*
OligoPclass is inherited from Pclass
with additional functions for reusing the partition function array in the
OligoWalk refolding.

*/
class OligoPclass:public Pclass {

public:
	PFPRECISION *copyw5,**copywca;
	pfunctionclass *copyw,*copyv,*copywmb,*copywl,*copywmbl,*copywcoax;
			//These arrays and classes store the informations to be copied in refolding
	OligoPclass(structure *CT, pfdatatable *DATA);//initiate the copy arrays
	~OligoPclass();

	//----------------------------------------------------------------------------------------------
	//folding the whole length of target sequence	
	void partition4refill(PFPRECISION *Q,char *save=NULL);
			//normal partition function plus a copy of the arrays of v,w,w5 ...
	void refill(structure *CT,PFPRECISION *Qc, int start, int stop,PFPRECISION &rescaleinrefill,char *save=NULL);
			//recalculate the partionfunction for constrained sequence
			//assisted with copying the stored information in copyw,copyv ...	
			
	//----------------------------------------------------------------------------------------------
	//refolding at different region site on target
	void scanfill(structure *CT,PFPRECISION *Q,int reverse=0,char *save=NULL);//refold different region on target
	void scanconstrain(structure *CT,PFPRECISION *Qc,int start, int stop,PFPRECISION &rescaleinrefill,char *save=NULL);
																//refold the region on target with constrain
	
	void reset4oligo(structure *CT);//reset the arrays to be refilled again for other sequence
};




#endif


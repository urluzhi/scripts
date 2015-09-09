
#if !defined (RNA_LIBRARY)
#define RNA_LIBRARY


#include "structure.h"
#include <cmath>





struct datatable //this structure contains all the info read from thermodynamic
							//data files
{
 short int poppen [5],maxpen,eparam[11],dangle[6][6][6][3],inter[31],bulge[31],
		hairpin[31],stack[6][6][6][6],tstkh[6][6][6][6],tstki[6][6][6][6],
		tloop[maxtloop+1][2],numoftloops,iloop22[6][6][6][6][6][6][6][6],
		iloop21[6][6][6][6][6][6][6],iloop11[6][6][6][6][6][6],
      coax[6][6][6][6],tstackcoax[6][6][6][6],coaxstack[6][6][6][6],
      tstack[6][6][6][6],tstkm[6][6][6][6],auend,gubonus,cint,cslope,c3,
      efn2a,efn2b,efn2c,triloop[maxtloop+1][2],numoftriloops,init,mlasym,strain,numofhexaloops,
	  singlecbulge,tstki23[6][6][6][6],tstki1n[6][6][6][6];
 int hexaloop[maxtloop+1][2];
float RT;                          //temperature*gas constant (used by multi-state bulge)
float prelog;
datatable();


};

struct stackstruct //this structure contains a stack of data, used by
						//	functions that analyze a structure piecewise
{
	int stk[STACKSIZE][4],sp;
};


class stackclass
{
	



	void allocate_stack();

public:
	short size,**stack,maximum;
	integersize *stackenergy;

	stackclass(short int stacksize = 50);
	
	


	bool pull(short int *i,short int *j, short int *open, integersize *energy, short int *pair);
	void push(short int i,short int j, short int open, integersize energy, short int pair);
	
	void delete_array();


	~stackclass();

};


int ergcoaxflushbases(int i, int j, int ip, int jp, datatable *data);
//this function calculates flush coaxial stacking
//it requires sequence in i,j,ip, and jp
int ergcoaxinterbases1(int i, int j, int ip, int jp, int k, int l, datatable *data); 
//this funtion calculates an intervening mismatch coaxial stack
int ergcoaxinterbases2(int i, int j, int ip, int jp, int k, int l, datatable *data); 
//this funtion calculates an intervening mismatch coaxial stack
int ergcoaxflushbases(int i, int j, int ip, int jp, structure *ct, datatable *data);
int ergcoaxinterbases1(int i, int j, int ip, int jp, structure *ct, datatable *data);
int ergcoaxinterbases2(int i, int j, int ip, int jp, structure *ct, datatable *data);
int decon1(int x);//used by ergmulti to find a nucleotide from a base pair
int decon2(int x);//used by ergmulti to find a nucleotide from a base pair
int ergmulti(int st, int ip, structure *ct, datatable *data, bool simplemb);
//calculate the multi branch loop free energy for a loop starting at nuc ip
//	in structure number st of ct
int ergexterior(int st, structure *ct, datatable *data);
//calculate the exterior loop free energy in structure number ip



//write is used to write data to a save file
void write(ofstream *out,short *i);
void write(ofstream *out,bool *i);
void write(ofstream *out,int *i);
void write(ofstream *out,char *i);
void write(ofstream *out,float *i);
void writesinglechar(ofstream *out,char *i); 

//read is used to read data from a save file
void read(ifstream *out,short *i);
void read(ifstream *out,bool *i);
void read(ifstream *out,int *i);
void read(ifstream *out,char *i);
void read(ifstream *out,float *i);
void readsinglechar(ifstream *out,char *i);

int opendat (char *loop2,char *stackf,char *tstackh,char *tstacki,
		char *tloop,char *miscloop, char *danglef, char *int22, char *int21,
      char *coax,char *tstackcoax,char *coaxstack,
      char *tstack,char *tstackm, char *triloop, char *int11, char *hexaloop,
	  char *tstacki23, char *tstacki1n, datatable* data);
		//gets thermodynamic data from data files

void push(stackstruct *stack,int a,int b,int c,int d);//push info onto the stack
void pull(stackstruct *stack,int *i,int *j,int *open,int *null,int *stz);
														//pull info from the stack

int erg1(int i,int j,int ip,int jp,structure *ct,datatable *data);
		//calculates energy of stacked base pairs
int erg2(int i,int j,int ip,int jp,structure *ct,datatable *data,char a,
	char b);
		//calculates energy of a bulge/internal loop
int erg3(int i,int j,structure *ct,datatable *data,char dbl);
		//calculates energy of a hairpin loop
int erg4(int i,int j,int ip,int jp,structure *ct,datatable *data,
	bool lfce);
		//calculates energy of a dangling base
inline int SHAPEend(int i, structure *ct);//calculate the SHAPE pseudo energy for a single
									//paired nucleotide

//this function calculates whether a terminal pair i,j requires the end penalty
inline int penalty(int i,int j,structure* ct, datatable *data) {
	
   if (ct->numseq[i]==4||ct->numseq[j]==4)
   	return data->auend;
	else return 0;//no end penalty
	


}
	
int penalty2(int i, int j, datatable *data);
int ergcoax(int i, int j, int ip, int jp, int k, structure *ct, datatable *data);
	//returns the free energy of coaxial stacking of i-j onto ip-jp






#endif


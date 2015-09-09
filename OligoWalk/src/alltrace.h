#include <math.h>
#include "stdafx.h"
#include "rna_library.h"
#include "pclass.h"
#include "platform.h"



#define MaxStructure 10000000 // the max number of structures explored
////////////////////////////////////////////////////////////////////////
//arrayclass encapsulates the large 2-d arrays of w and v, used by the dynamic
//	algorithm
class atarrayclass {
   int Size;

   public:
   	
      int k;
      integersize **dg;
      integersize infinite;
      

      //the constructor allocates the space needed by the arrays
   		atarrayclass(int size);

      //the destructor deallocates the space used
      ~atarrayclass();

      //f is an integer function that references the correct element of the array
   	inline integersize &f(int i, int j);
};

void alltrace(bool writect,int *ctenergy,structure* ct,datatable* data, short percentdelta, short absolutedelta, TProgressDialog* update, char* save);
void readalltrace(char *filename, structure *ct, 
			 short *w5,  
			 atarrayclass *v, atarrayclass *w, atarrayclass *wmb, atarrayclass *wmbl, atarrayclass *wl, atarrayclass *wcoax,
			 atarrayclass *w2, atarrayclass *wmb2, forceclass *fce, bool *lfce, bool *mod, datatable *data);

void realltrace(bool writect,int *ctenergy,char *savefilename, structure *ct, short percentdelta, short absolutedelta, char *ctname = NULL);


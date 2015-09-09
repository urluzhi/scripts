/*
 * Dynalign
 *  by David Mathews, copyright 2002, 2003,2004,2005
 *
 * Contributors: Chris Connett and Andrew Yohn, 2006
 */

#include "arrayclass.h"

#include "defines.h"

arrayclass::arrayclass(int size) {
	

	infinite = INFINITE_ENERGY;

  Size = size;
  register int i,j;
  dg = new integersize *[size+1];
    
	for (i=0;i<=(size);i++)  {
    dg[i] = new integersize [size+1];
  }
  for (i=0;i<=size;i++) {
    for (j=0;j<size+1;j++) {
      dg[i][j] = INFINITE_ENERGY;
    }
  }
}

// the destructor deallocates the space used
arrayclass::~arrayclass() {

      	
	int i;
       	
    for (i=0;i<=Size;i++) {
        delete[] dg[i];
    }
     delete[] dg;
}

/*
 * Dynalign
 *  by David Mathews, copyright 2002, 2003,2004,2005
 *
 * Contributors: Chris Connett and Andrew Yohn, 2006
 */

#ifndef ARRAYCLASS_H
#define ARRAYCLASS_H

#include "defines.h"

// arrayclass encapsulates the large 2-d arrays of w and v, used by
// the dynamic algorithm

class arrayclass {
private:
  int Size;

public:
  int k;
  integersize **dg;
  integersize infinite;

  // the constructor allocates the space needed by the arrays
  arrayclass(int size);
  
  // the destructor deallocates the space used
  ~arrayclass();

  // f is an integer function that references the correct element of
  // the array
  integersize &f(int i, int j);
};

inline integersize &arrayclass::f(int i, int j) {
   if (i > Size) {
     i -= Size;
     j -= Size;
   }

   if (i > j) {
        return infinite;
   }
   
   return dg[i][j-i];
}

#endif

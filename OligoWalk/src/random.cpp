#include "stdafx.h"

#define IM1 2147483563
#define IM2 2147483399
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define PRECISION 1.0e-300
#define WUP 12

#include "random.h"
 

//Random number generator of Numerical recipes in C, "ran2," period of about 2e18

randomnumber::randomnumber() {
	IM1INV = 1.0e0 / IM1;
	IMM1 = IM1-1;
	DIVISOR = (1+IMM1/TABSIZE);
	//RNMX=1.0e0 - PRECISION;

	seed(1234);  //Seed with a basic number in case the user does not seed

}

void randomnumber::seed(long seeddouble) {
	int i;
	long j;
	
	idum = seeddouble;
	

	if (idum<1) idum = 1; //prevent seeding with zero
	idum2 = idum;
	for (i=TABSIZE+WUP;i>=0;i--) {
		j = idum/IQ1;
		idum = IA1*(idum-j*IQ1)-j*IR1;
		if (idum < 0) idum = idum+IM1;
		if (i<TABSIZE) iv[i]=idum;

	}
	iy = iv[0];



}

double randomnumber::roll() {
	long j;
	int i;

	j = idum/IQ1;
	idum=IA1*(idum-j*IQ1)-j*IR1;
	if (idum<0) idum = idum + IM1;
	j = idum2/IQ2;
	idum2=IA2*(idum2-j*IQ2)-j*IR2;
	if (idum2<0) idum2 = idum2 + IM2;
	i = (int) (iy/DIVISOR);
	iy=iv[i]-idum2;
	iv[i] = idum;
	if (iy<1) iy=iy+IMM1;
	return IM1INV*iy;

}


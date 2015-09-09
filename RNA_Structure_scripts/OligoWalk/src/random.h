#if !defined(RANDOM_H)
#define RANDOM_H
#define TABSIZE 32

class randomnumber {

	private:
		long IMM1,DIVISOR,/*RNMX,*/idum2,idum,iy,iv[TABSIZE];
		double IM1INV;

	public:
		randomnumber();
		void seed(long seeddouble);
		double roll();



};


#endif


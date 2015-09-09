
#include "stdafx.h"
#include "platform.h"
#include "rna_library.h"
#include <cmath>

using namespace std;


//***********************************code for Structures:


datatable::datatable()
{
int a,b,c,d,e,f,g,h;


for (a=0;a<=5;a++) {
	for (b=0;b<=5;b++) {
		for (c=0;c<=5;c++) {
			for (d=0;d<=5;d++) {
				for (e=0;e<=5;e++) {
					for (f=0;f<=5;f++) {
               	iloop11[a][b][c][d][e][f] = infinity;
						for (g=0;g<=5;g++) {
			iloop21[a][b][c][d][e][f][g] = infinity;
         for (h=0;h<=5;h++) {
          	iloop22[a][b][c][d][e][f][g][h] = infinity;
         }
						}
					}
				}
			}
		}
	}
}

//Set RT equal to the defined rt (in defines.h).
//This assumes folding at 37 degrees C.
RT=(float) rt;


};

/*  Function opendat

		Function opens data files to read thermodynamic data
*/

int opendat (char *loop2,char *stackf,char *tstackh,char *tstacki,
		char *tloop,char *miscloop, char *danglef, char *int22, char *int21,
      char *coax,char *tstackcoax,char *coaxstack,
      char *tstack,char *tstackm, char *triloop, char *int11, char *hexaloop, 
	  char *tstacki23, char *tstacki1n,datatable* data)
{
char lineoftext[100],base[110];
int count,i,j,k,l, m, a, b, c, d,e,f,g;
float temp;
FILE *check;

//eparam[1] is a basepair bonus
//eparam[2] is a bulge loop bonus
//eparam[3] is an interior loop bonus

 ifstream ml1;
 ifstream lo1;
 ifstream st1;
 ifstream th1;
 ifstream ti1;
 ifstream tl1;
 ifstream da1;
 ifstream in1;
 ifstream in2;
 ifstream tri;
 ifstream co1;
 ifstream co2;
 ifstream co3;
 ifstream st2;
 ifstream tsm;
 ifstream i11;
 ifstream hl1;
 ifstream t23;
 ifstream t1n;

 //ofstream out("check.out");


//check that all the files exist with a C i/o function
if ((check = fopen(miscloop, "r"))
	== NULL) {
	return 0;
}

fclose(check);

if ((check = fopen(loop2, "r"))
	== NULL) {
	return 0;
}

fclose(check);

if ((check = fopen(stackf, "r"))
	== NULL) {
	return 0;
}

fclose(check);

if ((check = fopen(tstackh, "r"))
	== NULL) {
	return 0;
}

fclose(check);

if ((check = fopen(tstacki, "r"))
	== NULL) {
	return 0;
}

fclose(check);

if ((check = fopen(tloop, "r"))
	== NULL) {
	return 0;
}

fclose(check);

if ((check = fopen(danglef, "r"))
	== NULL) {
	return 0;
}

fclose(check);

if ((check = fopen(int22, "r"))
	== NULL) {
	return 0;
}

fclose(check);

if ((check = fopen(int21, "r"))
	== NULL) {
	return 0;
}

fclose(check);

if ((check = fopen(triloop, "r"))
	== NULL) {
	return 0;
}

fclose(check);

if ((check = fopen(coax, "r"))
	== NULL) {
	return 0;
}

fclose(check);

if ((check = fopen(tstackcoax, "r"))
	== NULL) {
	return 0;
}

fclose(check);

if ((check = fopen(coaxstack, "r"))
	== NULL) {
	return 0;
}

fclose(check);

if ((check = fopen(tstack, "r"))
	== NULL) {
	return 0;
}

fclose(check);

if ((check = fopen(tstackm, "r"))
	== NULL) {
	return 0;
}

fclose(check);

if ((check = fopen(hexaloop, "r"))
	== NULL) {
	return 0;
}

fclose(check);

if ((check = fopen(tstacki23, "r"))
	== NULL) {
	return 0;
}

fclose(check);


if ((check = fopen(tstacki1n, "r"))
	== NULL) {
	return 0;
}

fclose(check);


if ((check = fopen(int11,"r")) == NULL) return 0;

fclose(check);


/* Read information from miscloop */

// the key sequence "-->" now indicates a record


ml1.open(miscloop);

ml1>>lineoftext;
 while(strcmp(lineoftext,"-->")) ml1>>lineoftext;


 ml1 >> (data->prelog);

 data->prelog = (data->prelog)*conversionfactor;


 ml1>>lineoftext;
 while(strcmp(lineoftext,"-->")) ml1>>lineoftext;

 ml1 >> temp;
 data->maxpen = (int) (temp*conversionfactor + .5);


 ml1>>lineoftext;
 while(strcmp(lineoftext,"-->")) ml1>>lineoftext;

 for (count=1;count<= 4;count ++)
 { ml1 >> temp;
 (data->poppen[count])= (int) (temp*conversionfactor + .5);



 } 										//this reads float values, converts
											// 	them int and assigns them into
											//		array poppen
 ml1>>lineoftext;
 while(strcmp(lineoftext,"-->")) ml1>>lineoftext;
 

 data->eparam[1] = 0;						 // assign some variables that are
 data->eparam[2] = 0;						 //	"hard-wired" into code
 data->eparam[3] = 0;
 data->eparam[4] = 0;
 ml1 >> temp;
 data->eparam[5] = (short) (floor (temp*conversionfactor+.5));  //constant multi-loop penalty



 ml1 >> temp;
 data->eparam[6] = (short) (floor (temp*conversionfactor+.5));



 data->eparam[7] = 30;
 data->eparam[8] = 30;
 data->eparam[9] = -500;
 ml1 >> temp;
 data->eparam[10] = (short) (floor (temp*conversionfactor+.5));
 ml1>>lineoftext;
 while(strcmp(lineoftext,"-->")) ml1>>lineoftext;

 ml1 >> temp;

 if (ml1.peek()==EOF) {
  	//these are old energy rules -- treat the other constants properly
   data->efn2a = data->eparam[5];
   data->efn2b = data->eparam[6];
   data->efn2c = data->eparam[10];
   data->strain = 0;
   data->mlasym=0;
   data->auend=0;
   data->gubonus=0;
   data->cslope = 0;
   data->cint=0;
   data->c3=0;
   data->init=0;
   
   data->singlecbulge = 0;


 }

 else {



 	data->efn2a = (short) (floor (temp*conversionfactor+.5));  //constant multi-loop penalty for efn2

	 ml1 >> temp;
 	data->efn2b= (short) (floor(temp*conversionfactor+.5));

 	ml1 >> temp;
 	data->efn2c= (short) (floor(temp*conversionfactor+.5));

	ml1>>lineoftext;
	while(strcmp(lineoftext,"-->")) ml1>>lineoftext;
	ml1>>temp;
	data->mlasym = (short) (floor (temp*conversionfactor+.5));
	
	
	ml1>>lineoftext;
	while(strcmp(lineoftext,"-->")) ml1>>lineoftext;
	ml1>>temp;
	data->strain = (short) (floor (temp*conversionfactor+.5));


 	//now read the terminal AU penalty:
   ml1>>lineoftext;
 	while(strcmp(lineoftext,"-->")) ml1>>lineoftext;
   ml1>> temp;

 	data->auend = (short) (floor (temp*conversionfactor+.5));

 	//now read the GGG hairpin bonus:
   ml1>>lineoftext;
 	while(strcmp(lineoftext,"-->")) ml1>>lineoftext;

 	ml1 >> temp;
 	data->gubonus = (short) (floor (temp*conversionfactor+.5));

	//now read the poly c hairpin penalty slope:
   ml1>>lineoftext;
 	while(strcmp(lineoftext,"-->")) ml1>>lineoftext;

 	ml1 >> temp;
 	data->cslope = (short) (floor (temp*conversionfactor+.5));

 	//now read the poly c hairpin penalty intercept:
   ml1>>lineoftext;
 	while(strcmp(lineoftext,"-->")) ml1>>lineoftext;

 	ml1 >> temp;
 	data->cint = (short) (floor (temp*conversionfactor+.5));

 	//now read the poly c penalty for a loop of 3:
   ml1>>lineoftext;
 	while(strcmp(lineoftext,"-->"))ml1>>lineoftext;

 	ml1 >> temp;
 	data->c3 = (short) (floor (temp*conversionfactor+.5));
   ml1>>lineoftext;
 	while(strcmp(lineoftext,"-->")) ml1>>lineoftext;
 	
   ml1 >> temp;
 	data->init = (short) (floor (temp*conversionfactor+.5));

 	//now read the GAIL rule indicator
	//ml1>>lineoftext;
 	//while(strcmp(lineoftext,"-->")) ml1>>lineoftext;

 	//ml1>> temp;
 	//data->gail = (short) (floor (temp+.5));

	//now read the single C bulge bonus
	ml1>>lineoftext;
 	while(strcmp(lineoftext,"-->")) ml1>>lineoftext;

 	ml1>> temp;
 	data->singlecbulge = (short) (floor (temp*conversionfactor+.5));


 }

 ml1.close();

 /*	read info from dangle */
 //add to dangle the case where X (represented as 0) is looked up
da1.open(danglef);
for (l = 1;l <=2; l++){
	for (i = 0;i <=5; i++){
		if ((i!=0)&&(i!=5)) for (count=1;count <=60;count++) da1 >> lineoftext;
		for (j=0;j<=5; j++) {
			for (k=0;k<=5; k++) {
				if ((i==0)||(j==0)||(k==0)) {
				    data->dangle[i][j][k][l] = 0;
				}
            else if ((i==5)||(j==5)||(k==5)) {
             	data->dangle[i][j][k][l] = 0;
            }
				else {
				    da1 >> lineoftext;
                //cout << lineoftext<<"\n";
				    if (strcmp(lineoftext,".")){
					data->dangle[i][j][k][l] = (short) (floor (conversionfactor*(atof(lineoftext))+.5));
				    }
				    else data->dangle[i][j][k][l] = infinity;
				}
            //cout <<"dangle = "<<i<<" "<<j<<" "<<k<<" "<<l<<" "<<data->dangle[i][j][k][l]<<"\n";
			}
         //cin >> m;
		}
	}
}

 da1.close();

/*	read info from loop for internal loops, hairpin loops, and bulge loops */

lo1.open(loop2);
for (count = 1; count <=26; count++) lo1 >> lineoftext; //get past text in file
for (i=1;i <= 30; i++) {
	lo1 >> lineoftext;//get past the size column in table
	lo1 >> lineoftext;
					if (strcmp(lineoftext,".")){
					data->inter[i] = (short) (floor (conversionfactor*(atof(lineoftext))+.5));
					}
					else data->inter[i] = infinity;

               //cout <<"inter = "<<data->inter[i]<<"\n";
	lo1 >> lineoftext;
					if (strcmp(lineoftext,"."))
					data->bulge[i] = (short) (floor(conversionfactor*(atof(lineoftext))+.5));
					else data->bulge[i] = infinity;

               //cout <<"bulge = "<<data->bulge[i]<<"\n";
	lo1 >> lineoftext;
					if (strcmp(lineoftext,".")){
					data->hairpin[i] = (short) (floor(conversionfactor*(atof(lineoftext))+.5));
					}
					else data->hairpin[i] = infinity;

               //cout <<"hair = "<<data->hairpin[i]<<"\n";
}

 lo1.close();

/* Read info from stack */
//add to the stack table the case where X (represented as 0) is looked up:

st1.open(stackf);
for (count=1;count<=42;count++) st1 >> lineoftext;//get past text in file
for (i=0;i<=5;i++) {
	if ((i!=0)&&(i!=5)) for (count=1;count<=60;count++) st1 >> lineoftext;
	for (k=0;k<=5;k++) {
		for (j=0;j<=5;j++) {
			for (l=0;l<=5;l++) {
				if ((i==0)||(j==0)||(k==0)||(l==0)) {
					data->stack[i][j][k][l]=0;
				}
            else if ((i==5)||(j==5)||(k==5)||(l==5)) {
             	data->stack[i][j][k][l] = infinity;
            }
				else {
					st1 >> lineoftext;
					if (strcmp(lineoftext,".")){
						data->stack[i][j][k][l] =(short) (floor(conversionfactor*(atof(lineoftext))+.5));
					}
					else data->stack[i][j][k][l] = infinity;
				}

			}

		}
	}
}

 st1.close();




/* Read info from tstackh */
//add to the tstackh table the case where X (represented as 0) is looked up:

th1.open(tstackh);
for (count=1;count<=46;count++) th1 >> lineoftext;//get past text in file
for (i=0;i<=5;i++) {
	if ((i!=0)&&(i!=5)) for (count=1;count<=60;count++) th1 >> lineoftext;
	for (k=0;k<=5;k++) {
		for (j=0;j<=5;j++) {
			for (l=0;l<=5;l++) {
				if ((i==0)||(j==0)||(k==0)||(l==0)) {
					data->tstkh[i][j][k][l]=0;
				}
				else if ((i==5)||(j==5)) {
             		data->tstkh[i][j][k][l] = infinity;
				}
				else if ((k==5)||(l==5)) {
             		data->tstkh[i][j][k][l] = 0;
				}
				else {
				    th1 >> lineoftext;
				    if (strcmp(lineoftext,".")){
					    data->tstkh[i][j][k][l] =(short) (floor(conversionfactor*(atof(lineoftext))+.5));
				    }
				    else data->tstkh[i][j][k][l] = infinity;
				}
            //cout <<"stack "<<i<<" "<<j<<" "<<k<<" "<<l<<"  "<<data->tstkh[i][j][k][l]<<"\n";
			}
         //cin >>m;
		}
	}
}

 th1.close();
/* Read info from tstacki */
//add to the tstacki table the case where X (represented as 0) is looked up:

ti1.open(tstacki);
for (count=1;count<=46;count++) ti1 >> lineoftext;//get past text in file
for (i=0;i<=5;i++) {
	if ((i!=0)&&(i!=5)) for (count=1;count<=60;count++) ti1 >> lineoftext;
	for (k=0;k<=5;k++) {
		for (j=0;j<=5;j++) {
			for (l=0;l<=5;l++) {

				if ((i==0)||(j==0)||(k==0)||(l==0)) {
					data->tstki[i][j][k][l]=0;
				}
				else if ((i==5)||(j==5)) {
             		data->tstki[i][j][k][l] = infinity;
				}
				else if ((k==5)||(l==5)) {
             		data->tstki[i][j][k][l] = 0;
				}
				//else if ((k==5)||(l==5)) {
				//include "5", linker for intermolecular for case of flush ends
				//	data->tstki[i][j][k][l]=0;
				//}
				else {
				    ti1 >> lineoftext;
                //cout <<lineoftext<<"\n";

				    if (strcmp(lineoftext,".")){
					data->tstki[i][j][k][l] =(short) (floor (conversionfactor*(atof(lineoftext))+.5));
				    }
				    else data->tstki[i][j][k][l] = infinity;
				}
            //cout <<"stack "<<i<<" "<<j<<" "<<k<<" "<<l<<"  "<<data->tstki[i][j][k][l]<<"\n";
			}
         //cin >> m;
		}
	}
}

 ti1.close();
/* Read info from tstacki23 */
//add to the tstacki table the case where X (represented as 0) is looked up:

t23.open(tstacki23);
for (count=1;count<=46;count++) t23 >> lineoftext;//get past text in file
for (i=0;i<=5;i++) {
	if ((i!=0)&&(i!=5)) for (count=1;count<=60;count++) t23 >> lineoftext;
	for (k=0;k<=5;k++) {
		for (j=0;j<=5;j++) {
			for (l=0;l<=5;l++) {

				if ((i==0)||(j==0)||(k==0)||(l==0)) {
					data->tstki23[i][j][k][l]=0;
				}
				else if ((i==5)||(j==5)) {
             		data->tstki23[i][j][k][l] = infinity;
				}
				else if ((k==5)||(l==5)) {
             		data->tstki23[i][j][k][l] = 0;
				}
				//else if ((k==5)||(l==5)) {
				//include "5", linker for intermolecular for case of flush ends
				//	data->tstki[i][j][k][l]=0;
				//}
				else {
				    t23 >> lineoftext;
                //cout <<lineoftext<<"\n";

				    if (strcmp(lineoftext,".")){
					data->tstki23[i][j][k][l] =(short) (floor (conversionfactor*(atof(lineoftext))+.5));
				    }
				    else data->tstki23[i][j][k][l] = infinity;
				}
            //cout <<"stack "<<i<<" "<<j<<" "<<k<<" "<<l<<"  "<<data->tstki[i][j][k][l]<<"\n";
			}
         //cin >> m;
		}
	}
}

/* Read info from tstacki1n */
//add to the tstacki table the case where X (represented as 0) is looked up:

t1n.open(tstacki1n);
for (count=1;count<=46;count++) t1n >> lineoftext;//get past text in file
for (i=0;i<=5;i++) {
	if ((i!=0)&&(i!=5)) for (count=1;count<=60;count++) t1n >> lineoftext;
	for (k=0;k<=5;k++) {
		for (j=0;j<=5;j++) {
			for (l=0;l<=5;l++) {

				if ((i==0)||(j==0)||(k==0)||(l==0)) {
					data->tstki1n[i][j][k][l]=0;
				}
				else if ((i==5)||(j==5)) {
             		data->tstki1n[i][j][k][l] = infinity;
				}
				else if ((k==5)||(l==5)) {
             		data->tstki1n[i][j][k][l] = 0;
				}
				//else if ((k==5)||(l==5)) {
				//include "5", linker for intermolecular for case of flush ends
				//	data->tstki[i][j][k][l]=0;
				//}
				else {
				    t1n >> lineoftext;
                //cout <<lineoftext<<"\n";

				    if (strcmp(lineoftext,".")){
					data->tstki1n[i][j][k][l] =(short) (floor (conversionfactor*(atof(lineoftext))+.5));
				    }
				    else data->tstki1n[i][j][k][l] = infinity;
				}
            //cout <<"stack "<<i<<" "<<j<<" "<<k<<" "<<l<<"  "<<data->tstki[i][j][k][l]<<"\n";
			}
         //cin >> m;
		}
	}
}

/*	Read info from tloops */
tl1.open(tloop);
for (count=1;count<=3;count++)	tl1 >> lineoftext;//get past text in file
data->numoftloops=0;
tl1>>lineoftext;


for (count=1;count<=maxtloop&&!tl1.eof();count++){
	//cout << lineoftext;

	(data->numoftloops)++;
	strcpy(base,lineoftext);
	strcpy(base+1,"\0");
	data->tloop[data->numoftloops][0] = tonumi(base);
//cout << base << "\n";
//cout << data->tloop[data->numoftloops][0] << "\n";

	strcpy(base,lineoftext+1);
	strcpy(base+1,"\0");
	data->tloop[data->numoftloops][0] = data->tloop[data->numoftloops][0]+
		5*tonumi(base);
//cout << base << "\n";
//cout << data->tloop[data->numoftloops][0] << "\n";
	strcpy(base,lineoftext+2);
	strcpy(base+1,"\0");
	data->tloop[data->numoftloops][0] = data->tloop[data->numoftloops][0]+
		25*tonumi(base);
//cout << base << "\n";
//cout << data->tloop[data->numoftloops][0] << "\n";
	strcpy(base,lineoftext+3);
	strcpy(base+1,"\0");
	data->tloop[data->numoftloops][0] = data->tloop[data->numoftloops][0]+
		125*tonumi(base);
//cout << base << "\n";
//cout << data->tloop[data->numoftloops][0] << "\n";
	strcpy(base,lineoftext+4);
	strcpy(base+1,"\0");
	data->tloop[data->numoftloops][0] = data->tloop[data->numoftloops][0]+
		625*tonumi(base);
//cout << base << "\n";
//cout << data->tloop[data->numoftloops][0] << "\n";
	strcpy(base,lineoftext+5);
	strcpy(base+1,"\0");
	data->tloop[data->numoftloops][0] = data->tloop[data->numoftloops][0]+
		3125*tonumi(base);
//cout << base << "\n";
//cout << data->tloop[data->numoftloops][0] << "\n";
	tl1 >> temp;
	data->tloop[data->numoftloops][1] = (short) (floor (conversionfactor*temp+0.5));


	//cout << "key = "<<data->tloop[data->numoftloops][0]<<"\n";
	//cout << "bonus = "<<data->tloop[data->numoftloops][1]<<"\n";
//	cin >> j;

   tl1 >> lineoftext;
}

 tl1.close();

//Read the 2x2 internal loops
//key iloop22[a][b][c][d][j][l][k][m] =
//a j l b
//c k m d


in1.open(int22);
for (count=1;count<=340;count++) in1 >> lineoftext;//get past text in file

for (i=1;i<=36;i++) {//read each of 36 tables
	for (j=1;j<=39;j++) in1 >> lineoftext;//get past text in file
	strcpy(base,lineoftext);
	strcpy(base+1, "\0");
	a = tonumi(base);
	for (j=1;j<=3;j++) in1 >> lineoftext;
	strcpy(base, lineoftext);
	strcpy(base+1, "\0");
	b = tonumi(base);
	in1>>lineoftext;
	strcpy(base, lineoftext);
	strcpy(base+1, "\0");
	c = tonumi(base);
	for (j=1;j<=3;j++) in1 >> lineoftext;
	strcpy(base, lineoftext);
	strcpy(base+1, "\0");
	d = tonumi(base);
	for (j=1;j<=3;j++) in1 >> lineoftext;//get past text in file
	for (j=1;j<=4;j++) {
	    for (k=1;k<=4;k++) {
		for (l=1;l<=4;l++) {
		    for (m=1;m<=4;m++) {
			in1 >> temp;
			data->iloop22[a][b][c][d][j][l][k][m] = (short) (floor(conversionfactor*temp+0.5));

         //no longer need to store the reverse order at same time because
         //the tables contain redundancy:
			//data->iloop22[d][c][b][a][m][k][l][j] = floor(conversionfactor*temp+0.5);


			//cout << "a = "<<a<<" b= "<<b<<" c= "<<c<<" d = "<<d<<"\n";

			//cout << "w = "<<j<<" x= "<<l<<" y= "<<k<<" z= "<<m<<"\n";
			//cout << data->iloop22[a][b][c][d][j][l][k][m]<<"\n";

		    }
          //cin >> foo;
		}
	    }
	}
}

 in1.close();

//Read the 2x1 internal loop data
in2.open(int21);
for (i=1;i<=58;i++) in2 >> lineoftext; //get past text at top of file
for (i=1;i<=6;i++) { //read each row of tables
	for (e=1;e<=4;e++) {
		for (j=1;j<=66;j++) in2 >> lineoftext; //get past text in file
		in2 >> lineoftext;
		strcpy(base,lineoftext);
		strcpy(base+1,"\0");
		a = tonumi(base);
		for (j=1;j<=11;j++) in2 >> lineoftext; //get past text in file
		in2 >> lineoftext;
		strcpy(base,lineoftext);
		strcpy(base+1,"\0");
		b = tonumi(base);
		for (j=1;j<=35;j++) in2 >> lineoftext; //get past text in file
		for (c=1;c<=4;c++) {
			for (j=1;j<=6;j++) {
				switch (j) {
					case 1:
						f = 1;
						g = 4;
						break;
					case 2:
						f = 2;
						g = 3;
						break;
					case 3:
						f = 3;
						g = 2;
						break;
					case 4:
						f = 4;
						g = 1;
						break;
					case 5:
						f = 3;
						g = 4;
						break;
					case 6:
						f = 4;
						g = 3;
						break;
				}
				for (d=1;d<=4;d++) {
					in2 >> temp;
					data->iloop21[a][b][c][d][e][f][g]=short (floor(conversionfactor*temp+0.5));
					//cout << a<<" "<<b<<" "<<c<<" "<<d<<" "<<e<<" "<<f<<" "<<g<<"\n";
               //cout << temp<<"\n";
               //cout << "conversionfactor*temp = "<<conversionfactor*temp<<"\n";
					//cout << data->iloop21[a][b][c][d][e][f][g]<<"\n";
					//cin >> temp;
				}
            //cin >> temp;
			}
		}
	}

}
 in2.close();

/*	Read info from triloops */
tri.open(triloop);
for (count=1;count<=3;count++)	tri >> lineoftext;//get past text in file
data->numoftriloops=0;
tri>>lineoftext;


for (count=1;count<=maxtloop&&!tri.eof();count++){

	//cout << lineoftext<<"\n";

	(data->numoftriloops)++;
	strcpy(base,lineoftext);
	strcpy(base+1,"\0");
	data->triloop[data->numoftriloops][0] = tonumi(base);
	strcpy(base,lineoftext+1);
	strcpy(base+1,"\0");
	data->triloop[data->numoftriloops][0] = data->triloop[data->numoftriloops][0]+
		5*tonumi(base);
	strcpy(base,lineoftext+2);
	strcpy(base+1,"\0");
	data->triloop[data->numoftriloops][0] = data->triloop[data->numoftriloops][0]+
		25*tonumi(base);
	strcpy(base,lineoftext+3);
	strcpy(base+1,"\0");
	data->triloop[data->numoftriloops][0] = data->triloop[data->numoftriloops][0]+
		125*tonumi(base);
	strcpy(base,lineoftext+4);
	strcpy(base+1,"\0");
	data->triloop[data->numoftriloops][0] = data->triloop[data->numoftriloops][0]+
		625*tonumi(base);
	tri >> temp;
	data->triloop[data->numoftriloops][1] = (short) (floor (conversionfactor*temp+0.5));

   //cout << data->triloop[data->numoftriloops][1]<< "  "<<data->triloop[data->numoftriloops][0]<<"\n";

   tri >> lineoftext;
}

 tri.close();

/*	Read info from hexaloops */
hl1.open(hexaloop);
for (count=1;count<=3;count++)	hl1 >> lineoftext;//get past text in file
data->numofhexaloops=0;
hl1>>lineoftext;


for (count=1;count<=maxtloop&&!hl1.eof();count++){

	//cout << lineoftext<<"\n";

	(data->numofhexaloops)++;
	strcpy(base,lineoftext);
	strcpy(base+1,"\0");
	data->hexaloop[data->numofhexaloops][0] = tonumi(base);
	strcpy(base,lineoftext+1);
	strcpy(base+1,"\0");
	data->hexaloop[data->numofhexaloops][0] = data->hexaloop[data->numofhexaloops][0]+
		5*tonumi(base);
	strcpy(base,lineoftext+2);
	strcpy(base+1,"\0");
	data->hexaloop[data->numofhexaloops][0] = data->hexaloop[data->numofhexaloops][0]+
		25*tonumi(base);
	strcpy(base,lineoftext+3);
	strcpy(base+1,"\0");
	data->hexaloop[data->numofhexaloops][0] = data->hexaloop[data->numofhexaloops][0]+
		125*tonumi(base);
	strcpy(base,lineoftext+4);
	strcpy(base+1,"\0");
	data->hexaloop[data->numofhexaloops][0] = data->hexaloop[data->numofhexaloops][0]+
		625*tonumi(base);
	strcpy(base,lineoftext+5);
	strcpy(base+1,"\0");
	data->hexaloop[data->numofhexaloops][0] = data->hexaloop[data->numofhexaloops][0]+
		3125*tonumi(base);
	strcpy(base,lineoftext+6);
	strcpy(base+1,"\0");
	data->hexaloop[data->numofhexaloops][0] = data->hexaloop[data->numofhexaloops][0]+
		15625*tonumi(base);
	strcpy(base,lineoftext+7);
	strcpy(base+1,"\0");
	data->hexaloop[data->numofhexaloops][0] = data->hexaloop[data->numofhexaloops][0]+
		78125*tonumi(base);
	hl1 >> temp;
	data->hexaloop[data->numofhexaloops][1] = (int) (floor (conversionfactor*temp+0.5));

   

   hl1 >> lineoftext;
}

 hl1.close();

/* Read info from coax */
//add to the stack table the case where X (represented as 0) is looked up:

//this is the array that keeps track of flush coaxial stacking (no intervening nucs)

//data arrangement of coax: data->coax[a][b][c][d]
//5'b-c3'
//3'a d5'
//this means the helix backbone is continuous between nucs b and c


co1.open(coax);
for (count=1;count<=42;count++) co1 >> lineoftext;//get past text in file
for (i=0;i<=5;i++) {
	if ((i!=0)&&(i!=5)) for (count=1;count<=60;count++) co1 >> lineoftext;
	for (k=0;k<=5;k++) {
		for (j=0;j<=5;j++) {
			for (l=0;l<=5;l++) {
				if ((i==0)||(j==0)||(k==0)||(l==0)) {
					data->coax[j][i][k][l]=0;
				}
            else if ((i==5)||(j==5)||(k==5)||(l==5)) {
             	data->coax[j][i][k][l] = infinity;
            }
				else {
					co1 >> lineoftext;

               //cout << lineoftext <<"end\n";

					if (strcmp(lineoftext,".")){
						data->coax[j][i][k][l] =short (floor(conversionfactor*(atof(lineoftext))+.5));
					}
					else data->coax[j][i][k][l] = infinity;

               //cin >> a;
				}
            //cout << j << " "<<i<<" "<<k<<" "<<l<<"  "<<data->coax[j][i][k][l]<<"\n";
			}
		}
      //cin>>foo;
	}
}

 co1.close();

/* Read info from tstackcoax */
//add to the tstackh table the case where X (represented as 0) is looked up:

//data arrangement of tstackcoax:
//5'a-c -> strand continues into stack
//3'b-d -> strand does not continue to stack
//pair between a-b, c-d is a mismatch


co2.open(tstackcoax);
for (count=1;count<=46;count++) co2 >> lineoftext;//get past text in file
for (i=0;i<=5;i++) {
	if (!(i==0||i==5)) for (count=1;count<=60;count++) co2 >> lineoftext;
	for (k=0;k<=5;k++) {
		for (j=0;j<=5;j++) {
			for (l=0;l<=5;l++) {
				if ((i==0)||(j==0)||(k==0)||(l==0)||(i==5)||(j==5)||(k==5)||(l==5)) {
					data->tstackcoax[i][j][k][l]=0;
				}
				else {
				    co2 >> lineoftext;
				    if (strcmp(lineoftext,".")){
					    data->tstackcoax[i][j][k][l] =short (floor(conversionfactor*(atof(lineoftext))+.5));
				    }
				    else data->tstackcoax[i][j][k][l] = infinity;
				}
			}
		}
	}
}

 co2.close();

/* Read info from coaxstack */
//add to the tstackh table the case where X (represented as 0) is looked up:

//data arrangement of coaxstack:
//5'a-c ->strand contnues into stack
//3'b d ->strand does not continue to stack
//pair between a-b, mismatch between c-d
//backbone is discontinuous between b and d


co3.open(coaxstack);
for (count=1;count<=46;count++) co3 >> lineoftext;//get past text in file
for (i=0;i<=5;i++) {
	if (!(i==0||i==5)) for (count=1;count<=60;count++) co3 >> lineoftext;
	for (k=0;k<=5;k++) {
		for (j=0;j<=5;j++) {
			for (l=0;l<=5;l++) {
				if ((i==0)||(j==0)||(k==0)||(l==0)||(i==5)||(j==5)||(k==5)||(l==5)) {
					data->coaxstack[i][j][k][l]=0;
				}
				else {
				    co3 >> lineoftext;
				    if (strcmp(lineoftext,".")){
					    data->coaxstack[i][j][k][l] =short (floor(conversionfactor*(atof(lineoftext))+.5));
				    }
				    else data->coaxstack[i][j][k][l] = infinity;
				}
			}
		}
	}
}

 co3.close();


/* Read info from tstack */
//this is the terminal mismatch data used in intermolecular folding
//add to the tstack table the case where X (represented as 0) is looked up.
//also add the case where 5 (the intermolecular linker) is looked up,
//this is actually a dangling end, not a terminal mismatch.

st2.open(tstack);
for (count=1;count<=46;count++) st2 >> lineoftext;//get past text in file
for (i=0;i<=5;i++) {
	if ((i!=0)&&(i!=5)) for (count=1;count<=60;count++) st2 >> lineoftext;
	for (k=0;k<=5;k++) {
		for (j=0;j<=5;j++) {
			for (l=0;l<=5;l++) {

				if ((i==0)||(j==0)||(k==0)||(l==0)) {
					data->tstack[i][j][k][l]=0;
				}
            else if ((i==5)||(j==5)) {
             	data->tstack[i][j][k][l] = infinity;

            }
				else if ((k==5)||(l==5)) {
				//include "5", linker for intermolecular for case of flush ends
            	if ((k==5)&&(l==5)) {//flush end
						data->tstack[i][j][k][l]=0;
               }
               else if (k==5) {//5' dangling end
               	//look up number for dangling end
               	data->tstack[i][j][k][l] = data->dangle[i][j][l][2];
               }
               else if (l==5) {//3' dangling end
               	data->tstack[i][j][k][l] = data->dangle[i][j][k][1];
               }
				}
				else {
				    st2 >> lineoftext;
				    if (strcmp(lineoftext,".")){
					data->tstack[i][j][k][l] =short (floor (conversionfactor*(atof(lineoftext))+.5));
				    }
				    else data->tstack[i][j][k][l] = infinity;
				}
            //cout <<"stack "<<i<<" "<<j<<" "<<k<<" "<<l<<"  "<<data->tstki[i][j][k][l]<<"\n";
			}
         //cin >> m;
		}
	}
}

 st2.close();

/* Read info from tstackm */
//add to the tstackm table the case where X (represented as 0) is looked up:

tsm.open(tstackm);
for (count=1;count<=46;count++) tsm >> lineoftext;//get past text in file
for (i=0;i<=5;i++) {
	if (i!=0) for (count=1;count<=60;count++) tsm >> lineoftext;
	for (k=0;k<=5;k++) {
		for (j=0;j<=5;j++) {
			for (l=0;l<=5;l++) {
				if ((i==0)||(j==0)||(k==0)||(l==0)) {
					data->tstkm[i][j][k][l]=0;
				}
				else if ((i==5)||(j==5)) {
             	data->tstkm[i][j][k][l] = infinity;

				}
				else if ((k==5)||(l==5)) {
					//include "5", linker for intermolecular for case of flush ends
            		if ((k==5)&&(l==5)) {//flush end
						data->tstkm[i][j][k][l]=0;
					}
					else if (k==5) {//5' dangling end
               			//look up number for dangling end
               			data->tstkm[i][j][k][l] = data->dangle[i][j][l][2]+penalty2(i,j,data);

						
					}
					else if (l==5) {//3' dangling end
               			data->tstkm[i][j][k][l] = data->dangle[i][j][k][1]+penalty2(i,j,data);
					}
				}
				else {
				    tsm >> lineoftext;
				    if (strcmp(lineoftext,".")){
					    data->tstkm[i][j][k][l] =short (floor(conversionfactor*(atof(lineoftext))+.5));
				    }
				    else data->tstkm[i][j][k][l] = infinity;
				}
            //cout <<"stack "<<i<<" "<<j<<" "<<k<<" "<<l<<"  "<<data->tstkh[i][j][k][l]<<"\n";
			}
         //cin >>m;
		}
	}
}

 tsm.close();

//data arrangement for 1x1 loop tables iloop11[a][b][c][d][e][f]:
//abc
//def


//Read the 1x1 internal loop data
//encode the data like:  abc
//                       def where b-e is a mismatch
i11.open(int11);
for (i=1;i<=58;i++) i11 >> lineoftext; //get past text at top of file
for (i=1;i<=6;i++) { //read each row of table
	if (i==1) {
    	a = 1;
      d = 4;
   }
   else if (i==2) {
    	a = 2;
      d = 3;
   }
   else if (i==3) {
    	a = 3;
      d = 2;
   }
   else if (i==4) {
    	a = 4;
      d = 1;
   }
   else if (i==5) {
    	a = 3;
      d = 4;
   }
   else {
   	a = 4;
      d = 3;
   }
	for (j=1;j<=114;j++) i11 >> lineoftext;//get past text
   for (b=1;b<=4;b++) {
   	for (j=1;j<=6;j++) {
      	if (j==1) {
    			c = 1;
      		f = 4;
   		}
   		else if (j==2) {
    			c = 2;
      		f = 3;
   		}
   		else if (j==3) {
    			c = 3;
      		f = 2;
   		}
   		else if (j==4) {
    			c = 4;
      		f = 1;
   		}
   		else if (j==5) {
    			c = 3;
      		f = 4;
   		}
   		else {
   			c = 4;
      		f = 3;
   		}
   		for (e=1;e<=4;e++) {
         	i11 >> temp;
            data->iloop11[a][b][c][d][e][f]=short (floor(conversionfactor*temp+0.5));

         }

      }
   }
}

 i11.close();


return 1;


}


//add the factor from SHAPE calculation
//This pseudo-energy was calculated when the file was loaded (see structure.cpp).
//The pseudo energy is applied twice for each nuc in interior pair and once for each nuc in terminal pair.
inline int SHAPEend(int i, structure *ct) {
	if (ct->shaped) return (int) ct->SHAPE[i];
	return 0;
}


void push(stackstruct *stack,int a,int b,int c,int d)
{
(stack->sp)++;
if (stack->sp >= STACKSIZE)	{
	cout<<"\n\nThe STACKSIZE is overflowed!!!The program need to be terminated\n";
	return;
}
stack->stk[stack->sp][0]= a;
stack->stk[stack->sp][1]= b;
stack->stk[stack->sp][2]= c;
stack->stk[stack->sp][3]= d;

}


void pull(stackstruct *stack,int *i,int *j,int *open,int *null,int *stz)
{
	if (stack->sp==0) {
		*stz = 1;
		return;
	}
	else {
	*stz = 0;
	*i = stack->stk[stack->sp][0];
	*j = stack->stk[stack->sp][1];
	*open= stack->stk[stack->sp][2];
	*null= stack->stk[stack->sp][3];
	stack->sp--;

	}
}

//calculate the energy of stacked base pairs

//oriented by:
//5' i ip 3'
//   |  |
//   j jp

int erg1(int i,int j,int ip,int jp,structure *ct, datatable *data)
{

		int energy;

		 if ((i==(ct->numofbases))||(j==((ct->numofbases)+1))) {
      	//this is not allowed because n and n+1 are not cavalently attached
         energy = infinity;
      }
		else {
      	energy = data->stack[(ct->numseq[i])][(ct->numseq[j])]
				[(ct->numseq[ip])][(ct->numseq[jp])]+data->eparam[1];
		//if (ct->shaped) {
			energy+=SHAPEend(i,ct);
			energy+=SHAPEend(j,ct);
			energy+=SHAPEend(ip,ct);
			energy+=SHAPEend(jp,ct);
		//}

      }
		return energy;
}


//calculate the energy of a bulge/internal loop
//where i is paired to j; ip is paired to jp; ip > i; j > jp
int erg2(int i,int j,int ip,int jp,structure *ct, datatable *data,
	char a, char b)
{

	int energy,size,size1,size2,loginc, lopsid, energy2,count,k;
	/* size,size1,size2 = size of a loop
		energy = energy calculated
		loginc = the value of a log used in large hairpin loops
	*/
   	if (((i<=(ct->numofbases))&&(ip>(ct->numofbases)))||((
      	jp<=(ct->numofbases))&&(j>(ct->numofbases)))) {
         //A loop cannot contain the ends of the sequence
         
         return infinity;
      }


	
      size1 = ip-i-1;
		size2 = j - jp - 1;

      if ((a>0)||(b>0)) {
      	if ((a&DOUBLE)||(b&DOUBLE)) return infinity;//the loop contains a nuc that
      		//should be double stranded
      	else if ((a&INTER)) {
         	//the loop is actually between two strands (ie: intermolecular)

             	if (size2>1) {//free energy is that of two terminal mismatches
               	//and the intermolecular initiation
                  energy = data->init + data->tstack[ct->numseq[i]][ct->numseq[j]]
                  	[ct->numseq[i+1]][ct->numseq[j-1]] +
                     data->tstack[ct->numseq[jp]][ct->numseq[ip]]
                  	[ct->numseq[jp+1]][ct->numseq[ip-1]];
               }
               else if (size2==1) {//find the best terminal mismatch and terminal
               	//stack free energies combination

                  energy = data->init + data->tstack[ct->numseq[i]][ct->numseq[j]]
                  		[ct->numseq[i+1]][ct->numseq[j-1]] +
                     	erg4 (jp,ip,ip-1,2,ct,data,false)+penalty(jp,ip,ct,data);
                  energy2 = data->init + data->tstack[ct->numseq[jp]][ct->numseq[ip]]
                  		[ct->numseq[jp+1]][ct->numseq[ip-1]] +
                     	erg4 (i,j,i+1,1,ct,data,false)+penalty(i,j,ct,data);
                  energy = min (energy,energy2);
                  //if ((ct->numseq[i+1]!=5)&&(ct->numseq[ip-1]!=5)) {
                     //now consider if coaxial stacking is better:
                     energy2 = data->init + data->tstackcoax[ct->numseq[jp]]
                     	[ct->numseq[ip]][ct->numseq[jp+1]][ct->numseq[ip-1]]
                        + data->coaxstack[ct->numseq[jp+1]][ct->numseq[ip-1]]
                        [ct->numseq[j]][ct->numseq[i]]+penalty(i,j,ct,data)+penalty(jp,ip,ct,data);
                     energy = min(energy,energy2);
                     energy2 = data->init + data->tstackcoax[ct->numseq[jp]]
                     	[ct->numseq[ip]][ct->numseq[j-1]][ct->numseq[ip-1]]
                        + data->coaxstack[ct->numseq[j-1]][ct->numseq[ip-1]]
                        [ct->numseq[j]][ct->numseq[i]]+penalty(i,j,ct,data)+penalty(jp,ip,ct,data);
                     energy = min(energy,energy2);
                  //}
               }
               else if (size2==0) {//just have dangling ends or flush stacking
               	energy = data->init + erg4 (jp,ip,ip-1,2,ct,data,false) +
                    	erg4 (i,j,i+1,1,ct,data,false)+penalty(i,j,ct,data)+penalty(jp,ip,ct,data);
                  energy2 = data->init + data->coax[ct->numseq[ip]][ct->numseq[jp]]
                  	[ct->numseq[j]][ct->numseq[i]]+penalty(i,j,ct,data)+penalty(jp,ip,ct,data);
                  energy = min(energy,energy2);
               }


         		return energy;
      	}
         else if (b&INTER) {
                  	//the loop is actually between two strands (ie: intermolecular)

             	if (size1>1) {//free energy is that of two terminal mismatches
               	//and the intermolecular initiation
                  energy = data->init + data->tstack[ct->numseq[i]][ct->numseq[j]]
                  	[ct->numseq[i+1]][ct->numseq[j-1]] +
                     data->tstack[ct->numseq[jp]][ct->numseq[ip]]
                  	[ct->numseq[jp+1]][ct->numseq[ip-1]];
               }
               else if (size1==1) {//find the best terminal mismatch and terminal
               	//stack free energies combination

                  energy = data->init + data->tstack[ct->numseq[i]][ct->numseq[j]]
                  		[ct->numseq[i+1]][ct->numseq[j-1]] +
                        erg4 (ip,jp,jp+1,1,ct,data,false)+penalty(ip,jp,ct,data);
                  energy2 = data->init + data->tstack[ct->numseq[jp]][ct->numseq[ip]]
                  		[ct->numseq[jp+1]][ct->numseq[ip-1]] +
                        erg4 (i,j,j-1,2,ct,data,false)+penalty(i,j,ct,data);

                  energy = min (energy,energy2);
                  //if ((ct->numseq[i+1]!=5)&&(ct->numseq[ip-1]!=5)) {
                     //now consider if coaxial stacking is better:
                     energy2 = data->init + data->tstackcoax[ct->numseq[i]]
                     	[ct->numseq[j]][ct->numseq[i+1]][ct->numseq[j-1]]
                        + data->coaxstack[ct->numseq[i+1]][ct->numseq[j-1]]
                        [ct->numseq[ip]][ct->numseq[jp]]+penalty(i,j,ct,data)+penalty(jp,ip,ct,data);
                     energy = min(energy,energy2);
                     energy2 = data->init + data->tstackcoax[ct->numseq[i]]
                     	[ct->numseq[j]][ct->numseq[ip-1]][ct->numseq[j-1]]
                        + data->coaxstack[ct->numseq[ip-1]][ct->numseq[j-1]]
                        [ct->numseq[ip]][ct->numseq[jp]]+penalty(i,j,ct,data)+penalty(jp,ip,ct,data);
                     energy = min(energy,energy2);
                  //}
               }
               else if (size1==0) {//just have dangling ends or flush stacking
               	energy = data->init + erg4 (jp,ip,jp+1,1,ct,data,false) +
                    	erg4 (i,j,j-1,2,ct,data,false)+penalty(i,j,ct,data)+penalty(jp,ip,ct,data);
                  energy2 = data->init + data->coax[ct->numseq[j]][ct->numseq[i]]
                  	[ct->numseq[ip]][ct->numseq[j]]+penalty(i,j,ct,data)+penalty(jp,ip,ct,data);
                  energy = min(energy,energy2);
               }


         		return energy;

         }
      }


      //a typical internal or bulge loop:
		//size1 = ip-i-1;
		//size2 = j - jp - 1;
		if (size1==0||size2==0) {//bulge loop


			size = size1+size2;
			if (size==1) {
				count = 1;
				energy = data->stack[ct->numseq[i]][ct->numseq[j]]
						[ct->numseq[ip]][ct->numseq[jp]]
						+ data->bulge[size] + data->eparam[2];
				if (size1==1)  {
					
					//count the number of alternative bulges that exist:
					
					k = i;
					while (ct->numseq[k]==ct->numseq[i+1]) {
						count++;
						k--;
					}
					k=ip;
					while(ct->numseq[k]==ct->numseq[i+1]) {
						count++;
						k++;
					}
					//give bonus to C adjacent to single C bulge
					if (ct->numseq[i+1]==2&&count>1) energy+= data->singlecbulge;
					
				}
				
				else {
					//size2 == 1
					
					//count the number of alternative bulges that exist:
					
					k = jp;
					while (ct->numseq[k]==ct->numseq[jp+1]) {
						count++;
						k--;
					}
					k=j;
					while(ct->numseq[k]==ct->numseq[jp+1]) {
						count++;
						k++;
					}
					//give bonus to C adjacent to single C bulge
					if (ct->numseq[j-1]==2&&count>1) energy+= data->singlecbulge;
					
				}
				//apply a correction for the number of equivalent states because 
					//the bulge can move to adjacent sites
				energy-= (int) ( (data->RT)*conversionfactor* log ((double) count));
			}
			else if (size>30) {

				loginc = int((data->prelog)*log(double ((size)/30.0)));
				energy = data->bulge[30] + loginc + data->eparam[2];
				energy = energy + penalty(i,j,ct,data) + penalty(jp,ip,ct,data);

			}
			else {
         		energy = data->bulge[size] + data->eparam[2];
				energy = energy + penalty(i,j,ct,data) + penalty(jp,ip,ct,data);
			}
		}
		else {//internal loop
			size = size1 + size2;
			lopsid = abs(size1-size2);

			if (size>30) {

				loginc = int((data->prelog)*log((double ((size))/30.0)));
				if (size1==1||size2==1) {
            		energy = data->tstki1n[ct->numseq[i]][ct->numseq[j]]
						[ct->numseq[i+1]][ct->numseq[j-1]] +
						data->tstki1n[ct->numseq[jp]][ct->numseq[ip]]
						[ct->numseq[jp+1]][ct->numseq[ip-1]] +
						data->inter[30] + loginc + data->eparam[3] +
						min(data->maxpen,(lopsid*
						data->poppen[min(2,min(size1,size2))]));

				}

				else {
					energy = data->tstki[ct->numseq[i]][ct->numseq[j]]
						[ct->numseq[i+1]][ct->numseq[j-1]] +
						data->tstki[ct->numseq[jp]][ct->numseq[ip]]
						[ct->numseq[jp+1]][ct->numseq[ip-1]] +
						data->inter[30] + loginc + data->eparam[3] +
						min(data->maxpen,(lopsid*
						data->poppen[min(2,min(size1,size2))]));
				}
			}
			else if ((size1==2)&&(size2==2))//2x2 internal loop
			    energy = data->iloop22[ct->numseq[i]][ct->numseq[ip]]
					[ct->numseq[j]][ct->numseq[jp]]
					[ct->numseq[i+1]][ct->numseq[i+2]]
					[ct->numseq[j-1]][ct->numseq[j-2]];


			else if ((size1==1)&&(size2==2)) {//2x1 internal loop
				energy = data->iloop21[ct->numseq[i]][ct->numseq[j]][ct->numseq[i+1]]
					[ct->numseq[j-1]][ct->numseq[jp+1]][ct->numseq[ip]][ct->numseq[jp]];


			}
			else if ((size1==2)&&(size2==1)) {//1x2 internal loop
				energy = data->iloop21[ct->numseq[jp]][ct->numseq[ip]][ct->numseq[jp+1]]
					[ct->numseq[ip-1]][ct->numseq[i+1]][ct->numseq[j]][ct->numseq[i]];

			}

			else if (size==2) //a single mismatch
				
				energy = data->iloop11[ct->numseq[i]][ct->numseq[i+1]][ct->numseq[ip]]
					[ct->numseq[j]][ct->numseq[j-1]][ct->numseq[jp]];
			else if (size1==1||size2==1) { //this loop is lopsided
         	//this is a 1xn loop:
				energy = data->tstki1n[ct->numseq[i]][ct->numseq[j]][ct->numseq[i+1]][ct->numseq[j-1]] +
						data->tstki1n[ct->numseq[jp]][ct->numseq[ip]][ct->numseq[jp+1]][ct->numseq[ip-1]] +
						data->inter[size] + data->eparam[3] +
					min(data->maxpen,(lopsid*data->poppen[min(2,min(size1,size2))]));
			}
        

			else if ((size1==2&&size2==3)||(size1==3&&size2==2)) {
			//this is a 2x3 loop
				energy = data->tstki23[ct->numseq[i]][ct->numseq[j]][ct->numseq[i+1]][ct->numseq[j-1]] +
					data->tstki23[ct->numseq[jp]][ct->numseq[ip]][ct->numseq[jp+1]][ct->numseq[ip-1]] +
					data->inter[size] + data->eparam[3] +
					min(data->maxpen,(lopsid*data->poppen[min(2,min(size1,size2))]));


			}
			else {
         	


         		energy = data->tstki[ct->numseq[i]][ct->numseq[j]][ct->numseq[i+1]][ct->numseq[j-1]] +
					data->tstki[ct->numseq[jp]][ct->numseq[ip]][ct->numseq[jp+1]][ct->numseq[ip-1]] +
					data->inter[size] + data->eparam[3] +
					min(data->maxpen,(lopsid*data->poppen[min(2,min(size1,size2))]));
			}
		}

		return energy;
}


//calculate the energy of a hairpin loop:
int erg3(int i,int j,structure *ct, datatable *data,char dbl)
{
int energy,size,loginc,count,key,k;
	/* size,size1,size2 = size of a loop
		energy = energy calculated
		loginc = the value of a log used in large hairpin loops
	*/


		if (dbl&DOUBLE) return infinity;//the loop contains a base that should be
      										//double stranded

      else if (dbl&INTER) {//intermolecular interaction
      	//intermolecular "hairpin" free energy is that of intermolecular
         //	initiation plus the stacked mismatch

         energy = data->init + min(data->tstack[ct->numseq[i]][ct->numseq[j]]
         	[ct->numseq[i+1]][ct->numseq[j-1]],erg4(i,j,i+1,1,ct,data,false))+penalty(i,j,ct,data);

         return energy;
      }


   	if ((i<=(ct->numofbases))&&(j>(ct->numofbases))) {
      	//A hairpin cannot contain the ends of the sequence
         energy = infinity;
         return energy;
      }

		size = j-i-1;



		if (size>30) {

			loginc = int((data->prelog)*log((double ((size))/30.0)));



			energy = data->tstkh[ct->numseq[i]][ct->numseq[j]]
				[ct->numseq[i+1]][ct->numseq[j-1]]
				+ data->hairpin[30]+loginc+data->eparam[4];
		}
		else if (size<3) {
      		energy = data->hairpin[size] + data->eparam[4];
				if (ct->numseq[i]==4||ct->numseq[j]==4) energy = energy+6;
		}
		else if (size==4) {
			
			key = (ct->numseq[j])*3125 + (ct->numseq[i+4])*625 +
				(ct->numseq[i+3])*125 + (ct->numseq[i+2])*25+(ct->numseq[i+1])*5+(ct->numseq[i]);
			for (count=1;count<=data->numoftloops;count++) {
					if (key==data->tloop[count][0]) return data->tloop[count][1];
			}
			energy = data->tstkh[ct->numseq[i]][ct->numseq[j]]
				[ct->numseq[i+1]][ct->numseq[j-1]]
				+ data->hairpin[size] + data->eparam[4];
		}
		else if (size==3) {
			
			key = (ct->numseq[j])*625 +
				(ct->numseq[i+3])*125 + (ct->numseq[i+2])*25+(ct->numseq[i+1])*5+(ct->numseq[i]);
			for (count=1;count<=data->numoftriloops;count++) {
				if (key==data->triloop[count][0]) return data->triloop[count][1];
			}
			
			energy =	data->hairpin[size] + data->eparam[4]
         	+penalty(i,j,ct,data);
		}
		else if (size==6) {
			key = (ct->numseq[j])*78125 + (ct->numseq[i+6])*15625 + (ct->numseq[i+5])*3125 
				+ (ct->numseq[i+4])*625 +
				(ct->numseq[i+3])*125 + (ct->numseq[i+2])*25+(ct->numseq[i+1])*5+(ct->numseq[i]);
			for (count=1;count<=data->numofhexaloops;count++) {
				if (key==data->hexaloop[count][0]) return data->hexaloop[count][1];
			}

			energy = data->tstkh[ct->numseq[i]][ct->numseq[j]]
				[ct->numseq[i+1]][ct->numseq[j-1]]
				+ data->hairpin[size] + data->eparam[4];
		}

		else {
			energy = data->tstkh[ct->numseq[i]][ct->numseq[j]]
				[ct->numseq[i+1]][ct->numseq[j-1]]
				+ data->hairpin[size] + data->eparam[4];
		}




		//check for GU closeure preceded by GG
      if (ct->numseq[i]==3&&ct->numseq[j]==4) {
      	if ((i>2&&i<ct->numofbases)||(i>ct->numofbases+2))
       		if (ct->numseq[i-1]==3&&ct->numseq[i-2]==3) {

         		energy = energy + data->gubonus;
            	


         	}
      }

      //check for an oligo-c loop
      
      for (k=1;(k<=size);k++) {
       	if (ct->numseq[i+k] != 2) return energy;//this is not an oligo-C loop
      }
      //this is a poly c loop so penalize
      if (size==3) return (energy + data->c3);
      else return (energy + data->cint + size*data->cslope);
      

}



//calculate the energy of a dangling end:
int erg4(int i,int j,int ip,int jp,structure *ct, datatable *data, bool lfce)
{
int energy;
//dangling base
		// jp = 1 => 3' dangle
		// jp = 2 => 5' dangle



      if (lfce) return infinity;//stacked nuc should be double stranded

	  //commented out 11/8/99
      //if (ip==5) return 0;//dangling nuc is an intermolecular linker

		energy = data->dangle[ct->numseq[i]][ct->numseq[j]][ct->numseq[ip]][jp];
		return energy;
}



//this function calculates whether a terminal pair i,j requires the end penalty
int penalty2(int i,int j, datatable *data) {


   if (i==4||j==4)
   	return data->auend;
   else return 0;//no end penalty


}


//When considering mismatch at the end of a helix, consult this function to check
//	whether the nucs are required to pair
int checknp(bool lfce1,bool lfce2) {

	if (lfce1||lfce2) return infinity;
	else return 0;
}

//this function calculates the free energy for coaxial stacking of pair i-j onto ip-jp
//	require that the backbone continues directly between j and ip without a nucleotide in
//	a canonical pair
//k indicates the intervening nuc in intervening mismatch that is not sandwiched between
//	j and ip	
int ergcoax(int i, int j, int ip, int jp, int k, structure *ct, datatable *data) {
	

	if (ip==j+1) {
		//flush stacking

		return data->coax[ct->numseq[i]][ct->numseq[j]][ct->numseq[ip]][ct->numseq[jp]];

	}

	else if (k>0) {
		//coaxial stacking with an intervening mismatch
		if (k==i-1) {
			return data->tstackcoax[ct->numseq[j]][ct->numseq[i]][ct->numseq[j+1]][ct->numseq[i-1]] + 
				data->coaxstack[ct->numseq[j+1]][ct->numseq[k]][ct->numseq[ip]][ct->numseq[jp]];	
			
		}
		else { //if (k==jp+1) {
			return data->tstackcoax[ct->numseq[jp]][ct->numseq[ip]][ct->numseq[jp+1]][ct->numseq[ip-1]] +
				data->coaxstack[ct->numseq[j]][ct->numseq[i]][ct->numseq[j+1]][ct->numseq[k]];	
		}
		//else {
			//some error -- give message
			//errmsg(100,1);

		//}

	}
	else return infinity;
	


}
//these functions calculate the free energy for coaxial stacking of pair i-j onto ip-jp
//	require that the backbone continues directly between j and ip without a nucleotide in
//	a canonical pair
//k (determined by the function called) indicates the intervening nuc in intervening mismatch that is not sandwiched between
//	j and ip	
int ergcoaxflushbases(int i, int j, int ip, int jp, structure *ct, datatable *data) {
	//flush stacking, k==0

		return data->coax[ct->numseq[i]][ct->numseq[j]][ct->numseq[ip]][ct->numseq[jp]];

}
int ergcoaxinterbases1(int i, int j, int ip, int jp, structure *ct, datatable *data) {
	//k==i-1
	return data->tstackcoax[ct->numseq[j]][ct->numseq[i]][ct->numseq[j+1]][ct->numseq[i-1]] + 
				data->coaxstack[ct->numseq[j+1]][ct->numseq[i-1]][ct->numseq[ip]][ct->numseq[jp]];

}

int ergcoaxinterbases2(int i, int j, int ip, int jp, structure *ct, datatable *data) {
	//k==jp+1
	return data->tstackcoax[ct->numseq[jp]][ct->numseq[ip]][ct->numseq[jp+1]][ct->numseq[ip-1]] +
				data->coaxstack[ct->numseq[j]][ct->numseq[i]][ct->numseq[j+1]][ct->numseq[jp+1]];

}
	

//this function calculates the free energy for coaxial stacking of pair i-j onto ip-jp
//	require that the backbone continues directly between j and ip without a nucleotide in
//	a canonical pair


//Note that this form of the function takes the sequences, not the number of the nuc
//	in i,j,ip,and jp.	
int ergcoaxflushbases(int i, int j, int ip, int jp, datatable *data) {
	
		

		return data->coax[i][j][ip][jp];

}



//this function calculates the free energy for coaxial stacking of pair i-j onto ip-jp
//	require that the backbone continues directly between j and ip without a nucleotide in
//	a canonical pair
//k indicates the intervening nuc in intervening mismatch that is not sandwiched between
//	j and ip
//l is the nuc sandwiched between j and ip


//this requires a backbone like:
// j-l-ip
// | . |
// i-k jp
//i.e. a discontinuity between k and jp


//Note that this form of the function takes the sequences, not the number of the nuc
//	in i,j,ip,jp,k and l.	
int ergcoaxinterbases1(int i, int j, int ip, int jp, int k, int l, datatable *data) {
	
		//coaxial stacking with an intervening mismatch
		
			return data->tstackcoax[j][i][l][k] + 
				data->coaxstack[l][k][ip][jp];	
			


}


//this function calculates the free energy for coaxial stacking of pair i-j onto ip-jp
//	require that the backbone continues directly between j and ip without a nucleotide in
//	a canonical pair
//k indicates the intervening nuc in intervening mismatch that is not sandwiched between
//	j and ip
//l is the nuc sandwiched between j and ip


//this requires a backbone like:
// j-l-ip
// | . |
// i k-jp
//i.e. a discontinuity between i and jp


//Note that this form of the function takes the sequences, not the number of the nuc
//	in i,j,ip,jp,k and l.	
int ergcoaxinterbases2(int i, int j, int ip, int jp, int k, int l, datatable *data) {
	
		//coaxial stacking with an intervening mismatch
		
			return data->tstackcoax[jp][ip][k][l] + 
				data->coaxstack[j][i][l][k];	
			


}

//These functions are used by ergmulti to deconvolute the base pair code back into nucs
int decon1(int x) {
		
	return (int) (floor(((float)(x))/10)-1);

}
int decon2(int x) {

	return x - (decon1(x)+1)*10-1;
}



//This function will calculate the free energy of a multiloop starting at nuc i for 
//	structure #st
//This uses a recursive algorithm

#define js true //this switch can turn on and off the logarithmic dep. on 
				//unpaired nucleotides.  This is helpful for debugging

//simplemb, when true, indicates that the logarithmic dependence should be off.
//This sets the energy function equal to that used by the dynamic programming
//algorithms.

int ergmulti(int st, int ip, structure *ct, datatable *data, bool simplemb) {
	short int *element,**energy;
	short int i,count,b,c,size,j,minimum,au,current,first,recent,biggest;
	bool intermolecular;
	float average;

	//trace out the loop to learn the size:
	count = 0;
	i = ip;




	while(i!=ip||count==0) {
		i++;
		if (ct->basepr[st][i]) {
			i=ct->basepr[st][i];
		}

		count++;


	}
	
	
	//allocate element to store this info:
	element = new short int [count+3];
	energy = new short int *[count+3];
	for (i=0;i<count+3;i++) {
		energy[i]=new short int [count+3];
		for (j=0;j<count+3;j++)
			energy[i][j]=0;

	}
	biggest = 0;
	average = 0;
	b=0;//keep track of the number of unpaired nucs 
	c=0;//keep track of the number of helixes
	au=0;//keep track of the number of terminal AU/GU pairs
	intermolecular = false; //catch whether this multiloop contains the
							//intermolecular linker
							//if so, there are no multiloop penalties

	//record the info:

	element[0] = 10*(ct->numseq[ct->basepr[st][ip]]+1);
	element[0] = element[0] + ct->numseq[ip]+1;

	//if (ct->numseq[ct->basepr[st][ip]]==4) au++;
	current = 0;
	i = ip;
	count = 1;
	while (i!=ip||count==1) {
		i++;
		if (ct->basepr[st][i]) {
			
			if(c>0) {
				if(abs(current-recent)>biggest) biggest=abs(current-recent);	
				average = average + (float) (abs(current-recent));
			}
			
			else first = current;
			recent = current;
			current=0;
			
			

			c++;
			//store pairs in a code:
			element[count] = 10*(ct->numseq[i]+1);
			if (ct->numseq[i]==4) au++;
			i=ct->basepr[st][i];
			element[count]=element[count]+ct->numseq[i]+1;
			if (ct->numseq[i]==4) au++;
			
			



		}
		else {
			current++;
			b++;
			element[count] = ct->numseq[i];
			if (element[count]==5) intermolecular=true;


		}
		
		count++;



	}

	if(abs(current-first)>biggest) biggest=abs(current-first);
	average = average + (float) (abs(current-recent));

	average = average/((float) (c));
	

	element[count] = element[1];
	element[count+1] = element[2];
	count--;  //correct the count elements 
				//bcs the first helix was counted twice


	//start the minimization routine:
	
	for (size=2;size<=count;size++) {
		for (i=0;i+size<count+4;i++) {

			
			if (size==2) {
				//the interactions are either a nuc stacking 
				//onto a helix or flush coaxial stacking;

				//check for flush stacking first:
				if (element[i]>10&&element[i+1]>10) {
					//calculate flush stacking:
					
					energy[i][i+1]=ergcoaxflushbases(decon1(element[i]),
						decon2(element[i]),decon1(element[i+1]),
						decon2(element[i+1]),data);


				}
				else if (element[i]>10) {
					//3' dangling end:
					energy[i][i+1] = data->dangle[decon2(element[i])]
						[decon1(element[i])][element[i+1]][1];

				}
				else if (element[i+1]>10) {
					//5' dangling end:
					energy[i][i+1] = data->dangle[decon2(element[i+1])]
						[decon1(element[i+1])][element[i]][2];


				}
				else 
					energy[i][i+1]=0;

				



			}
			else if (size==3) {
				energy[i][i+2] = min(energy[i][i+1],energy[i+1][i+2]);
				
			


				if (element[i]<10&&element[i+1]>10&&element[i+2]<10) {
					//consider mismatch stack ono helix:
					energy[i][i+2]=min(energy[i][i+2],
						data->tstkm[decon2(element[i+1])][decon1(element[i+1])]
							[element[i+2]][element[i]]);

				}


				

			}

			else if (size==4) {


				energy[i][i+3]=min(energy[i][i+1]+energy[i+2][i+3],
					energy[i][i+2]);
				energy[i][i+3]=min(energy[i][i+3],energy[i+1][i+3]);

				//now the fragment is big enough for coaxial stacking with
				//	an intervening mismatch:
				if (element[i]>10&&element[i+1]<10&&element[i+2]>10&&element[i+3]<10) {

					energy[i][i+3] = min(energy[i][i+3],
						ergcoaxinterbases2(decon1(element[i]),decon2(element[i]),
							decon1(element[i+2]), decon2(element[i+2]), 
							element[i+3], element[i+1], data)); 


				}
				else if (element[i]<10&&element[i+1]>10&&element[i+2]<10&&
					element[i+3]>10) {

					energy[i][i+3] = min(energy[i][i+3],
						ergcoaxinterbases1(decon1(element[i+1]),decon2(element[i+1]),
							decon1(element[i+3]),decon2(element[i+3]),
							element[i],element[i+2],data));



				}


				

			}
			else {
				//now size > 4
				//energy[i][i+size-1]=0;
				for (j=i;j<i+size-1;j++) {

					energy[i][i+size-1]=min(energy[i][i+size-1],
						energy[i][j]+energy[j+1][i+size-1]);


				}

				




			}
		}
	}


	//now find the minimum combination: 
	//There are 4 phases of length = count to be checked:


	minimum = min(energy[0][count-1],energy[1][count]);
	minimum = min(minimum,energy[2][count+1]);
	minimum = min(minimum,energy[3][count+2]);


	
	
	//deallocate memory use:
	
	for (i=0;i<count+3;i++) delete[] energy[i];
	delete[] energy;
	delete[] element;


	//return the energy:
	
	if (intermolecular)

		return minimum+au*data->auend+data->init;
	
	
	
	
	//do not add strin term for simplemb==true
	if (((c%2)!=0)&&((b==0)||(b==1))&&(!simplemb)) minimum = minimum+data->strain;

	//do not use the assymetry for simplemb==true
	if (simplemb) average = 0;
	else if (average>2) average =2; 

	minimum = minimum+(short) (((float) data->mlasym)*average);
	
	

	//note that this function, used by efn2, has a logarithmic dependence in
	//unpaired nucleotides after 8.  This is hardwired below.
	if (b>8&&js&&!intermolecular&&!simplemb) 
		return minimum + data->efn2a + c*data->efn2c +
		8*(data->efn2b) 
			+ int(11.*log(double(((double (b))/8.))) + 0.5) +
au*data->auend;

       else 

	
		minimum = minimum + data->efn2a + b*data->efn2b + c*data->efn2c
			+au*data->auend;

	return minimum;
}

// Energyout: writes to file a list of energys calculated by efn2
void energyout(structure *ct,char *energyfile) {
int i;
ofstream out(energyfile);

for (i=1;i<=ct->numofstructures;i++)
	out << "Structure: "<<i<<"   Energy = "<<(float (ct->energy[i])/10)<<"   \n";

}


//This function will calculate the free energy of an exterior loop 
//  for structure #st
//This uses a recursive algorithm



int ergexterior(int st, structure *ct, datatable *data) {
	short int *element,**energy;
	short int i,count,size,j,minimum,au,helices;
	bool intermolecular;

	//trace out the loop to learn the size:
	count = 0;
	//start counting at 1
	i = 0;

	intermolecular = false; //keep track as to whether the bimolecular linker is found
							//	if so, add the intermolecular initiation free energy

	helices  = 0;

	while(i!=ct->numofbases) {
		i++;
		if (ct->basepr[st][i]) {
			i=ct->basepr[st][i];
			helices++;
		}

		count++;


	}
	
	//check for empty structure and return 0 if empty
	if (helices==0) return 0;
	
	//allocate element to store this info:
	element = new short int [count];
	energy = new short int *[count];
	for (i=0;i<count;i++) {
		energy[i]=new short int [count];
		for (j=0;j<count;j++)
			energy[i][j]=0;

	}

	
	au=0;//keep track of the number of terminal AU/GU pairs
	

	//record the info:

	

	

	i = 0;
	count = 0;
	while (i!=ct->numofbases) {
		i++;
		if (ct->basepr[st][i]) {
			
			
			//store pairs in a code:
			element[count] = 10*(ct->numseq[i]+1);
			if (ct->numseq[i]==4) au++;
			i=ct->basepr[st][i];
			element[count]=element[count]+ct->numseq[i]+1;
			if (ct->numseq[i]==4) au++;
			


		}
		else {
			
			element[count] = ct->numseq[i];
			if (element[count]==5) intermolecular=true;
			


		}
		
		count++;



	}

	



	//start the minimization routine:
	
	for (size=2;size<=count;size++) {
		for (i=0;i+size<=count;i++) {

			
			if (size==2) {
				//the interactions are either a nuc stacking 
				//onto a helix or flush coaxial stacking;

				//check for flush stacking first:
				if (element[i]>10&&element[i+1]>10) {
					//calculate flush stacking:
					
					energy[i][i+1]=ergcoaxflushbases(decon1(element[i]),
						decon2(element[i]),decon1(element[i+1]),
						decon2(element[i+1]),data);


				}
				else if (element[i]>10) {
					//3' dangling end:
					energy[i][i+1] = data->dangle[decon2(element[i])]
						[decon1(element[i])][element[i+1]][1];

				}
				else if (element[i+1]>10) {
					//5' dangling end:
					energy[i][i+1] = data->dangle[decon2(element[i+1])]
						[decon1(element[i+1])][element[i]][2];


				}
				else 
					energy[i][i+1]=0;

				



			}
			else if (size==3) {
				energy[i][i+2] = min(energy[i][i+1],energy[i+1][i+2]);
				
				if (element[i]<10&&element[i+1]>10&&element[i+2]<10) {
					//consider mismatch stack ono helix:
					energy[i][i+2]=min(energy[i][i+2],
						data->tstack[decon2(element[i+1])][decon1(element[i+1])]
							[element[i+2]][element[i]]);

				}


				

			}

			else if (size==4) {

			

				energy[i][i+3]=min(energy[i][i+1]+energy[i+2][i+3],
					energy[i][i+2]);
				energy[i][i+3]=min(energy[i][i+3],energy[i+1][i+3]);

				//now the fragment is big enough for coaxial stacking with
				//	an intervening mismatch:
				if (element[i]>10&&element[i+1]<10&&element[i+2]>10&&element[i+3]<10) {

					energy[i][i+3] = min(energy[i][i+3],
						ergcoaxinterbases2(decon1(element[i]),decon2(element[i]),
							decon1(element[i+2]), decon2(element[i+2]), 
							element[i+3], element[i+1], data));


				}
				else if (element[i]<10&&element[i+1]>10&&element[i+2]<10&&
					element[i+3]>10) {

					energy[i][i+3] = min(energy[i][i+3],
						ergcoaxinterbases1(decon1(element[i+1]),decon2(element[i+1]),
							decon1(element[i+3]),decon2(element[i+3]),
							element[i],element[i+2],data));



				}


				

			}
			else {
				//now size > 4
				//energy[i][i+size-1]=0;
				for (j=i;j<i+size-1;j++) {

					energy[i][i+size-1]=min(energy[i][i+size-1],
						energy[i][j]+energy[j+1][i+size-1]);


				}

				




			}
		}
	}



	
	minimum = energy[0][count-1];

	if (intermolecular) minimum = minimum + data->init;
	//deallocate memory use:
	
	for (i=0;i<count;i++) delete[] energy[i];
	delete[] energy;
	delete[] element;


	//return the energy:
	

	return minimum+au*data->auend;
}




void de_allocate (int **v,int i) {//this function deallocates the memory used
												//in an array
int a;

	for (a=0;a<i;a++) {
   	delete[] v[a];
   }
   delete[] v;
}


void de_allocate (short int **v,int i) {//this function deallocates the memory used
												//in an array
int a;

	for (a=0;a<i;a++) {
   	delete[] v[a];
   }
   delete[] v;
}

void de_allocate (bool **v,int i) {//this function deallocates the memory used
												//in an array
int a;

	for (a=0;a<i;a++) {
   	delete[] v[a];
   }
   delete[] v;
}

void stackclass::allocate_stack() {
	short i;
		
	stackenergy =new integersize [maximum];
	stack=new short int *[maximum];
	for (i=0;i<maximum;i++) stack[i] = new short int [4];

}



stackclass::stackclass(short int stacksize) {
	maximum = stacksize;
	size = 0;
	allocate_stack();

}
	
	


bool stackclass::pull(short int *i,short int *j, short int *open, 
		integersize *energy, short int *pair) {
		
	if (size==0) return false;
	else {
		size--;
		*i = stack[size][0];
		*j = stack[size][1];
		*open = stack[size][2];
		*energy = stackenergy[size];
		*pair = stack[size][3];
		return true;
			
	}

}
	
void stackclass::push(short int i,short int j, short int open, 
		integersize energy, short int pair){

	short k;

	if (size == maximum) {
		//allocate more space:
		stackclass *temp;
		temp = new stackclass(maximum);
		for (k=0;k<maximum;k++) {
			temp->push(stack[k][0],stack[k][1],stack[k][2],stackenergy[k],stack[k][3]);
		}
		delete_array();
		maximum = 2*maximum;

		allocate_stack();
		for (k=0;k<(maximum/2);k++) {
			temp->pull(&(stack[k][0]),&(stack[k][1]),&(stack[k][2]),&(stackenergy[k]),&(stack[k][3]));
		}

		delete temp;
		
	}
		
	stack[size][0] = i;
	stack[size][1] = j;
	stack[size][2] = open;
	stackenergy[size] = energy;
	stack[size][3] = pair;
	size++;

}
	
void stackclass::delete_array() {
	short i;
		
	for (i=0;i<maximum;i++) delete[] stack[i];
	delete[] stack;

	delete[] stackenergy;


}


stackclass::~stackclass() {
		
	delete_array();

}

//read is used to read data from a save file
void read(ifstream *out,short *i) {
	
	out->read((char *) i,sizeof(*i));
}

void read(ifstream *out,bool *i) {

	out->read((char *) i,sizeof(*i));
}

void read(ifstream *out,int *i) {

	out->read((char *) i,sizeof(*i));
}

void read(ifstream *out,char *i) {
	int length;

	out->read((char *) (&length), sizeof(length));
	out->read(i,length);
}

void readsinglechar(ifstream *out,char *i) {
	
	out->read(i,sizeof(char));
}

void read(ifstream *out,float *i) {

	out->read((char *) i,sizeof(*i));
}

//write is used to write data to a save file
void write(ofstream *out,short *i) {
	
	out->write((char *) i,sizeof(*i));
}

void write(ofstream *out,bool *i) {

	out->write((char *) i,sizeof(*i));
}

void write(ofstream *out,int *i) {

	out->write((char *) i,sizeof(*i));
}

void write(ofstream *out,char *i) {
	int length;

	length = strlen(i)+1;

	out->write((char *) (&length), sizeof(length));
	out->write(i,strlen(i)+1);
}
void write(ofstream *out,float *i) {

	out->write((char *) i,sizeof(*i));
}

void writesinglechar(ofstream *out, char *i) {
	out->write(i,sizeof(char));

}

void writehelixfile(char *filename,structure *ct,int StructureNumber) {
	//write a helix file that can be read by XRNA
	int i,count;
	ofstream out;

	out.open(filename);

	i=1;
	while (i<=ct->numofbases) {
		if (ct->basepr[StructureNumber][i]>i) {
			//found a base pair
			out << i << " " << ct->basepr[StructureNumber][i] << " ";

			//determine length of helix
			count = 1;
			while (ct->basepr[StructureNumber][i+1]==ct->basepr[StructureNumber][i]-1) {
				i++;
				count++;

			}
			out << count << "\n";
			i++;
			

		}
		else i++;


	}

	

}


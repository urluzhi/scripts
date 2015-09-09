/**********************************************************************
* Making tables of free energy and enthalpy  for version 3.? and 2.3
* revised from Function opendat
* Function opens data files to read thermodynamic data
* Created: Jul. 2005				Modified: Dec. 2005
* Copy Right: 	Zhi John Lu
************************************************************************/


/*
1.tstack.dat23/dh23: 
	(from dangle.dat23/dh23)

*/
#include<cstdio>
#include<fstream>
#include<iostream>
#include<cstring>
using namespace std;

#define infinity 140000;
int tonumi(char *base);	
void outspace(ifstream &input, ofstream &output);
void tstack(char *in1, char *in2);
	//a.out tstack.dat/dh dangle.dat23/dh23
void tstackh(char *in);
    //a.out -h tstack.dat/dh
void tstacki(char *in);
	//a.out tstacki.datc/tstacki.dhc/tstacki23.dhc/datc
void tloops(char *in1, char *in2, char *out);
	// a.out tstackh.dat23/dh23 tloop.dat/dh tloop.dat23/dh23
void int11(char *outf);
	//a.out int11.dhc
void int21(char *output);
	//a.out int21.dhc/datc
void int22(char *input1, char *input2);
	//a.out tstacki.dat23/dh23  int22.dat/dh
int main (int argc, char *argv[])
{
	if (strcmp(argv[1],"int11.dhc")==0) int11(argv[1]);
	if (strcmp(argv[1],"-h")==0)	tstackh(argv[2]);
	if (strcmp(argv[1],"tstacki.dhc")==0  || strcmp(argv[1],"tstacki.datc")==0 || strcmp(argv[1],"tstacki23.datc")==0  || strcmp(argv[1],"tstacki23.dhc")==0 )	tstacki(argv[1]);
	if (strcmp(argv[1], "int21.dhc")==0  || strcmp(argv[1],"int21.datc")==0 )	int21(argv[1]);
	if (argc==4)
	if (strcmp(argv[3],"tloop.dat23")==0 || strcmp (argv[3],"tloop.dh23")==0 )
		{
    	cout<< argv[3]<<flush;	
		tloops(argv[1], argv[2],argv[3]);
		}
	 if(argc==3)
		{
		tstack(argv[1], argv[2]);
		if(strcmp(argv[2], "int22.dat")==0 || strcmp(argv[2], "int22.dh")==0)
		{
		cout<< argv[2]<<"\n"<<flush;
		int22(argv[1],argv[2]);
		}
		}
		
		return 1;


}

void int11(char *outf)
{
	int count,i,j,a,b,c,d,e,f;
	char base[7],lineoftext[100],inf[20];
	char numstr[10];
	float temp,asym;
	ifstream in;
	ofstream out;
	if (strcmp(outf,"int11.dhc") == 0) strcpy (inf, "int11.dh");
	out.open(outf);
	in.open(inf);
	float iloop11[5][5][5][5][5][5]; 
	
//data arrangement for 1x1 loop tables iloop11[a][b][c][d][e][f]:
//abc
//def

//encode the data like:  abc
//                       def where b-e is a mismatch
for (i=1;i<=58;i++)
{
 	outspace(in,out);
	in >> lineoftext; //get past text at top of file
	out << lineoftext;	
}
for (i=1;i<=6;i++) { //read each row of table
 	outspace(in,out);
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
	for (j=1;j<=114;j++)
	{
 	outspace(in,out);
 	in >> lineoftext;//get past text
	out << lineoftext;	
    }
	for (b=1;b<=4;b++) {
 	outspace(in,out);
   	for (j=1;j<=6;j++) {
 	outspace(in,out);
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
 			outspace(in,out);
         	in >> temp;
			temp= -10.53;
			
			if (a==4||d==4)    temp+=5.025;
			if (c==4||f==4)    temp+=5.025; 
 			if ((a==1||a==3)&&( b==4&&e==4)) 	temp-=3.39;
			if ((f==1||f==3)&&( b==4&&e==4)) 	temp-=3.39;
			if ( b==3&&e==3)   temp-=7.87;
			sprintf(numstr,"%.1f",temp);
			out << numstr;
         }

      }
   }
}

	out.close();
	in.close();

}

void int22(char *input1, char *input2)
{
	int count,a,b,c,d,m,i,j,k,l;
	char base[7],lineoftext[100],output[20];
	char numstr[10];
	float temp,asym;
	ifstream in1;
	ifstream in2;
	ofstream out;
	in1.open(input1);
	in2.open(input2);
	if (strcmp(input2,"int22.dat") == 0) strcpy (output, "int22.dat23");
 	if (strcmp(input2,"int22.dh")==0 ) strcpy(output,"int22.dh23");
	out.open(output);
	float iloop22[5][5][5][5][5][5][5][5]; 
	float tstacki[5][5][5][5];

 /*	read info from stacki */
 /*key:
5'   i  k  3'
3'   j  l  5'
A=1, C=2, G=3, U=4
*/
for (count=1;count<=46;count++)
{
	in1 >> lineoftext;//get past text in file
}
for (i=1;i<=4;i++) {
	for (count=1;count<=60;count++) 
	{
		in1 >> lineoftext;
	}

	for (k=1;k<=4;k++) {
		for (j=1;j<=4;j++) {
			for (l=1;l<=4;l++) {
					    in1 >> lineoftext;
						if (strcmp(lineoftext,".")){
									tstacki[i][j][k][l]=atof(lineoftext);
				    		}
				    	else 
						{
							tstacki[i][j][k][l]=infinity;
						}
	
}
}
}
}
		


//Read and write a the 2x2 internal loop data
		//Read the 2x2 internal loops
//key iloop22[a][b][c][d][j][l][k][m] =
//a j l b
//c k m d


for (count=1;count<=340;count++) 
	{
	outspace(in2,out);
	in2 >> lineoftext;//get past text in file
	out << lineoftext;
	}
for (i=1;i<=36;i++) {//read each of 36 tables
	outspace(in2,out);
	for (j=1;j<=39;j++)
    {
		outspace(in2,out);
		in2 >> lineoftext;//get past text in file
		out << lineoftext;	
	}
	strcpy(base,lineoftext);
	strcpy(base+1, "\0");
	a = tonumi(base);
	for (j=1;j<=3;j++)
	{
	 outspace(in2,out);
	 in2 >> lineoftext;
	 out << lineoftext;
	}
	strcpy(base, lineoftext);
	strcpy(base+1, "\0");
	b = tonumi(base);
	outspace(in2,out);
	in2>>lineoftext;
	out<<lineoftext;
	strcpy(base, lineoftext);
	strcpy(base+1, "\0");
	c = tonumi(base);
	for (j=1;j<=3;j++)
	{
	 outspace(in2,out);
	 in2 >> lineoftext;
	 out << lineoftext;
	}
	strcpy(base, lineoftext);
	strcpy(base+1, "\0");
	d = tonumi(base);
	for (j=1;j<=3;j++) 
	{
	outspace(in2,out);
	in2 >> lineoftext;//get past text in file
	out << lineoftext;
	}
	for (j=1;j<=4;j++) {
		outspace(in2,out);
	    for (k=1;k<=4;k++) {
		outspace(in2,out);
		for (l=1;l<=4;l++) {
		outspace(in2,out);
		    for (m=1;m<=4;m++) {
			outspace(in2,out);
			in2 >> temp;
			if (strcmp(input2, "int22.dat")==0) 	{temp=4.9; }
			if (strcmp(input2, "int22.dh")==0)		{temp=0;  }
			iloop22[a][b][c][d][j][l][k][m] = tstacki[a][c][j][k] + tstacki[d][b][m][l]+ temp;
			sprintf(numstr,"%.1f",iloop22[a][b][c][d][j][l][k][m]);
			out<< numstr;
		

    		    }
		}
	    }
	}
}
			in1.close();
			in2.close();
			out.close();

}


void int21(char *output)
{
	int count,a,b,c,d,e,f,g,i,j,k,l;
	char base[7],lineoftext[100],input[20];
	char numstr[10];
	float temp,asym;
	ifstream in2;
	ofstream out;
	if(strcmp(output, "int21.dhc")==0)		strcpy(input, "int21.dh");
	else if(strcmp(output, "int21.datc")==0) strcpy(input, "int21.dat");
	out.open(output);
	in2.open(input);
	float tint21[5][5][5][5][5][5][5]; 

	



 /*key:

       c
5'   a   f  3'
3'   b   g  5'
      d e

A=1, C=2, G=3, U=4
*/

//Read and write a the 2x1 internal loop data
for (i=1;i<=58;i++)
{
	outspace(in2,out);
 	in2 >> lineoftext; //get past text at top of file
	out << lineoftext;
}
for (i=1;i<=6;i++) { //read each row of tables
	outspace(in2,out);
	for (e=1;e<=4;e++) {
		outspace(in2,out);
		for (j=1;j<=66;j++)
		{
			outspace(in2,out);
			in2 >> lineoftext; //get past text in file
			out << lineoftext;
		}
		outspace(in2,out);
		in2 >> lineoftext;
		out << lineoftext;
		strcpy(base,lineoftext);
		strcpy(base+1,"\0");
		a = tonumi(base);
		for (j=1;j<=11;j++)
		{
		outspace(in2,out);
		in2 >> lineoftext; //get past text in file
		out << lineoftext;
		}
		outspace(in2, out);
		in2 >> lineoftext;
		out << lineoftext;
		outspace(in2, out);
		strcpy(base,lineoftext);
		strcpy(base+1,"\0");
		b = tonumi(base);
		for (j=1;j<=35;j++)
		{
			outspace(in2,out);
			in2 >> lineoftext; //get past text in file
			out << lineoftext;
		}
		for (c=1;c<=4;c++) {
			outspace(in2,out);
			for (j=1;j<=6;j++) {
				outspace(in2,out);
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
					outspace(in2,out);
					in2 >> temp;
					
					if(strcmp(output, "int21.dhc")==0)
					{	
						//initiation and asymmetry
						temp=0.29+3.205;
						//AU/GU closure penalty
						if (a==4||b==4) temp+=5.025;
						if (f==4||g==4) temp+=5.025;
						//UU bonus
						if (c==4&& (d==4||e==4))    temp-=10.15;
						//G bonus
						if ( (c==3&&(d==3||d==1)) ||(e==3&&(c==3||c==1)) )  temp-=5.76;
						else if ( ((a==2||a==4)&&c==1&&d==3)|| ((g==2||g==4)&&e==1&&c==3) ) temp-=5.76;
		
						//special experimental data
						if (a==3 && f==3)
						{
							if(c==1 && d==1 && e==1) 		temp=6.9;
							else if(c==1 && d==1 && e==3) 	temp=-5.0;
							else if(c==1 && d==3 && e==1)	temp=2.4;
							else if(c==2 && d==2 && e==4)	temp=1.3;
						}
						else if(a==2 && f==2)
						{
							if(c==1 && d==1 && e==1)		temp=2.6;
							else if(c==1 && d==1 && e==3)	temp=-2.3;
							else if(c==1 && d==2 && e==2)	temp=-2.5;
							else if(c==1 && d==3 && e==1)	temp=2.5;
							else if(c==2 && d==1 && e==1)	temp=1.6;
							else if(c==2 && d==1 && e==4)	temp=6.4;
							else if(c==2 && d==2 && e==2)	temp=7.8;
							else if(c==3 && d==1 && e==1)	temp=-1.6;
							else if(c==3 && d==1 && e==3)	temp=-3.1;
							else if(c==3 && d==3 && e==1)	temp=-4.1;
							else if(c==4 && d==2 && e==2)	temp=5.8;
							else if(c==4 && d==2 && e==4)	temp=-6.9;
							else if(c==4 && d==4 && e==2)	temp=-5.8;
							else if(c==4 && d==4 && e==4)	temp=-7.0;
							
						}
					}
				
					
					else if(strcmp(output, "int21.datc")==0)
					{
						//initiation and asymmetry
						temp=1.68+0.58;
						//AU/GU closure penalty
						if (a==4||b==4) temp+=0.74;
						if (f==4||g==4) temp+=0.74;
						//UU bonus
						if (c==4&& (d==4||e==4))    temp-=0.77;
						//G bonus
						if ( (c==3&&(d==3||d==1)) ||(e==3&&(c==3||c==1)) )  temp-=1.15;
						else if ( ((a==2||a==4)&&c==1&&d==3)|| ((g==2||g==4)&&e==1&&c==3) ) temp-=1.15;
						//special experimental data
						if (a==3 && f==3)
						{
							if(c==1 && d==1 && e==1) 		temp=2.5;
							else if(c==1 && d==1 && e==3) 	temp=1.2;
							else if(c==1 && d==3 && e==1)	temp=2.1;
							else if(c==2 && d==2 && e==4)	temp=1.9;
						}
						else if(a==2 && f==2)
						{
							if(c==1 && d==1 && e==1)		temp=2.5;
							else if(c==1 && d==1 && e==3)	temp=0.8;
							else if(c==1 && d==2 && e==2)	temp=1.7;
							else if(c==1 && d==3 && e==1)	temp=1.1;
							else if(c==2 && d==1 && e==1)	temp=2.3;
							else if(c==2 && d==1 && e==4)	temp=2.5;
							else if(c==2 && d==2 && e==2)	temp=2.5;
							else if(c==3 && d==1 && e==1)	temp=1.7;
							else if(c==3 && d==1 && e==3)	temp=1.2;
							else if(c==3 && d==3 && e==1)	temp=0.8;
							else if(c==4 && d==2 && e==2)	temp=2.2;
							else if(c==4 && d==2 && e==4)	temp=1.7;
							else if(c==4 && d==4 && e==2)	temp=1.5;
							else if(c==4 && d==4 && e==4)	temp=1.4;
							
						}
			
					}
					temp=(floor(10*temp+0.5))/10;
					sprintf(numstr,"%.1f",temp);
					out<< numstr;
				}
			}
		}
	}

}
			outspace(in2,out);
			in2.close();
			out.close();
}


/*	Read and output info from tloops */

/*
key:
    c d 
   b   e
    a-f
*/


void tloops(char *in1,char *in2, char *out)
{
    int count,a,b,c,d,e,f,i,j,k,l;
	char base[7],lineoftext[100];
	float temp;
	char numstr[10];	
	ifstream tl1;
	ofstream tl2;
	ifstream tsh;
	tl1.open(in2);
	tl2.open(out);
	tsh.open(in1);
 	
	float tloop[5][5][5][5][5][5]; 
	float tstackh[5][5][5][5];

 /*	read info from stackh */
 /*key:
5'   i  k  3'
3'   j  l  5'
A=1, C=2, G=3, U=4
*/
for (count=1;count<=46;count++)
{
	tsh >> lineoftext;//get past text in file
}
for (i=1;i<=4;i++) {
	for (count=1;count<=60;count++) 
	{
		tsh >> lineoftext;
	}

	for (k=1;k<=4;k++) {
		for (j=1;j<=4;j++) {
			for (l=1;l<=4;l++) {
					    tsh >> lineoftext;
						if (strcmp(lineoftext,".")){
									tstackh[i][j][k][l]=atof(lineoftext);
				    		}
				    	else 
						{
							tstackh[i][j][k][l]=infinity;
						}
	
}
}
}
}
		
		tsh.close();







/* Read tloops.dat/dh and write tloops.dat23/dh23*/

	for (count=1;count<=2;count++)
	{
		tl1 >> lineoftext;//get past text in file
		tl2 << lineoftext<<" \t";
	}
	tl2<<"\n";
	tl1>>lineoftext;
	tl2<<lineoftext<<"\n";

for (count=1;!tl1.eof();count++){
	tl1>>lineoftext;
	tl2<<lineoftext<<" \t";
	
	strcpy(base,lineoftext);
	strcpy(base+1,"\0");
	a = tonumi(base);
	strcpy(base,lineoftext+1);
	strcpy(base+1,"\0");
	b = tonumi(base);
	strcpy(base,lineoftext+2);
	strcpy(base+1,"\0");
	c = tonumi(base);
	strcpy(base,lineoftext+3);
	strcpy(base+1,"\0");
	d = tonumi(base);
	strcpy(base,lineoftext+4);
	strcpy(base+1,"\0");
	e = tonumi(base);
	strcpy(base,lineoftext+5);
	strcpy(base+1,"\0");
	f = tonumi(base);

	tl1 >> temp;
	
	if (strcmp(out,"tloop.dat23")==0) 
	temp = tstackh[a][f][b][e]-2.0 +4.9 ;
	if (strcmp(out,"tloop.dh23")==0)
	temp = tstackh[a][f][b][e]-4.0 ; 
    tl2 << temp<<"\n"; 
}

	tl1.close();
	tl2.close();
	
}

int tonumi(char *base)	{
int	a;
if (!strcmp(base,"A")||!strcmp(base,"B")) (a = 1);
else if(!strcmp(base,"C")||!strcmp(base,"Z")) (a = 2);
else if(!strcmp(base,"G")||!strcmp(base,"H")) (a = 3);
else if(!strcmp(base,"U")||!strcmp(base,"V")) (a = 4);
else if(!strcmp(base,"T")||!strcmp(base,"W")) (a = 4);
else (a=0);  //this is for others, like X
return a;
}



/* Read info from tstack */
//add to the tstacki table :
void tstacki(char *in)
{
	int count,i,j,k,l;
	char lineoftext[100];
	float temp;
	char numstr[10];	
	char tstackallf[20], tstackhf[20];
	strcpy (tstackhf,in );
    if (strcmp(in,"tstacki.dhc")==0) strcpy (tstackallf, "tstacki.dh");
    if (strcmp(in,"tstacki23.dhc")==0) strcpy (tstackallf, "tstacki23.dh");
    if (strcmp(in,"tstacki.datc")==0) strcpy (tstackallf, "tstacki.dat");
	if (strcmp(in,"tstacki23.datc")==0) strcpy (tstackallf, "tstacki23.dat");
	ifstream tsa;
	ofstream tsh;
	tsa.open(tstackallf);
	tsh.open(tstackhf);
 
float tstackh[6][6][6][6]; 
 /*	read info from stackall */
 /*key:
5'   i  k  3'
3'   j  l  5'
A=1, C=2, G=3, U=4
*/
for (count=1;count<=46;count++)
{
 	outspace(tsa,tsh);
	tsa >> lineoftext;//get past text in file
	tsh << lineoftext;
}
outspace(tsa,tsh);
for (i=1;i<=4;i++) {
outspace(tsa,tsh);
	for (count=1;count<=60;count++) 
	{
		outspace(tsa,tsh);
		tsa >> lineoftext;
		tsh << lineoftext;
	}

	for (k=1;k<=4;k++) {
		outspace(tsa,tsh);
		for (j=1;j<=4;j++) {
			outspace(tsa,tsh);
			for (l=1;l<=4;l++) {
						outspace(tsa,tsh);
					    tsa >> lineoftext;
						if (strcmp(lineoftext,".")){
									tstackh[i][j][k][l]=atof(lineoftext);
									if (strcmp(in,"tstacki23.datc")==0 || strcmp(in,"tstacki.datc")==0) {
										if (i==4 ||j==4) 	tstackh[i][j][k][l] = 0.74 ;
										else                tstackh[i][j][k][l] = 0 ;
									}
									else {
										if (i==4 ||j==4) 	tstackh[i][j][k][l] = 5.03 ;
										else                tstackh[i][j][k][l] = 0 ;
									}
									if (strcmp(in,"tstacki.datc")==0)
									{
										if ( (k==1&&l==3)) 	tstackh[i][j][k][l] = tstackh[i][j][k][l]-0.83;
										if ( (k==3&&l==1)) 	tstackh[i][j][k][l] = tstackh[i][j][k][l]-1.03;
										if ( (k==3&&l==3)) 	tstackh[i][j][k][l] = tstackh[i][j][k][l]-1.03;
										if ( (k==4&&l==4)) 	tstackh[i][j][k][l] = tstackh[i][j][k][l]-0.61;

									}
									if (strcmp(in,"tstacki.dhc")==0)
									{
										if ( (k==1&&l==3)) 	tstackh[i][j][k][l] = tstackh[i][j][k][l]-3.38;
										if ( (k==3&&l==1)) 	tstackh[i][j][k][l] = tstackh[i][j][k][l]-7.63;
										if ( (k==3&&l==3)) 	tstackh[i][j][k][l] = tstackh[i][j][k][l]+2.84;
										if ( (k==4&&l==4)) 	tstackh[i][j][k][l] = tstackh[i][j][k][l]-5.80;

									}
									
									if (strcmp(in,"tstacki23.datc")==0)
									{
										if ( i==2||i==4)
										{
										if ( k==1&&l==3) 	tstackh[i][j][k][l] = tstackh[i][j][k][l]-0.51;
										if ( k==3&&l==1) 	tstackh[i][j][k][l] = tstackh[i][j][k][l]-1.13;
										}
										if ( (i==1||i==3)&&(k==3&&l==1))	tstackh[i][j][k][l] = tstackh[i][j][k][l]-1.18;
										if ( k==3&&l==3) 	tstackh[i][j][k][l] = tstackh[i][j][k][l]-0.74;
										if ( k==4&&l==4) 	tstackh[i][j][k][l] = tstackh[i][j][k][l]-0.35;
									}

    								if (strcmp(in,"tstacki23.dhc")==0)
									{
										if ( i==2||i==4)
										{
										if ( k==1&&l==3) 	tstackh[i][j][k][l] = tstackh[i][j][k][l]-5.67;
										if ( k==3&&l==1) 	tstackh[i][j][k][l] = tstackh[i][j][k][l]-8.64;
										}
										if ( (i==1||i==3)&&(k==3&&l==1))	tstackh[i][j][k][l] = tstackh[i][j][k][l]-10.87;
										if ( k==3&&l==3) 	tstackh[i][j][k][l] = tstackh[i][j][k][l]-9.03;
										if ( k==4&&l==4) 	tstackh[i][j][k][l] = tstackh[i][j][k][l]-6.42;
									}
									sprintf (numstr,"%.1f",tstackh[i][j][k][l]);
									tsh << numstr;
				    		}
				    	else tsh << lineoftext;
outspace(tsa,tsh);
	
}
}
}
}
		tsa.close();
		tsh.close();
}





/* Read info from tstack */
//add to the tstackh table the case where X (represented as 0) is looked up:
void tstackh(char *in)
{
	if (strcmp(in,"tstack.dat")==0 || strcmp(in,"tstack.dh")==0)
	{
	int count,i,j,k,l;
	char lineoftext[100];
	float temp;
	char numstr[10];	
	char tstackallf[20], tstackhf[20];
	strcpy (tstackallf, in);
    if (strcmp(in,"tstack.dat")==0)	strcpy (tstackhf, "tstackh.datc");
    if (strcmp(in,"tstack.dh")==0)	strcpy (tstackhf, "tstackh.dhc");
	ifstream tsa;
	ofstream tsh;
	tsa.open(tstackallf);
	tsh.open(tstackhf);
 
float tstackh[6][6][6][6]; 
 /*	read info from stackall */
 /*key:
5'   i  k  3'
3'   j  l  5'
A=1, C=2, G=3, U=4
*/
for (count=1;count<=46;count++)
{
 	outspace(tsa,tsh);
	tsa >> lineoftext;//get past text in file
	tsh << lineoftext;
}
outspace(tsa,tsh);
for (i=1;i<=4;i++) {
outspace(tsa,tsh);
	for (count=1;count<=60;count++) 
	{
		outspace(tsa,tsh);
		tsa >> lineoftext;
		tsh << lineoftext;
	}

	for (k=1;k<=4;k++) {
		outspace(tsa,tsh);
		for (j=1;j<=4;j++) {
			outspace(tsa,tsh);
			for (l=1;l<=4;l++) {
						outspace(tsa,tsh);
					    tsa >> lineoftext;
						if (strcmp(lineoftext,".")){
									tstackh[i][j][k][l]=atof(lineoftext);
									if (i==4 ||j==4)
										tstackh[i][j][k][l] = tstackh[i][j][k][l] +3.7 ;
									if ( (k==3&&l==1)||(k==4&&l==4) )
										tstackh[i][j][k][l] = tstackh[i][j][k][l]-5.8;
									sprintf (numstr,"%.1f",tstackh[i][j][k][l]);
									tsh << numstr;
				    		}
				    	else tsh << lineoftext;
outspace(tsa,tsh);
	
}
}
}
}
		tsa.close();
		tsh.close();
}
}


void tstack(char *in1, char *in2)
	{

	
//cout<< in1;cout.flush();
//cout<< in2<<flush;
	if (strcmp(in1,"tstack.dat")==0 || strcmp(in1,"tstack.dh")==0)
	{
	/* Write info to tstack */
	//this is the terminal mismatch data used in intermolecular folding
	//add to the tstack table the case where X (represented as 0) is looked up.
	//also add the case where 5 (the intermolecular linker) is looked up,
	//this is actually a dangling end, not a terminal mismatch.
	int count,i,j,k,l;
	char lineoftext[100];
	float temp;
	char numstr[10];	
	
	char tstackf[20], danglef[20], tstack23[20];
	strcpy (tstackf, in1);
	strcpy (danglef, in2);
	strcpy (tstack23, in1);
	strcat (tstack23,"23");
	ifstream tst;
	ifstream dag;
	ofstream ts2;
	tst.open(tstackf);
	dag.open(danglef); 
	ts2.open(tstack23);
 
float dangle[6][6][6][3], tstack[6][6][6][6]; 
 /*	read info from dangle */
 //add to dangle the case where X (represented as 0) is looked up
/*key:
l=1:   i  k  3'
   	   j

l=2:   i  
       j k 5'
*/

for (l = 1;l <=2; l++){
	for (i = 0;i <=5; i++){
		if ((i!=0)&&(i!=5)) for (count=1;count <=60;count++) dag >> lineoftext;
		for (j=0;j<=5; j++) {
			for (k=0;k<=5; k++) {
				if ((i==0)||(j==0)||(k==0)) {
				    dangle[i][j][k][l] = 0;
				}
            else if ((i==5)||(j==5)||(k==5)) {
             	dangle[i][j][k][l] = 0;
            }
				else {
				    dag >> lineoftext;
                //cout << lineoftext<<"\n";
				    if (strcmp(lineoftext,".")){
					dangle[i][j][k][l] = atof(lineoftext);
				    }
				    else 
					{
					dangle[i][j][k][l]=infinity;
					}
				}
            //cout <<"dangle = "<<i<<" "<<j<<" "<<k<<" "<<l<<" "<<dangle[i][j][k][l]<<"\n";
			}
         //cin >> m;
		}
	}
}

/*read tstack.dat/dh and write tstack.dat23/tstack.dh23*/
/*key:

5' i k 3'
3' j l 5'
*/
for (count=1;count<=46;count++)
{
 	outspace(tst,ts2);
	tst >> lineoftext;//get past text in file
	ts2 << lineoftext;
}
outspace(tst,ts2);
for (i=1;i<=4;i++) {
outspace(tst,ts2);
	if ((i!=0)&&(i!=5)) for (count=1;count<=60;count++) 
	{
		outspace(tst,ts2);
		tst >> lineoftext;
		ts2 << lineoftext;
	}

		for (k=1;k<=4;k++) {
			outspace(tst,ts2);
		for (j=1;j<=4;j++) {
			outspace(tst,ts2);
			for (l=1;l<=4;l++) {
							outspace(tst,ts2);
						    tst >> lineoftext;
							if (strcmp(lineoftext,".")){
										tstack[i][j][k][l] = dangle[i][j][k][1] +dangle[i][j][l][2] ;
										sprintf (numstr,"%.1f",tstack[i][j][k][l]);
										ts2 << numstr;
				    			}
				    		else ts2 << lineoftext;
outspace(tst,ts2);
				
            //cout <<"stack "<<i<<" "<<j<<" "<<k<<" "<<l<<"  "<<tstki[i][j][k][l]<<"\n";
			}
         //cin >> m;
		}
	}
}
	
	tst.close();
	dag.close();
	ts2.close();

	}
}



void outspace(ifstream &input, ofstream &output)
{
	char ctemp;
	ctemp=(char)input.peek();
	while (ctemp=='\t' || ctemp=='\n' || ctemp==' ')
	{
		input.get(ctemp);
		output.put(ctemp);
		ctemp=(char)input.peek();
	}
}

//This program will read a 2x2 internal loop datatable with measured parameters
//	the user hardwires delta values
//	then the program outputs a 2x2 internal loop table with each position calculated


#include "Stdafx.h"
#include "algorithm.cpp"
#define conversionfactor 10



void main() {
 	int i,j,k,l,m,n,iloop22[6][6][6][6][6][6][6][6],a,b,c,d,len,count;
   int tenth,first,integer,size1,size2,stable,delta,gg;
   char input[100],output[100],temp[20],base[10],lineoftext[500];
   int penalty,bonus;
   int inc[6][6]={{0,0,0,0,0,0},{0,0,0,0,1,0},{0,0,0,1,0,0},{0,0,1,0,1,0},
	{0,1,0,1,0,0},{0,0,0,0,0,0}};
   int np,lp,mp,kp,writeinteger;
   int tstack[6][6][6][6];

   



   
   cout << "output 2x2 file: ";
   cin >> output;

   
   ofstream out(output);
   ifstream st2;
   st2.open("dnatstack.dat");


for (count=1;count<=46;count++) st2 >> lineoftext;//get past text in file
for (i=0;i<=5;i++) {
	if ((i!=0)&&(i!=5)) for (count=1;count<=60;count++) st2 >> lineoftext;
	for (k=0;k<=5;k++) {
		for (j=0;j<=5;j++) {
			for (l=0;l<=5;l++) {

				if ((i==0)||(j==0)||(k==0)||(l==0)) {
					tstack[i][j][k][l]=0;
				}
            else if ((i==5)||(j==5)) {
             	tstack[i][j][k][l] = infinity;

            }
				else if ((k==5)||(l==5)) {
				//include "5", linker for intermolecular for case of flush ends
            	if ((k==5)&&(l==5)) {//flush end
						tstack[i][j][k][l]=0;
               }
               else if (k==5) {//5' dangling end
               	//look up number for dangling end
               	tstack[i][j][k][l] = 0;
               }
               else if (l==5) {//3' dangling end
               	tstack[i][j][k][l] = 0;
               }
				}
				else {
				    st2 >> lineoftext;
				    if (strcmp(lineoftext,".")){
					tstack[i][j][k][l] =short (floor (conversionfactor*(atof(lineoftext))+.5));
				    }
				    else tstack[i][j][k][l] = infinity;
				}
            //cout <<"stack "<<i<<" "<<j<<" "<<k<<" "<<l<<"  "<<data->tstki[i][j][k][l]<<"\n";
			}
         //cin >> m;
		}
	}
}





//output a 2x2 internal loop table -- keep all of redundancies for Zuker
	out << "Data tables for symetric interior loops of size 4\n";
   out << "Free energies at 37 degrees for DNA\n";
   out << "Data arrangement:\n\n";
   out << "                                Y\n";
   out << "     ----------------------------------------------------------------\n";
   out << "(X)   A   A   A   A   C   C   C   C   G   G   G   G   T   T   T   T\n";
   out << "      A   C   G   T   A   C   G   T   A   C   G   T   A   C   G   T\n";
   out << "     ----------------------------------------------------------------\n";
   out << "                           5' ------> 3'\n";
   out << "                            A \\/ \\_/ A\n";
   out << "                            T /\\  |  T\n";
   out << "                           3' <------ 5'\n\n";
   out << "(AA)  .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .  \n";
   out << "(AC)  .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .  \n";
   out << "(AG)  .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .  \n";
   out << "(AT)  .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .  \n";
   out << "(CA)  .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .  \n";
   out << "(CC)  .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .  \n";
   out << "(CG)  .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .  \n";
   out << "(CT)  .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .  \n";
   out << "(GA)  .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .  \n";
   out << "(GC)  .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .  \n";
   out << "(GG)  .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .  \n";
   out << "(GT)  .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .  \n";
   out << "(TA)  .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .  \n";
   out << "(TC)  .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .  \n";
   out << "(TG)  .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .  \n";
   out << "(TT)  .   .   .   .   .   .   .   .   .   .   .   .   .   .   .   .  \n\n";

   for (i=1;i<=6;i++) {
    	for (j=1;j<=6;j++) {

      	if (i==1) {
         	a = 1;
            c = 4;
         }
         else if (i==2) {
          	a = 2;
            c = 3;
         }
         else if (i==3) {
          	a = 3;
            c = 2;
         }
         else if (i==4) {
          	a = 3;
            c = 4;
         }
         else if (i==5) {
          	a = 4;
            c = 1;
         }
         else if (i==6) {
          	a = 4;
            c = 3;
         }
         if (j==1) {
         	b = 1;
            d = 4;
         }
         else if (j==2) {
          	b = 2;
            d = 3;
         }
         else if (j==3) {
          	b = 3;
            d = 2;
         }
         else if (j==4) {
          	b = 3;
            d = 4;
         }
         else if (j==5) {
          	b = 4;
            d = 1;
         }
         else if (j==6) {
          	b = 4;
            d = 3;
         }

      	//output the header
         out << "\n                                       Y\n";
         out << " ------------------------------------------------------------------------------\n";
         out << "   A    A    A    A    C    C    C    C    G    G    G    G    T    T    T    T  \n";
         out << "   A    C    G    T    A    C    G    T    A    C    G    T    A    C    G    T  \n";
         out << " ------------------------------------------------------------------------------\n";
         out << "                                  5' ------> 3'\n";
         out << "                                   ";
         if (i==1) out << "A";
         else if (i==2) out << "C";
         else if (i==3||i==4) out << "G";
         else out << "T";
         out << " \\/ \\_/ ";
         if (j==1) out << "A\n";
         else if (j==2) out << "C\n";
         else if (j==3||j==4) out << "G\n";
         else out << "T\n";          
         out << "                                   ";
         if (i==5) out << "A";
         else if (i==3) out << "C";
         else if (i==2||i==6) out << "G";
         else out << "T";
         out << " /\\  |  ";
         if (j==5) out << "A\n";
         else if (j==3) out << "C\n";
         else if (j==2||j==6) out << "G\n";
         else out << "T\n";
         out << "                                  3' <------ 5'\n";
       	for (n=1;n<=4;n++) {
          	for (m=1;m<=4;m++) {
            	for (l=1;l<=4;l++) {
               	for (k=1;k<=4;k++) {
                     //first check for a w-c pair in the loop:

					np=n;
					lp=l;
					mp=m;
					kp=k;
                     //if (inc[n][m]||inc[l][k]) {
						// if (inc[n][m]) {
						//	if (n==3) np = 1;
						//	if (n==4) np = 2;
						//	if (m==3) mp = 1;
						//	if (m==4) mp=2;
						// }
						// if (inc[l][k]) {
						//	if (l==3) lp=1;
						//	if (l==4)  lp=2;
						//	if (k==3) kp=1;
						//	if (k==4) kp=2;

						// }
						 //check for symmetric loop with caonical pairs
						 //if (inc[n][m]&&inc[l][k]) {
						//	 if (n==k&&m==l&&a==d&&b==c&&a==b&&c==d) {


						//		iloop22[a][b][c][d][n][l][m][k]=iloop22[a][b][c][d][np][lp][mp][kp];

						//	 }						

						// }
					 //}



                  	//calculate the value of iloop[...] if it is infinity now
                    /* if (iloop22[a][b][c][d][n][l][m][k]==infinity) {


						iloop22[a][b][c][d][n][l][m][k]=(iloop22[a][c][c][a][np][mp][mp][np]
                              + iloop22[d][b][b][d][kp][lp][lp][kp])/2;




						//check for deltas:

						if (np==3&&mp==3) {

							if (!(lp==4&&kp==4)&&!(lp==1&&kp==3)&&!(lp==3&&kp==1)&&!(lp==3&&kp==3)) {
								
								iloop22[a][b][c][d][n][l][m][k]+=gg;	

							}


						}
						else if (lp==3&&kp==3) {


							if (!(np==4&&mp==4)&&!(np==1&&mp==3)&&!(np==3&&mp==1)) {
								
								iloop22[a][b][c][d][n][l][m][k]+=gg;	

							}


						}


						if ((np==4&&mp==4)) {
							if ((lp==1&&kp==1)) iloop22[a][b][c][d][n][l][m][k]+=delta;

						}
							
							
						else if ((np==3&&mp==1)||(np==1&&mp==3)) {
							if ((lp==2&&kp==2)||(lp==2&&kp==4)||(lp==4&&kp==2))

								iloop22[a][b][c][d][n][l][m][k]+=delta;
							

						}

						else if ((kp==4&&lp==4)) {
							if ((np==1&&mp==1)) iloop22[a][b][c][d][n][l][m][k]+=delta;

						}
							
							
						else if ((kp==3&&lp==1)||(kp==1&&lp==3)) {
							if ((np==2&&mp==2)||(np==2&&mp==4)||(np==4&&mp==2))

								iloop22[a][b][c][d][n][l][m][k]+=delta;
							

						}

						



					 }*/

					iloop22[a][b][c][d][n][l][m][k]=0+tstack[a][c][n][m]+tstack[d][b][k][l];
					if (a==4||c==4) iloop22[a][b][c][d][n][l][m][k]+=32;
					if (b==4||d==4) iloop22[a][b][c][d][n][l][m][k]+=32;
                  	if (iloop22[a][b][c][d][n][l][m][k]>500) {
                     	strcpy(temp,".");
                        out << "   . ";
                     }
                     else {
                     	if (iloop22[a][b][c][d][n][l][m][k]>=0) {
                         	out << "  ";
							writeinteger = iloop22[a][b][c][d][n][l][m][k];

                        }
                        else {
                         	out << " -";
                           writeinteger =
                           	-iloop22[a][b][c][d][n][l][m][k];
                        }
                        tenth = writeinteger%10;
                        first = (writeinteger-tenth)/10;
                        itoa(first,temp,10);
                        out << temp << ".";
                        itoa(tenth,temp,10);
                        out << temp;
                     }
                     


                  }
               }
               out << "\n";
            }
           
         }
      }
   }



   out.close();

}

void errmsg(int i,int j) {

 //here is an error message thing for algorithm.cpp

}

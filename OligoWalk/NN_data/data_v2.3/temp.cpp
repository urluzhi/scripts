
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

//open the files using the C++ method for reading
 ml1.open(miscloop);
 lo1.open(loop2);
 st1.open(stackf);
 th1.open(tstackh);
 ti1.open(tstacki);
 tl1.open(tloop);
 da1.open(danglef);
 in1.open(int22);
 in2.open(int21);
 tri.open(triloop);
 co1.open(coax);
 co2.open(tstackcoax);
 co3.open(coaxstack);
 st2.open(tstack);
 tsm.open(tstackm);
 i11.open(int11);
 hl1.open(hexaloop);
 t23.open(tstacki23);
 t1n.open(tstacki1n);

/* Read information from miscloop */

// the key sequence "-->" now indicates a record


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


 /*	read info from dangle */
 //add to dangle the case where X (represented as 0) is looked up
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


/*	read info from loop for internal loops, hairpin loops, and bulge loops */

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


/* Read info from stack */
//add to the stack table the case where X (represented as 0) is looked up:

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





/* Read info from tstackh */
//add to the tstackh table the case where X (represented as 0) is looked up:

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

/* Read info from tstacki */
//add to the tstacki table the case where X (represented as 0) is looked up:

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

/* Read info from tstacki23 */
//add to the tstacki table the case where X (represented as 0) is looked up:

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


//Read the 2x2 internal loops
//key iloop22[a][b][c][d][j][l][k][m] =
//a j l b
//c k m d


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

//Read the 2x1 internal loop data
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
/*	Read info from triloops */
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


/*	Read info from hexaloops */
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


/* Read info from coax */
//add to the stack table the case where X (represented as 0) is looked up:

//this is the array that keeps track of flush coaxial stacking (no intervening nucs)

//data arrangement of coax: data->coax[a][b][c][d]
//5'b-c3'
//3'a d5'
//this means the helix backbone is continuous between nucs b and c


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

/* Read info from tstackcoax */
//add to the tstackh table the case where X (represented as 0) is looked up:

//data arrangement of tstackcoax:
//5'a-c -> strand continues into stack
//3'b-d -> strand does not continue to stack
//pair between a-b, c-d is a mismatch


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
/* Read info from coaxstack */
//add to the tstackh table the case where X (represented as 0) is looked up:

//data arrangement of coaxstack:
//5'a-c ->strand contnues into stack
//3'b d ->strand does not continue to stack
//pair between a-b, mismatch between c-d
//backbone is discontinuous between b and d


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



/* Read info from tstack */
//this is the terminal mismatch data used in intermolecular folding
//add to the tstack table the case where X (represented as 0) is looked up.
//also add the case where 5 (the intermolecular linker) is looked up,
//this is actually a dangling end, not a terminal mismatch.

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


/* Read info from tstackm */
//add to the tstackm table the case where X (represented as 0) is looked up:

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

//data arrangement for 1x1 loop tables iloop11[a][b][c][d][e][f]:
//abc
//def


//Read the 1x1 internal loop data
//encode the data like:  abc
//                       def where b-e is a mismatch
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



 ml1.close();
 lo1.close();
 st1.close();
 th1.close();
 ti1.close();
 tl1.close();
 da1.close();
 in1.close();
 in2.close();
 tri.close();
 co1.close();
 co2.close();
 co3.close();
 st2.close();
 tsm.close();
 i11.close();
 hl1.close();
return 1;


}

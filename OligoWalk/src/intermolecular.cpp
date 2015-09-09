/*=======================================================================
intermolecular.h and intermolecular.cpp inculde funcions calculating and report
different free energy for binding in OligoWalk.
intermolecular_test.cpp fold the whole sequence, save partition function in a file 
and reuse it.
stochastic sampling method were also embedded.

They are revised based on Mathews' code from RNAStructure.
olig() generate the energy data;
report save the generated data;
siprefileter and filterbysirna() were used to prefilter and postfileter the 
functional siRNA respectively, using criterias other than thermodynamics


															----Aug. 6, 2006
															John(Zhi Lu)

=======================================================================*/
//#define debugmode true
#undef debugmode
#define TESTMODE false   //read the save file instead of calculating if in TEST mode
#include "intermolecular.h"


//define the parameters for alltrace()
#define percent 100	//100% of the suboptimal structure windows
#define delta 6	//only the suboptimal structures having energy less than 2.8kcal/mol higher than optimal one
                    //not store structures in ct file for alltrace() to save memory
// #define MaxStructure 10000000 (10M) is defined in alltrace.h



void dgetdat(char *loop, char *stackf, char *tstackh, char *tstacki,
		char *tloop, char *miscloop, char *danglef, char *int22,
      char *int21,char *coax, char *tstackcoax,
      char *coaxstack, char *tstack, char *tstackm, char *triloop,
      char *int11, char *path);

int readrd (rddata* data,char* dnarna);

inline void scancopy(OligoPclass *region, OligoPclass *copyregion);
inline void scancopyend(OligoPclass *region, OligoPclass *copyregion) ;


/*=======================================================================
oligo fills the array table with thermodynamic data for each oligo (2nd dimension)
in the second dimension:
table[i][0] = overall DG
table[i][1] = duplex DG
table[i][2] = free energy of breaking target structure
table[i][3] = intramolecular oligo free energy
table[i][4] = intermolecular oligo free energy
table[i][5] = Tm of duplex (x10)

option 1 - break local target structure to bind oligo 
option 2 - refold target RNA after oligo binding
option 3 - no target structure considered

Usesub 0 - only consider optimal structure
Usesub 1 - like choice 3,using suboptimal structures,but the whole set from alltarce() function prediction
Usesub 2 - using partition funcion considering every possible structure  of target
			   suboptimal 2 can only used with mold 2
Usesub 3 - using suboptimal structures (heuristic method) for both oligo-free and oligo-bound target RNA
Usesub 4 - using stochasit sampling method to sample 1000 structures

prefileter 1 - using creteria to prefill functional siRNA
foldsize >0  - only folding a fragment with size=foldsize+binding length, 
			   which is centered on the siRNA binding region
			   when foldsize>1, only option 2 plus Usesub 2 is the availabe option
=======================================================================*/
void olig(int isdna, int option, structure *ct, int length,double c, int **table,int **numofsubstructures,
		  datatable& data,datatable& ddata, rddata *hybriddata, int Usesub,TProgressDialog *update,
		  thermo* helixstack,int start, int stop, siPREFILTER *prefilter,int foldsize,int distance,char *shapefile,int *TEST,bool WRITE) {
	
	int i,j,k;
	int ip,jp;
	int foldstart, foldstop;
	int dh,ds,dgeff;
	int *energyarray, *ctenergy;//ctenergy is the energy array to store structure energies from alltrace()
	int numofstructures;
	long double k1b,k1u;
	long double fn,sn,sum;
	char savefile[50],pfsfile[50],pos[50];
	int energy=infinity;//to store the free energy of lowest free energy structure
	int *temp,**temp2;//temp will store basepairing info in option == 1
	

	PFPRECISION Q,Qc,Qolig,pftemp=310.15;//store partion function without and with constrains
	PFPRECISION rescaleinrefill;//rescaling facor for Qc relative to Q
	pfdatatable *pfdata,*dpfdata;  //store the data tables for partionfunction
	OligoPclass *interoligo,*intraoligo;
	OligoPclass *target,*targetcopy,*targettemp;
	structure oligo,oligo1;//to store the structural info for the oligo
						   //oligo is intermolecular structure; oligo1 is intramolecular
	structure fracct;// to store folded fraction of target centered at siRNA binding site
	
	//define shape information for the whole sequnence (ct)
	if ( shapefile[0] != '\0'  ) {
		ct->SHAPEslope = 32;
		ct->SHAPEintercept = -10;
		ct->ReadSHAPE(shapefile) ;
		if (Usesub==2) pfshape(ct,pftemp);
	}


	
	//----------------------------------------------------------------------------------------------
	//allocate space in oligo for the sequence
	//difine for inter oligo sequence, oligo structure can be used both by intra and inter molecular
	//oligo:  intermolecular structure   aagucXXXggcaa
	//oligo1: intramolecular structure   aaguc
	oligo.allocate(2*length+4);
	oligo1.allocate(length+1);
	strcpy(oligo.ctlabel[1],"Oligo_inter");
	strcpy(oligo1.ctlabel[1],"Oligo_intra");
	oligo.numofbases = 2*length+3;
	oligo1.numofbases = length;
	for (j=1;j<=3;j++) {		oligo.inter[j-1] = length + j;	}

	foldsize= foldsize/2*2; // a trick to make foldsize even
	//----------------------------------------------------------------------------------------------
	//define the base pair distance constrain
	if (distance > 0) {
		ct->limitdistance =true;
		ct->maxdistance = distance;
	}
	//define the size of fracct, which is the region to be folded
	if (foldsize >0) {
		fracct.allocate(foldsize+length+1);
		strcpy(fracct.ctlabel[1],"fraction_of_target");
		fracct.numofbases=foldsize+length;
		//define the base pair distance constrain
		if (distance > 0) {
			fracct.limitdistance =true;
			fracct.maxdistance = distance;
		}
	}	
	
	//----------------------------------------------------------------------------------------------
	//some variables declared for each option 
	if (option ==1) {
		if (Usesub==0)	{
			temp = new int [length];
			ct->numofstructures = 1;//only interested in lowest free energy structure
			if (foldsize == 0) {	
				ct->nnopair=0;
				dynamic(ct,&data,1,10,0,0);
				efn2(&data,ct);
		   		energy = ct->energy[1];
			}
		}
		else if (Usesub==3) {
			if (foldsize == 0) {	
				ct->nnopair=0;
				dynamic(ct,&data,1000,10,0);
				efn2(&data,ct);
				energyarray = new int [ct->numofstructures+1];
				for (k=1;k<=ct->numofstructures;k++) energyarray[k] = ct->energy[k];
				temp2 = new int *[ct->numofstructures+1];
				for (i=1;i<=ct->numofstructures;i++)   		temp2[i] = new int [ct->numofbases];
				numofstructures=ct->numofstructures;
			}
      	}
		
	}
	else if (option==2) {
   		temp = new int [ct->numofbases+1];
		for (j=1;j<=ct->numofbases;j++) {
      		temp[j] = ct->basepr[1][j];
		}
		if (Usesub ==0 ) {
			ct->numofstructures = 1;//only interested in lowest free energy structure
			if (foldsize == 0) {	
				ct->nnopair=0;
				strcpy(savefile, ct->ctlabel[1]);
				savefile[strlen(savefile)-1]='\0'; //get off the new line charactor in the string
				strcat(savefile,"_s0.sav");
				//check if a sav file of partition function result exist with a C i/o function
				ifstream sav(savefile,ios::binary);
				if(TESTMODE && sav) {
					//read information from file if it was found
					read(&sav,&energy);
					sav.close();
				}
				else {
					//close the readonly file sav
					sav.close();
					//write the save file information so that the fold need not to be done again
					dynamic(ct,&data,1000,10,0,0,true);
				
					energy = ct->energy[1];
					if (WRITE) {
 						ofstream sav(savefile,ios::binary);
						write(&sav,&energy);
						sav.close();
					}
				}
			}
		}
		else if	(Usesub==1)	 {	
			ctenergy=new int [MaxStructure+1];
			if (foldsize ==0) {	
				ct->nnopair=0;
				strcpy(savefile, ct->ctlabel[1]);
				savefile[strlen(savefile)-1]='\0'; //get off the new line charactor in the string
				strcat(savefile,"_s1.sav");
				//check if a sav file of partition function result exist with a C i/o function
				ifstream sav(savefile,ios::binary);
				if(TESTMODE && sav) {
					//read information from file if it was found
					read(&sav,&numofstructures);
					energyarray = new int [numofstructures+1];
					for (j=1;j<=numofstructures;j++)	read(&sav, energyarray+j);
					sav.close();
				}
				else {
					//close the readonly file sav
					sav.close();
					//calculate the whole set of suboptimal structures
					alltrace(false,ctenergy,ct,&data, percent, delta, NULL,NULL);
					//write the save file information so that the fold need not to be done again
					numofstructures=ct->numofstructures;
					energyarray = new int [numofstructures+1];
					for (j=1;j<=numofstructures;j++)	{energyarray[j]=ctenergy[j];}
					if (WRITE) {
 						ofstream sav(savefile,ios::binary);
						write(&sav,&numofstructures);
						for (j=1;j<=numofstructures;j++)		write(&sav, &energyarray[j]);
						sav.close();
					}
				}

			}
		}
		else if	(Usesub==3)	 {	
			if (foldsize ==0) {	
				ct->nnopair=0;
				strcpy(savefile, ct->ctlabel[1]);
				savefile[strlen(savefile)-1]='\0'; 
				strcat(savefile,"_s3.sav");
				//check if a sav file of partition function result exist with a C i/o function
				ifstream sav(savefile,ios::binary);
				if(TESTMODE && sav ) {
					//read information from file if it was found
					read(&sav,&numofstructures);
					energyarray = new int [numofstructures+1];
					for (j=1;j<=numofstructures;j++)	read(&sav, energyarray+j);
					sav.close();
				}
				else {
					//close the readonly file sav
					sav.close();
					dynamic(ct,&data,1000,10,0);
					//write the save file information so that the fold need not to be done again
					numofstructures=ct->numofstructures;
					energyarray = new int [numofstructures+1];
					for (j=1;j<=numofstructures;j++)	{energyarray[j]=ct->energy[j];}
					if (WRITE) {
 						ofstream sav(savefile,ios::binary);
						write(&sav,&numofstructures);
						for (j=1;j<=numofstructures;j++)		write(&sav, &energyarray[j]);
						sav.close();
					}
				}

			}
		}
		else if (Usesub==2) {
			if (isdna !=0 ) {//oligo is a DNA
				dpfdata = new pfdatatable (&ddata,scalingdefinition,pftemp);
				pfdata = new pfdatatable (&data,scalingdefinition,pftemp);
				intraoligo = new OligoPclass(&oligo1,dpfdata);
				interoligo = new OligoPclass(&oligo,dpfdata);
			}
			else {//oligo is a RNA
				pfdata = new pfdatatable (&data,scalingdefinition,pftemp);
				intraoligo = new OligoPclass(&oligo1,pfdata);
				interoligo = new OligoPclass(&oligo,pfdata);
			}
			//calculate partion function for the whole target without constrain
			ct->nnopair=0;
			if (foldsize==0) {//folding the whole sequence at one time
				target = new OligoPclass(ct,pfdata);
				strcpy(savefile, ct->ctlabel[1]);
				savefile[strlen(savefile)-1]='\0'; //get off the new line charactor in the string
				strcat(savefile,"_s2.sav");
				strcpy(pfsfile, ct->ctlabel[1]);
				pfsfile[strlen(pfsfile)-1]='\0'; //get off the new line charactor in the string
				strcat(pfsfile,"_s2.pfs");
				//check if a sav file of partition function result exist with a C i/o function
				ifstream sav(savefile,ios::binary);
				ifstream pfs(savefile,ios::binary);
				if (TESTMODE && sav && pfs  ){
					pfs.close();
					//read information from file if it was found
					read(&sav,&Q);
					for (i=0;i<=ct->numofbases;i++) {
						read(&sav,&(target->copyw5[i]));
						for (j=0;j<=ct->numofbases;j++) {
							read(&sav,&(target->copyv->f(i,j)));
							read(&sav,&(target->copyw->f(i,j)));
							read(&sav,&(target->copywca[i][j]));
							read(&sav,&(target->copywmb->f(i,j)));
							read(&sav,&(target->copywl->f(i,j)));
							read(&sav,&(target->copywmbl->f(i,j)));
							read(&sav,&(target->copywcoax->f(i,j)));
						}	
					}
					sav.close();
				}
				else {
					//close the readonly file sav
					sav.close();
					pfs.close();
					//calculate the partition function if no file found
					if (WRITE) {
						target->partition4refill(&Q,pfsfile);
					}
					else target->partition4refill(&Q);

				
					//write the save file information so that the partition function can be re-folded,
					if (WRITE) {
						ofstream sav(savefile,ios::binary);
						write(&sav,&Q);
						for (i=0;i<=ct->numofbases;i++) {
							write(&sav,&(target->copyw5[i]));
							for (j=0;j<=ct->numofbases;j++) {
								write(&sav,&(target->copyv->f(i,j)));
								write(&sav,&(target->copyw->f(i,j)));
								write(&sav,&(target->copywca[i][j]));
								write(&sav,&(target->copywmb->f(i,j)));
								write(&sav,&(target->copywl->f(i,j)));
								write(&sav,&(target->copywmbl->f(i,j)));
								write(&sav,&(target->copywcoax->f(i,j)));
							}		
						}
						sav.close();
					}
						
				}
			}
			else if(prefilter->useit!=0) {//using prefilter and fold different region each time
				target = new OligoPclass(&fracct,pfdata);

			}
		}
		else if (Usesub==4) {
			if (isdna) {//oligo is a DNA
				dpfdata = new pfdatatable (&ddata,scalingdefinition,pftemp);
				pfdata = new pfdatatable (&data,scalingdefinition,pftemp);
				intraoligo = new OligoPclass(&oligo1,dpfdata);
				interoligo = new OligoPclass(&oligo,dpfdata);
			}
			else {//oligo is a RNA
				pfdata = new pfdatatable (&data,scalingdefinition,pftemp);
				intraoligo = new OligoPclass(&oligo1,pfdata);
				interoligo = new OligoPclass(&oligo,pfdata);
			}
			//calculate partion function for the whole target without constrain
			ct->nnopair=0;
			if (foldsize==0) {//folding the whole sequence at one time
				target = new OligoPclass(ct,pfdata);
				strcpy(pfsfile, ct->ctlabel[1]);
				pfsfile[strlen(pfsfile)-1]='\0'; //get off the new line charactor in the string
				strcat(pfsfile,"_s2.pfs");
				//check if a sav file of partition function result exist with a C i/o function
				ifstream pfs(pfsfile,ios::binary);
				if(TESTMODE && pfs) {//read information from file if it was found
					pfs.close();				
				}
				else{
					//close the readonly file sav
					pfs.close();
					//calculate the partition function if no file found
					target->partition4refill(&Q,pfsfile); 
					
				}
				stochastic(ct,pfsfile,1000,5);
				efn2(&data,ct);
				//ofstream sout("RL_sample_oligo.out");
				//for (j=1;j<=1000;j++)	{ sout << ct->energy[j]<<"\n";	}
				//sout.close();
				numofstructures=ct->numofstructures;
				energyarray = new int [numofstructures+1];
				for (j=1;j<=numofstructures;j++)	{energyarray[j]=ct->energy[j];}
			}
			else if(prefilter->useit!=0) {//using prefilter and fold different region each time
				//allocate target only, if prefilter==0, allocate both target and targetcopy later
				target = new OligoPclass(&fracct,pfdata);

			}
		}
		
	}
	//having single-strand constrained fro the calculation of constrained energy
	ct->nnopair=length;
	ct->checknopair();	

	//------------------------------------------------------------------------------------------
	//------------------------------------------------------------------------------------------
	//-------------------------------------------------------------------------------------------
	//begin to scan the target from start to stop
	for (i=start; i<=stop; i++) {

		if (update!=0) {
   			update->update (int((double (i-start))*100/(double (stop-start))));
		}
		//define the oligo sequence
		oligo1.intermolecular = false;
		oligo1.numofbases = length;
		//set sequence
		 for (j=1;j<=length;j++) {
			oligo1.numseq[j] = complement(i+length-j,ct);
			oligo.numseq[j] = complement(i+length-j,ct);
		}
   		//-------------------------------------------------------------------------------------------
		//prefiltering the functional siRNA
		if (prefilter->useit != 0) {		
			prefilter->count(&oligo1,i,TEST[i]);
			if (prefilter->score[i] < FILTER_PASS) {
				continue;
			}
		}
	
	//-------------------------------------------------------------------------------------------
	//-------------------------------------------------------------------------------------------
	//calculate the free energy of breaking target structure

	//Option 1: break local structure only
    
	if (option==1) {
		if (foldsize ==0) {
			//-------------------------------------------------------------------------------------------
			//option 1 + consider only the first structure:
			   if (Usesub ==0){
				//store basepairing info
		     	for (j=0;j<length;j++) {
					temp[j] = ct->basepr[1][i+j];
		        	if (ct->basepr[1][i+j] > 0) {
		        		ct->basepr[1][ct->basepr[1][i+j]] = 0;
					   	ct->basepr[1][i+j] = 0;
			       	}
		     	}
         	
				efn2(&data,ct);
		    	table[i][2] = ct->energy[1] - energy;
         	
				//restore basepairing
		    	for (j=0;j<length;j++) {
		      		if (temp[j]>0) {
		         		ct->basepr[1][i+j] = temp[j];
		           		ct->basepr[1][temp[j]] = i+j;
		        	}
		    	}
			}	
			//-------------------------------------------------------------------------------------------
			//option 1 + consider all subotimal structures:
			else if (Usesub==3) { 
				sn = 0;
		        sum = 0;
	          //store the basepairing info and force the region of hybridization
	         //to have no structure
	          for (k=1;k<=ct->numofstructures;k++) {
		          	for (j=0;j<length;j++) {//store basepairing info
		     			temp2[k][j] = ct->basepr[k][i+j];
		           		if (ct->basepr[k][i+j] != 0) {
		         			ct->basepr[k][ct->basepr[k][i+j]] = 0;
	             			ct->basepr[k][i+j] = 0;
		          		}
		      		}

				}
				efn2(&data,ct);
		      for (k=1;k<=ct->numofstructures;k++) {

		          	fn = expl(-(((long double)energyarray[k]-energyarray[1])/(rt*conversionfactor)));
		           sn = sn + fn;
		   			sum = sum+fn*(((long double)(ct->energy[k]- energyarray[k])));
		       }
		        //restore the basepairing
		        for (k=1;k<=ct->numofstructures;k++) {
			       	for (j=0;j<length;j++) {
						if (temp2[k][j]>0) {
           					ct->basepr[k][i+j] = temp2[k][j];
               				ct->basepr[k][temp2[k][j]] = i+j;
           				}
        			}
				}
         		
	            table[i][2] = (int)(sum/sn);
				numofsubstructures[i][0]=numofstructures; //report in report() function
				numofsubstructures[i][1]=ct->numofstructures; 

			}
		}
		if (foldsize > 0 ) {
			//define the region for refolding 
			if ( (i-foldsize/2) > 1 && (i+length-1+foldsize/2) < ct->numofbases ) {
				for (j=1;j<=foldsize+length;j++) {
      				fracct.numseq[j] = ct->numseq[i-foldsize/2+j-1];
				}
			}
			else if ( (i-foldsize/2)<=1 ) 
			for (j=1;j<=foldsize+length;j++) {
      			fracct.numseq[j] = ct->numseq[j];
			}
			else if( (i-1+length+foldsize/2)>=ct->numofbases )
			for (j=1;j<=foldsize+length;j++) {
      			fracct.numseq[j] = ct->numseq[(ct->numofbases)+j-foldsize-length];
			}
			fracct.nnopair=0;
			//calculate the strucutre without constrains
			if (Usesub==0)	{
				fracct.numofstructures = 1;//only interested in lowest free energy structure
				dynamic(&fracct,&data,1,10,0,0);
				efn2(&data,&fracct);
		   		energy = fracct.energy[1];
			}
			else if (Usesub==3) {
				dynamic(&fracct,&data,1000,10,0);
				efn2(&data,&fracct);
				energyarray = new int [fracct.numofstructures+1];
				for (k=1;k<=fracct.numofstructures;k++) energyarray[k] = fracct.energy[k];
				temp2 = new int *[fracct.numofstructures+1];
				for (k=1;k<=fracct.numofstructures;k++)   		temp2[k] = new int [fracct.numofbases];
				numofstructures=fracct.numofstructures;
			
      		}
			//set the constrained position
			if (i-foldsize/2<=1) {//constrained positions are different at two ends, as the folded region did not move
				foldstart=i;
				foldstop=i+length-1;
			} 
			else if(i+length-1+foldsize/2>= (ct->numofbases) ) {
				foldstart=foldsize+length-(ct->numofbases) +i;
				foldstop=foldsize+length+length-1-(ct->numofbases) +i;
			}
			else {//folded region begin to move, constrained position will be always in the middle of this region
				foldstart=foldsize/2+1;
				foldstop=foldsize/2+length;
			}
			
			//recalc the energy again
			if (Usesub==0) {
				//store basepairing info
		     	for (j=foldstart;j<=foldstop;j++) {
					temp[j-foldstart] = fracct.basepr[1][j];
		        	if (fracct.basepr[1][j] > 0) {
		        		fracct.basepr[1][fracct.basepr[1][j]] = 0;
					   	fracct.basepr[1][j] = 0;
			       	}
		     	}
         	
				efn2(&data,&fracct);
		    	table[i][2] = fracct.energy[1] - energy;
         	
				//restore basepairing
		    	for (j=foldstart;j<=foldstop;j++) {
		      		if (temp[j-foldstart]>0) {
		         		fracct.basepr[1][j] = temp[j-foldstart];
		           		fracct.basepr[1][temp[j-foldstart]] = j;
		        	}
		    	}
			}
			else if (Usesub==3) {
				sn = 0;
		        sum = 0;
	            //store the basepairing info and force the region of hybridization
	            //to have no structure
				for (k=1;k<=fracct.numofstructures;k++) {
		          	for (j=foldstart;j<=foldstop;j++) {//store basepairing info
		     			temp2[k][j] = fracct.basepr[k][j];
		           		if (fracct.basepr[k][j] != 0) {
		         			fracct.basepr[k][fracct.basepr[k][j]] = 0;
	             			fracct.basepr[k][j] = 0;
		          		}
		      		}

				}
				efn2(&data,&fracct);
				for (k=1;k<=fracct.numofstructures;k++) {

		          	fn = expl(-(((long double)energyarray[k]-energyarray[1])/(rt*conversionfactor)));
		            sn = sn + fn;
		   			sum = sum+fn*(((long double)(fracct.energy[k]- energyarray[k])));
				}
		        //restore the basepairing
		        for (k=1;k<=fracct.numofstructures;k++) {
			       	for (j=foldstart;j<=foldstop;j++) {
						if (temp2[k][j]>0) {
           					fracct.basepr[k][j] = temp2[k][j];
               				fracct.basepr[k][temp2[k][j]] = j;
           				}
        			}
				}
			
			    table[i][2] = (int)(sum/sn);
				numofsubstructures[i][0]=numofstructures; //report in report() function
				numofsubstructures[i][1]=fracct.numofstructures; 

				delete[] energyarray;
				for (k=1;k<=fracct.numofstructures;k++) {	delete[] temp2[k];	}
				delete[] temp2;

			}
		}
	}


	//-------------------------------------------------------------------------------------------
	//option 2: refold  rna
	else if (option == 2) {
		//not scan: refolding the whole sequence
		if (foldsize ==0) {
			//reset the nopair constrains
			for (j=0;j<length;j++)		{	ct->nopair[j+1] = i+j;	}		
			
			//option 2 +notrefold+ consider only the first suboptimal structure
			if (Usesub ==0) {
				dynamic(ct,&data,1000,10,0,0,true);
			
				table[i][2] = ct->energy[1] - energy;
			}
			//option 2 +notrefold+ consider all suboptimal structures with ensemble energy
			else if (Usesub==1) {
				//sum of unconstrained energy from energyarray[]
				Q = 0;
				for (k=1;k<=numofstructures;k++) {
	          		fn = expl(-(((long double)energyarray[k]-energyarray[1])/(rt*conversionfactor)));
					Q = Q + fn;
         		}
								
				//calculate the whole set of suboptimal structures with constrain
				alltrace(false,ctenergy,ct,&data, percent, delta, NULL,NULL);
				//sum of constrained energy now
				Qc = 0;
				for (k=1;k<=ct->numofstructures;k++) {
	          		fn = expl(-(((long double)ctenergy[k]-energyarray[1])/(rt*conversionfactor)));
					Qc = Qc + fn;
         		}
				//final free energy difference					
				table[i][2] = (int)( conversionfactor*rt*logl( (long double)Q/(long double)Qc ) );

				numofsubstructures[i][0]=numofstructures; //report in report() function
				numofsubstructures[i][1]=ct->numofstructures; 
			}
			//option 2 +notrefold+ consider heuristic suboptimal structures with average energy
			else if (Usesub==3) {
				sum = 0;
				sn = 0;
				for (k=1;k<=numofstructures;k++) {
	          		fn = expl(-(((long double)energyarray[k]-energyarray[1])/(rt*conversionfactor)));
					sn = sn + fn;
         			sum = sum + fn*((long double)energyarray[k]);
				}
				table[i][2] = (int)(sum/sn);
							
				//save the energies in a save file
				strcpy(savefile, ct->ctlabel[1]);
				savefile[strlen(savefile)-1]='\0'; 
				sprintf(pos,"_%d",i);
				strcat(savefile,pos);
				strcat(savefile,"_s3_0_constrain.sav");
				//fold the constrained structure
				dynamic(ct,&data,1000,10,0);
				if (WRITE) {
					ofstream sav(savefile,ios::binary);
					write(&sav,&(ct->numofstructures));
					for (j=1;j<=ct->numofstructures;j++)		write(&sav, & (ct->energy[j]));
					sav.close();
				}
				sum = 0;
				sn = 0;
				for (k=1;k<=ct->numofstructures;k++) {
	          		fn = expl(-(((long double)ct->energy[k]-ct->energy[1])/(rt*conversionfactor)));
					sn = sn + fn;
         			sum = sum + fn*((long double)ct->energy[k]);
					//cout<< ct->energy[k]<<"\n";
				}
				table[i][2] = (int)(sum/sn) - table[i][2] ;
				numofsubstructures[i][0]=numofstructures; //report in report() function
				numofsubstructures[i][1]=ct->numofstructures; 
		   	}
			//option 2 +notrefold+ partionfunction
			else if (Usesub ==2) {
				//save the energies in a save file
				strcpy(pfsfile, ct->ctlabel[1]);
				pfsfile[strlen(pfsfile)-1]='\0'; 
				sprintf(pos,"_%d",i);
				strcat(pfsfile,pos);
				strcat(pfsfile,"_s2_0_constrain.pfs");
				if (WRITE) {
					target->refill(ct,&Qc,i, i+length-1,rescaleinrefill,pfsfile);
				}
				else target->refill(ct,&Qc,i, i+length-1,rescaleinrefill);


				table[i][2] = (int)(conversionfactor*rt* 
					( (long double)(ct->numofbases)*logl(rescaleinrefill) 
					+ logl((long double)Q/(long double)Qc) ) );
				
			}
			//option 2 +notrefold+ stochastic
			else if (Usesub ==4) {
				sum = 0;
				for (k=1;k<=numofstructures;k++) {
	          		sum = sum + energyarray[k];
				}
				table[i][2] = (int)(sum/numofstructures);
			
				strcpy(pfsfile, ct->ctlabel[1]);
				pfsfile[strlen(pfsfile)-1]='\0'; 
				sprintf(pos,"_%d",i);
				strcat(pfsfile,pos);
				strcat(pfsfile,"_s2_0_constrain.pfs");
				ifstream pfs(pfsfile,ios::binary);
				if(TESTMODE && pfs) {
					pfs.close();
				}
				else{
					pfs.close();
					target->refill(ct,&Qc,i, i+length-1,rescaleinrefill,pfsfile);
					
				}
				stochastic(ct,pfsfile,1000,5);
				efn2(&data,ct);
				sum = 0;
				for (k=1;k<=ct->numofstructures;k++) {
	          		sum = sum + ct->energy[k];
				}
				table[i][2] = (int)(sum/ct->numofstructures) - table[i][2] ;
				numofsubstructures[i][0]=numofstructures; //report in report() function
				numofsubstructures[i][1]=ct->numofstructures; 			
			}
			
		}
		//----------------------------------------------------------
		//----------------------------------------------------------
		//scan: only folding the region close to the binding site
		else if (foldsize >0) {
			
			//define the region for refolding 
			if ( (i-foldsize/2) > 1 && (i+length-1+foldsize/2) < ct->numofbases ) {
				for (j=1;j<=foldsize+length;j++) {
      				fracct.numseq[j] = ct->numseq[i-foldsize/2+j-1];
				}
			}
			else if ( (i-foldsize/2)<=1 ) 
			for (j=1;j<=foldsize+length;j++) {
      			fracct.numseq[j] = ct->numseq[j];
			}
			else if( (i-1+length+foldsize/2)>=ct->numofbases )
			for (j=1;j<=foldsize+length;j++) {
      			fracct.numseq[j] = ct->numseq[(ct->numofbases)+j-foldsize-length];
			}
			fracct.nnopair=0;

			//----------------------------------------------------------
			//fold the scanned region without any constrain:
			//refold for different Usesub options
			if (Usesub ==0) {
				dynamic(&fracct,&data,1000,10,0,0,true);
				fracct.numofstructures = 1;//only interested in lowest free energy structure
   				energy = fracct.energy[1];
				
			}
			else if(Usesub==1){
				strcpy(savefile, ct->ctlabel[1]);
				savefile[strlen(savefile)-1]='\0'; 
				sprintf(pos,"_%d_s1_%d.sav",i,foldsize+length);
				strcat(savefile,pos);
				ifstream sav(savefile,ios::binary);
				if (TESTMODE && sav) {
					//read information from file if it was found
					read(&sav,&numofstructures);
					energyarray = new int [numofstructures+1];
					for (j=1;j<=numofstructures;j++)	read(&sav, energyarray+j);
					sav.close();
				}
				else {
					sav.close();
					alltrace(false,ctenergy,&fracct,&data, percent, delta, NULL,NULL);
					numofstructures=fracct.numofstructures;
					energyarray = new int [numofstructures+1];
					for (j=1;j<=numofstructures;j++)		energyarray[j]=ctenergy[j];
				}

				
   			}
			else if(Usesub==3){
				dynamic(&fracct,&data,1000,10,0);
				//save the energies in a save file
				strcpy(savefile, ct->ctlabel[1]);
				savefile[strlen(savefile)-1]='\0'; 
				sprintf(pos,"_%d_s3_%d.sav",i,foldsize+length);
				strcat(savefile,pos);
				if (WRITE) {
					ofstream sav(savefile,ios::binary);
					write(&sav,&(fracct.numofstructures));
					for (j=1;j<=fracct.numofstructures;j++)		write(&sav, & (fracct.energy[j]));
					sav.close();
				}
				numofstructures=fracct.numofstructures;
				energyarray = new int [numofstructures+1];
				for (j=1;j<=numofstructures;j++)		energyarray[j]=fracct.energy[j];
   			}
			else if(Usesub==2) {
				strcpy(pfsfile, ct->ctlabel[1]);
				pfsfile[strlen(pfsfile)-1]='\0'; 
				sprintf(pos,"_%d_s2_%d.pfs",i,foldsize+length);
				strcat(pfsfile,pos);
		
				//not using prefilter, so arrays can be reused when region move to the right
				//fold the first region, the folded region begin to move to right in the middle of target
				if (prefilter->useit == 0) {
					if (i==start) {
						target=new OligoPclass(&fracct,pfdata);
						targetcopy=new OligoPclass(&fracct,pfdata);
						target->partition(true,&Q,NULL);
						//char *report1="report1.out";				
						//target->partition(true,&Q,NULL,report1);
					}
					//when folded region begin moving,reuse some arrays overlapped expect for those on the edges
					else if((i-foldsize/2)>1 && (i+length-1+foldsize/2)<= (ct->numofbases) ){
				   	//char *report2="report2.out";		
					target->scanfill(&fracct,&Q,0);
					}
					//copy the arrays to be reused for next scan region
					//folded region is not moving at two ends, so copy the array outside for next folding without constrain
					if ( (i-foldsize/2)<1 || (i+length-1+foldsize/2)>=(ct->numofbases) ) {
						scancopyend(target,targetcopy);
					}
					//folded region begin moving next, copy the overlapped region outside with different index
					else	scancopy(target,targetcopy);
			
				}
				//using prefilter, not reuse arrays,refold the new region every time
				else {
					
					//refold the new region without using any information from previous folding 
					target->reset4oligo(&fracct);
					if (WRITE) {
						target->partition(true,&Q,NULL,pfsfile);
						
					}
					else target->partition(true,&Q,NULL);
				}
			}
			else if(Usesub==4) {
				strcpy(pfsfile, ct->ctlabel[1]);
				pfsfile[strlen(pfsfile)-1]='\0'; 
				sprintf(pos,"_%d_s2_%d.pfs",i,foldsize+length);
				strcat(pfsfile,pos);
				//not using prefilter, so arrays can be reused when region move to the right
				//fold the first region, the folded region begin to move to right in the middle of target
				if (prefilter->useit == 0) {
					cout <<"Must use prefilter!\n";
					return;
				}
				//using prefilter, not reuse arrays,refold the new region every time
				else {
					ifstream pfs(pfsfile,ios::binary);
					if(TESTMODE && pfs) {
						pfs.close();
					}
					else{
						pfs.close();
						//refold the new region without using any information from previous folding 
						target->reset4oligo(&fracct);
						target->partition(true,&Q,NULL,pfsfile);
						
					}
					stochastic(&fracct,pfsfile,1000,5);
					efn2(&data,&fracct);
				/*	ofstream sout("RL_sample_oligo_scan.out");
					for (j=1;j<=1000;j++)	{ sout << fracct.energy[j]<<"\n";	}
					sout.close();
				*/
					numofstructures=fracct.numofstructures;
					energyarray = new int [numofstructures+1];
					for (j=1;j<=numofstructures;j++)		energyarray[j]=fracct.energy[j];
				}
			}

			
			
			
			//set the single-stranded constrain 
			fracct.nnopair=length;
			fracct.checknopair();
			if (i-foldsize/2<=1) {//constrained positions are different at two ends, as the folded region did not move
				for (j=0;j<length;j++) {
					fracct.nopair[j+1] = i+j;
				}
				foldstart=i;
				foldstop=i+length-1;
			} 
			else if(i+length-1+foldsize/2>= (ct->numofbases) ) {
				for (j=0;j<length;j++) {
					fracct.nopair[j+1] = foldsize+length-(ct->numofbases)+i+j;
				}
				foldstart=foldsize+length-(ct->numofbases) +i;
				foldstop=foldsize+length+length-1-(ct->numofbases) +i;
			}
			else {//folded region begin to move, constrained position will be always in the middle of this region
				for (j=0;j<length;j++) {
					fracct.nopair[j+1] = foldsize/2+1+j;
				}
				foldstart=foldsize/2+1;
				foldstop=foldsize/2+length;
			}

			//---------------------------------------------------------------------
			//refold with constrain:
			//option 2 + refold + consider only the first suboptimal structure
			if (Usesub==0) {
				dynamic(&fracct,&data,1000,10,0,0,true);
				table[i][2] = fracct.energy[1] - energy;			
			}
			//---------------------------------------------------------------------------------------------
			//option 2 + refold + consider all suboptimal structures ensemble energy
			else if (Usesub==1) {
				
				//sum of unconstrained energy from energyarray[]
				Q = 0;
				for (k=1;k<=numofstructures;k++) {
	          		fn = expl(-(((long double)energyarray[k]-energyarray[1])/(rt*conversionfactor)));
					Q = Q + fn;
         		}
				//calculate the whole set of suboptimal structures with constrain
				alltrace(false,ctenergy,&fracct,&data, percent, delta, NULL,NULL);
				//sum of constrained energy now
				Qc = 0;
				for (k=1;k<=fracct.numofstructures;k++) {
	          		fn = expl(-(((long double)ctenergy[k]-energyarray[1])/(rt*conversionfactor)));
					Qc = Qc + fn;
         		}
				//final free energy difference					
				table[i][2] = (int)( conversionfactor*rt*logl( (long double)Q/(long double)Qc ) );

				numofsubstructures[i][0]=numofstructures; //report in report() function
				numofsubstructures[i][1]=fracct.numofstructures; 
								      			
				delete[] energyarray;
		    }
			//---------------------------------------------------------------------------------------------
			//option 2 + refold + consider heuristic suboptimal structures average free energy
			else if (Usesub==3) {
      			sum = 0;
				sn = 0;
				for (k=1;k<=numofstructures;k++) {
	          			fn = expl(-(((long double)energyarray[k]-energyarray[1])/(rt*conversionfactor)));
						sn = sn + fn;
         				sum = sum + fn*((long double)(energyarray[k]));
				}
				table[i][2] = (int)(sum/sn);
				delete[] energyarray;

				//dynamic(&fracct,&data,1000,10,0,0,true);
      			dynamic(&fracct,&data,1000,10,0);
				//save the energies in a save file
				strcpy(savefile, ct->ctlabel[1]);
				savefile[strlen(savefile)-1]='\0'; 
				sprintf(pos,"_%d_s3_%d_constrain.sav",i,foldsize+length);
				strcat(savefile,pos);
				if (WRITE) {
					ofstream sav(savefile,ios::binary);
					write(&sav,&(fracct.numofstructures));
					for (j=1;j<=fracct.numofstructures;j++)		write(&sav, & (fracct.energy[j]));
					sav.close();
				}
				
   			    sum = 0;
				sn = 0;
				for (k=1;k<=fracct.numofstructures;k++) {
	          			fn = expl(-(((long double)fracct.energy[k]-fracct.energy[1])/(rt*conversionfactor)));
						sn = sn + fn;
         				sum = sum + fn*((long double)fracct.energy[k]);
						//cout<< fracct.energy[k]<<"\n";
				}
				table[i][2] = (int)(sum/sn) - table[i][2];
				numofsubstructures[i][0]=numofstructures; //report in report() function
				numofsubstructures[i][1]=fracct.numofstructures; 
		    }
			//---------------------------------------------------------------------------------------------
			//option 2 + refold + partionfunction
			else if (Usesub==2) {
				strcpy(pfsfile, ct->ctlabel[1]);
				pfsfile[strlen(pfsfile)-1]='\0'; 
				sprintf(pos,"_%d_s2_%d_constrain.pfs",i,foldsize+length);
				strcat(pfsfile,pos);
				//reuse the arrays filled without constrained 
				if (WRITE) {
					target->scanconstrain(&fracct,&Qc,foldstart,foldstop,rescaleinrefill,pfsfile);

				}
				else target->scanconstrain(&fracct,&Qc,foldstart,foldstop,rescaleinrefill);
				//exchange arrays to be used for next folding site without constrain
				if (prefilter->useit==0) {//not using targetcopy when prefilter is used
				targettemp=target;
				target=targetcopy;
				targetcopy=targettemp;
				}
			    table[i][2] = (int)(conversionfactor*rt*
				( (long double)(fracct.numofbases)*logl(rescaleinrefill) 
					+ logl((long double)Q/(long double)Qc) ) );
			}
			//---------------------------------------------------------------------------------------------
			//option 2 + refold + stochastic sample
			else if (Usesub==4) {			
				
   			    sum = 0;
				for (k=1;k<=numofstructures;k++) {
	          			sum = sum + energyarray[k];
				}
				table[i][2] = (int)(sum/numofstructures);
				delete[] energyarray;

				strcpy(pfsfile, ct->ctlabel[1]);
				pfsfile[strlen(pfsfile)-1]='\0'; 
				sprintf(pos,"_%d_s2_%d_constrain.pfs",i,foldsize+length);
				strcat(pfsfile,pos);
				ifstream pfs(pfsfile,ios::binary);
				if(TESTMODE && pfs) {
					pfs.close();
				}
				else{
					pfs.close();
					//reuse the arrays filled without constrained 
					target->scanconstrain(&fracct,&Qc,foldstart,foldstop,rescaleinrefill,pfsfile);
					
				}
				stochastic(&fracct,pfsfile,1000,5);
				efn2(&data,&fracct);
				sum = 0;
				for (k=1;k<=fracct.numofstructures;k++) {
	          		sum = sum + fracct.energy[k];
						
				}
				table[i][2] = (int)(sum/fracct.numofstructures) - table[i][2];
				numofsubstructures[i][0]=numofstructures; //report in report() function
				numofsubstructures[i][1]=fracct.numofstructures; 

				//exchange arrays to be used for next folding site without constrain
				if (prefilter->useit==0) {//not using targetcopy when prefilter is used
				/*targettemp=target;
				target=targetcopy;
				targetcopy=targettemp;
				*/
					cout << "must use prefilter for stochasitic folding\n";
					return;
				}
			}
		
		}
	}

	//-----------------------------------------------------------------------------------------------------
	//option 3: no local structure considered

	else {
		table[i][2] = 0;
    }

	
	
	//-----------------------------------------------------------------------------------------------------
	//-----------------------------------------------------------------------------------------------------
	//calculate the stability of the hybrid duplex
		
	//DNA-RNA:
	if (isdna ==1) {
		
		table[i][1] = hybriddata->init;//initiation
        for (j=0;j<(length-1);j++) {
		
			table[i][1] += 
				hybriddata->stack[ct->numseq[i+j]][complement(i+j,ct)][ct->numseq[i+j+1]][complement(i+j+1,ct)];
        }
    }

	//RNA-RNA:
    else if(isdna ==0) {
		table[i][1] = data.init;//initiation
      	for (j=0;j<(length-1);j++) {
          	table[i][1] += 
				data.stack[ct->numseq[i+j]][complement(i+j,ct)][ct->numseq[i+j+1]][complement(i+j+1,ct)];
         }
		//consider AU end effects for RNA/RNA duplexes
		if (ct->numseq[i]==1||ct->numseq[i]==4) table[i][1]+=data.auend;
		if (ct->numseq[i+length-1]==1||ct->numseq[i+length-1]==4) table[i][1]+=data.auend;

    }
	//DNA-DNA
	else if(isdna == 2) {
		table[i][1] = ddata.init;//initiation
      	for (j=0;j<(length-1);j++) {
          	table[i][1] += 
				ddata.stack[ct->numseq[i+j]][complement(i+j,ct)][ct->numseq[i+j+1]][complement(i+j+1,ct)];
         }
	}
	
	//calculate the Tm of the duplex
	if (isdna <=1) { //no dh ds for DNA-DNA yet
	ds = helixstack->dsi;
    dh = helixstack->dhi;
    for (j=0;j<(length-1);j++) {
	  	dh +=helixstack->dh[ct->numseq[i+j]][complement(i+j,ct)][ct->numseq[i+j+1]][complement(i+j+1,ct)];
        ds +=helixstack->ds[ct->numseq[i+j]][complement(i+j,ct)][ct->numseq[i+j+1]][complement(i+j+1,ct)];
   	}
    if (ct->numseq[i]==1||ct->numseq[i]==4) {
		dh = dh + helixstack->dha;
        ds = ds + helixstack->dsa;
    }
    if (ct->numseq[i+length-1]==1||ct->numseq[i+length-1]==4) {
		dh = dh + helixstack->dha;
        ds = ds + helixstack->dsa;
    }
      
	table[i][5] = (int) ((conversionfactor*((double)dh*1000)/
							((double)ds+conversionfactor*R*log(c)))-273.15*conversionfactor);

	}
	else { //set Tm for DNA-DNA to zero
	table[i][5] = 0;
	}

	//-----------------------------------------------------------------------------------------------------
	//--------------------------------------------------------------------------------------------
    //calculate the free energy of intramolecular folding of oligo

    
	if (Usesub==2) { //-2 is unavailabe, so not using partition function as interoligo cannot use it
		//reuse some arrays when scanning along the sequence for intramolecule
		//only calculate the partition function for the first one
		//cannot reuse array when prefilter is used
		if (prefilter==0) {
			if (i==start) {	
				intraoligo->reset4oligo(&oligo1);
				intraoligo->partition(true,&Qolig);
			
			}
			//reuse some arrays overlapped, only change the index of left binding site
			//be careful that the oligo is complementary to the target, copy direction is reversed
			else {
			
				for (ip=length-1;ip>2;ip--) {
					for (jp=length-1;jp>=ip;jp--) {
							
						intraoligo->wca[ip][jp]=intraoligo->wca[ip-1][jp-1];
						intraoligo->w->f(ip,jp)=intraoligo->w->f(ip-1,jp-1);
						intraoligo->v->f(ip,jp)=intraoligo->v->f(ip-1,jp-1);
						intraoligo->wmb->f(ip,jp)=intraoligo->wmb->f(ip-1,jp-1);
						intraoligo->wl->f(ip,jp)=intraoligo->wl->f(ip-1,jp-1);
						intraoligo->wmbl->f(ip,jp)=intraoligo->wmbl->f(ip-1,jp-1);
						intraoligo->wcoax->f(ip,jp)=intraoligo->wcoax->f(ip-1,jp-1);
					}
				}
				//oligo is complementary to target, set reverse=1 for scanfill()
				intraoligo->scanfill(&oligo1,&Qolig,1);

			}
		}
		else {
			intraoligo->reset4oligo(&oligo1);
			intraoligo->partition(true,&Qolig);
		}
		oligo1.energy[1]= (int)( conversionfactor*rt*
			((long double)length*logl( (intraoligo->data)->scaling) - logl((long double)Qolig)) );

	}
	else {//use the lowest free energy for Usesub 0 and 1 since oligo is small
    	if (isdna !=0)	dynamic(&oligo1,&ddata,1000,10,0,0,true);
		else	dynamic(&oligo1,&data,1000,10,0,0,true);
	}
    	
    table[i][3] = oligo1.energy[1];
	//if the intra structure is unfavorable , not consider it in the total energy
	if (table[i][3]>0) {
		table[i][3]=0;
    }

	//-----------------------------------------------------------------------------------------
    //calculate the free energy of intermolecular folding of oligo
    oligo.intermolecular = true;
    oligo.numofbases = 2*length + 3;
    for (j=1;j<=length;j++) {
		oligo.numseq[j+length+3] = oligo.numseq[j];
    }
    for (j=1;j<=3;j++) {
		oligo.numseq[length+j] = 5;
    }
    
	if (Usesub==2) {	//-2 is unavailabe, as partition function is not availabe for inter oligo	
			
		interoligo->reset4oligo(&oligo);
		interoligo->partition(true,&Qolig);
				
		oligo.energy[1]= (int)( conversionfactor*rt*
			( (long double)(oligo.numofbases)*logl((interoligo->data)->scaling)  - logl((long double)Qolig) ) );
	}
	else {
		if (isdna !=0 )	dynamic(&oligo,&ddata,1000,10,0,0,true);
		else	dynamic(&oligo,&data,1000,10,0,0,true);
	}
		      	
	table[i][4] = oligo.energy[1];
	//if the intra structure is unfavorable , not consider it in the total energy
    if (table[i][4]>0) {
		table[i][4]=0;
    }
	


	//-----------------------------------------------------------------------------------------
	//-----------------------------------------------------------------------------------------
    //calculate the overall free energy for the binding
    k1b = expl(((-1)*(long double)(table[i][4]))/((1.9872)*(3.1)));
    k1u = expl(((-1)*(long double)(table[i][3]))/((1.9872)*(3.1)));

	//I don't know what's going on here.:)  ---John , Nov.9,2005
	if (i==867) {
		table[0][0]=56;
    }

    if ((k1u/(c*k1b))<100) {

      	dgeff = (int)(-((1.9872)*(3.1)*logl(((4*k1b*(long double)(c))/(-1-k1u+sqrtl(powl(1+k1u,2)
						+8*k1b*(long double)(c))))-1)));
    }
    else {
		dgeff = table[i][3];

	}
    //this line converts the free energy to the convention explained in
    //the written version of the algorithm:
    table[i][2] = -table[i][2];

    if ((table[i][2]<0)&&(dgeff<0)) {
	   	table[i][0] = table[i][1] + (int)((rt*conversionfactor)*logl( (expl(-(long double)(table[i][2])/
					  (rt*conversionfactor))+1.0)*  (expl(-(long double)(dgeff)/(rt*conversionfactor))+1.0) )+0.5);
    }
    else if (table[i][2]<0) {
		table[i][0] = table[i][1] + (int)((rt*conversionfactor)*logl( (expl(-(long double)(table[i][2])/
					  (rt*conversionfactor))+1.0))+0.5);
    }
    else {
		table[i][0] = table[i][1] + (int)((rt*conversionfactor)*logl(  (expl(-(long double)(dgeff)/
					  (rt*conversionfactor))+1.0) )+0.5);
    }

    //table[i][0] = table[i][1] + table[i][2] - dgeff;

}
//The scan is finished now
//-------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------  



	//clean up the ct file
	ct->energy[1] = energy;

	//clean up memory use
	if (option==1) {
		if (Usesub==3 && foldsize==0) {
			ct->nnopair=0;
      		for (i=1;i<=ct->numofstructures;i++) delete[] temp2[i];
			delete[] temp2;
			delete[] energyarray;
		}
		else if (Usesub==0) delete[] temp;
		
	}
    else if (option ==2) {
   		ct->nnopair=0;
   		for (j=1;j<=ct->numofbases;j++) {
    		ct->basepr[1][j] = temp[j];
		}
		delete[] temp;
		if (Usesub==1)	delete[] ctenergy;
		if ( (Usesub==1||Usesub==3||Usesub==4) && foldsize==0 )	delete[] energyarray;
		else if (Usesub==2 || Usesub ==4) {
			if(isdna !=0 )	delete dpfdata;
			delete pfdata;
			delete target;
			if (foldsize>0 && prefilter->useit == 0)	delete targetcopy;
			if (Usesub==2) {
				delete intraoligo ;
				delete interoligo ;
			}
		}
	}
//cout <<"\n"<< numofstructures<<"\n";	
}


//=======================================================================
int readrd (rddata* data,char* dnarna) {
	
	int count,i,k,j,l;
	ifstream dr;
	char lineoftext[100];

	//make sure the file exists
	FILE *check;
	if ((check = fopen(dnarna, "r"))== NULL) {
	return 0;
	}

	fclose(check);
	dr.open(dnarna);
	/* Read info from stackdr */
	//add to the stack table the case where X (represented as 0) is looked up:
	for (count=1;count<=2;count++) dr >> lineoftext;//get past text in file
	dr >> lineoftext;
	data->init =(int)floor(conversionfactor*(atof(lineoftext)));

	for (count=1;count<=42;count++) dr >> lineoftext;//get past text in file
	for (i=0;i<=4;i++) {
		if (i!=0) for (count=1;count<=60;count++) dr >> lineoftext;
		for (k=0;k<=4;k++) {
			for (j=0;j<=4;j++) {
				for (l=0;l<=4;l++) {
					if ((i==0)||(j==0)||(k==0)||(l==0)) {
						data->stack[i][j][k][l]=0;
					}
					else {
						dr >> lineoftext;
						if (strcmp(lineoftext,".")){
							data->stack[i][j][k][l] =(int)floor(conversionfactor*(atof(lineoftext))+.5);
						}
						else data->stack[i][j][k][l] = infinity;
					}
					//cout <<"stack "<<i<<" "<<j<<" "<<k<<" "<<l<<"  "<<data->stack[i][j][k][l]<<"\n";
				}
				//cin >> m;
			}
		}
	}

	return 1;
}



//=======================================================================
//return the numerical equivalent of the complementary base to nucleotide i
int complement(int i, structure *ct) {
	
	int a;

	if (ct->numseq[i] == 0) return 0;
	else {
	   	a = 5 - ct->numseq[i];
		return a;
    }
}



//=======================================================================
//Function gets the names of data files to open for DNA-DNA parameters
void dgetdat(char *loop, char *stackf, char *tstackh, char *tstacki,
		     char *tloop, char *miscloop, char *danglef, char *int22,
			 char *int21,char *coax, char *tstackcoax,
		     char *coaxstack, char *tstack, char *tstackm,
			 char *triloop, char *int11, char *Path) {
	strcpy (loop,Path);
	strcpy (stackf,Path);
	strcpy (tstackh,Path);
	strcpy (tstacki,Path);
	strcpy (tloop,Path);
	strcpy (miscloop,Path);
	strcpy (danglef,Path);
	strcpy (int22,Path);
	strcpy (int21,Path);
	strcpy (triloop,Path);
	strcpy (coax,Path);
	strcpy (tstackcoax,Path);
	strcpy (coaxstack,Path);
	strcpy (tstack,Path);
	strcpy (tstackm,Path);
	strcpy (int11,Path);

	strcat (loop,"dnaloop.dat");
	strcat (stackf,"dnastack.dat");
	strcat (tstackh,"dnatstackh.dat");
	strcat (tstacki,"dnatstacki.dat");
	strcat (tloop,"dnatloop.dat");
	strcat (miscloop,"dnamiscloop.dat");
	strcat (danglef,"dnadangle.dat");
	strcat (int22,"dnaint22.dat");
	strcat (int21,"dnaint21.dat");
	strcat (triloop,"dnatriloop.dat");
	strcat (coax,"dnacoaxial.dat");
	strcat (tstackcoax,"dnatstackcoax.dat");
	strcat (coaxstack,"dnacoaxstack.dat");
	strcat (tstack,"dnatstack.dat");
	strcat (tstackm,"dnatstackm.dat");
	strcat (int11,"dnaint11.dat");
}



//=======================================================================
//This is a older post-filter of siRNA 
//use siRNA selection criteria to filter the output
//mask will contain true for those that meet the criteria
void filterbysirna ( structure *ct, int **table, int length, datatable *data, 
				     bool *mask, double asuf, double tofe, double fnnfe) {
	
	int iasuf,itofe,ifnnfe;
	int i,j,*k;
	k = new int [length];

	iasuf = (int) (conversionfactor*asuf);
	itofe = (int) (conversionfactor*tofe);
	ifnnfe = (int) (conversionfactor*fnnfe);


	for (i=1;i<=(ct->numofbases-length+1); i++) {
		//filter all oligos
		mask[i]=true;
		if (iasuf>table[i][3]) {
			mask[i]=false;
		}
		if (itofe>table[i][2]) {
			mask[i]=false;
		}

		//filter out those that have a repeat of A, G, or U of more than 4. 
		if (length>4) {
			for (j=length-1;j>=0;j--) {
   				k[j] = complement(i+j,ct);
      		
			}
			for (j=0;j<length-3;j++) {
				if (k[j]==1) {
					if (k[j+1]==1&&k[j+2]==1&&k[j+3]==1)	mask[i]=false;
				}
				else if (k[j]==3) {
					if (k[j+1]==3&&k[j+2]==3&&k[j+3]==3)	mask[i]=false;
				}
				else if (k[j]==4) {
					if (k[j+1]==4&&k[j+2]==4&&k[j+3]==4)	mask[i]=false;
				}
			}
		}
		table[i][6]= data->stack[complement(i+length-1,ct)][ct->numseq[i+length-1]]
								[complement(i+length-2,ct)][ct->numseq[i+length-2]];
		//account for change in AU end if necessary:
		if (ct->numseq[i+length-1]==1||ct->numseq[i+length-1]==4) {
			if (ct->numseq[i+length-1]==2||ct->numseq[i+length-1]==3) {
				table[i][6]+=data->auend;
			}
		}
		else {
			if (ct->numseq[i+length-1]==1||ct->numseq[i+length-1]==4) {
				table[i][6]-=data->auend;
			}
		}
		if (ifnnfe>table[i][6])
			mask[i]=false;

	}

	delete[] k;
}


//=======================================================================
//output tab delimited report
void report(char* filename, structure *ct, int **table,int **numofsubstructures, int length, int isdna,
			double conc, int Usesub,int start,int stop,siPREFILTER *prefilter,int foldsize, 
			bool *mask, double asuf, double tofe, double fnnfe) {
	
	int i,j,k;
	bool use;

	#if defined (debugmode)
	ofstream out(filename);
	#endif
	
	if (mask!=NULL) {
		cout << "Oligonucleotides filtered by siRNA selection criteria: ";
		cout << "\tAntisense strand unimolecular folding free energy >= "<<asuf<<"<br>\n";
		cout << "\tTarget opening free energy >= "<<tofe<<"<br>\n";
		cout << "\tFirst nearest neighbor free energy >= "<<fnnfe<<"<br>\n";
		cout << "\tSequences with more than three G's, U's, or A's in a row removed.<br>\n";
		cout << "Antisense strand shown 5' to 3'.<br>\n";
	}
	cout <<"<br><br>\n";
 	cout <<"<h3>Energy table:</h3><br>\n";
	cout <<"<table>\n";
	cout <<"<tr> <td>Pos.</td>\t<td>Oligo(5\'->3\')</td><td> </td>\t\t\t\t\t<td>Overall</td>\t";
	cout <<"<td>Duplex</td>\t<td>Tm-Dup</td>\t<td>Break-targ.</td>";
	cout <<"\t<td>Intraoligo</td>\t<td>Interoligo</td>\t<td>End_diff</td>\t<td>prefilter_score</td>";
	if (Usesub==1 || Usesub==3||Usesub==4)	cout <<"\t<td>structures #</td>\t<td>constrained structures #</td>";
	if (mask!=NULL) cout << "\t<td>First NN Opening</td>";
	cout <<"</tr><tr> </tr>\n\n";
	cout <<"\t\t\t\t\t\t<tr> <td> </td><td> </td> <td><td>kcal/mol</td>\t<td>kcal/mol</td>\t";
	cout <<"<td>degC</td>\t\t<td>kcal/mol</td>\t\t<td>kcal/mol</td>\t<td>kcal/mol</td>\t<td>kcal/mol</td>";
	if (mask!=NULL) "\t <td>kcal/mol</td>"; 
	cout << "</tr>\n";

	#if defined (debugmode)
			
		out <<"<tr> <td>Pos.</td>\t<td>Oligo</td><td> </td>\t\t\t\t\t<td>Overall</td>\t";
		out <<"<td>Duplex</td>\t<td>Tm-Dup</td>\t<td>Break-targ.</td>";
		out <<"\t<td>Intraoligo</td>\t<td>Interoligo</td>\t<td>End_diff</td>";
		if (mask!=NULL) cout << "\t<td>First NN Opening</td>";
		out <<"</tr><tr> </tr>\n\n";
		out <<"\t\t\t\t\t\t<tr> <td> </td><td> </td> <td><td>kcal/mol</td>\t<td>kcal/mol</td>\t";
		out <<"<td>degC</td>\t\t<td>kcal/mol</td>\t\t<td>kcal/mol</td>\t<td>kcal/mol</td>\t<td>kcal/mol</td>";
		if (mask!=NULL) "\t <td>kcal/mol</td>"; 
		out << "</tr>\n";
	#endif

	for (i=start;i<=stop; i++) {

		if (prefilter->useit != 0) {
			if (prefilter->score[i] < FILTER_PASS) {
				continue;
			}
		}

 
		//output each possible oligo
		if (mask!=NULL) use = mask[i];
		else use = true;
		
		if (use) {
			cout << "<tr><td>"<<i <<"</td>\t";
			cout << "<td>";
			for (j=length-1;j>=0;j--) {
   				k = complement(i+j,ct);
      			if (isdna != 0) {
         			if (k==0) cout << "X";
					else if (k==1) cout << "A";
					else if (k==2) cout << "C";
					else if (k==3) cout << "G";
					else if (k==4) cout << "T";
				}
				else {
		     		if (k==0) cout << "X";
					else if (k==1) cout << "A";
					else if (k==2) cout << "C";
					else if (k==3) cout << "G";
					else if (k==4) cout << "U";
				}
			}
			cout << "</td><td> </td>\t\t";
		
			if (prefilter->useit != 0) {
		
				if (prefilter->score[i] < FILTER_PASS) {
					cout<<"<td>-</td>\t\t<td>-</td>\t\t<td>-</td>\t\t<td>-</td>\t\t\t<td>-</td>\t\t<td>-</td>\t\t</tr>\n";
					continue;
				}
			}

			cout <<"<td>"<< ((double) (table[i][0]))/conversionfactor << "</td>\t\t";
			cout <<"<td>"<< ((double) (table[i][1]))/conversionfactor << "</td>\t\t";
			cout <<"<td>"<< ((double) (table[i][5]))/conversionfactor << "</td>\t\t";
			cout <<"<td>"<< ((double) (table[i][2]))/conversionfactor << "</td>\t\t\t";
			cout <<"<td>"<< ((double) (table[i][3]))/conversionfactor << "</td>\t\t";
			cout <<"<td>"<< ((double) (table[i][4]))/conversionfactor << "</td>\t\t";
			
			#if defined (debugmode)
				out <<"<td>"<< ((double) (table[i][0]))/conversionfactor << "</td>\t\t";
				out <<"<td>"<< ((double) (table[i][1]))/conversionfactor << "</td>\t\t";
				out <<"<td>"<< ((double) (table[i][5]))/conversionfactor << "</td>\t\t";
				out <<"<td>"<< ((double) (table[i][2]))/conversionfactor << "</td>\t\t\t";
				out <<"<td>"<< ((double) (table[i][3]))/conversionfactor << "</td>\t\t";
				out <<"<td>"<< ((double) (table[i][4]))/conversionfactor << "</td>\t\t";
			#endif



			if (prefilter->useit != 0)	{
				cout <<"<td>"<<  prefilter->enddiff[i] <<"</td>\t\t";
				cout <<"<td>"<<  prefilter->score[i] <<"</td>\t\t";
			}
			if (Usesub==1 || Usesub==3||Usesub==4)	 {
				cout <<"<td>"<< numofsubstructures[i][0]  <<"</td>\t\t";
				cout <<"<td>"<< numofsubstructures[i][1]  <<"</td>\t\t";
			}
			
			if (mask!=NULL) {
		
				cout <<"<td>"<< ((double) (table[i][6])/conversionfactor) << "</td>\t\t";
		
			}
			cout << "</tr>\n";

			#if defined (debugmode)
			out << "\n";
			#endif
		}
	}
	cout<< "</table>\n";
//	cout<< "</body></html>\n";

	#if defined (debugmode)
	out.close();
	#endif



}

//=================================================================================================
//copy the arrays to be reused with different index when the folded region move one nucleotide to the right
inline void scancopy(OligoPclass *region, OligoPclass *copyregion) {

	int j,i;
	//only copy the overlapped region which is in the middle of sequence
	for (i=2;i<=(copyregion->number)-2;i++) {
		for (j=i;j<=(copyregion->number)-2;j++) {
			//w5 is not copied as it will be recalc when i==1						
			copyregion->wca[i][j]=region->wca[i+1][j+1];
			copyregion->w->f(i,j)=region->w->f(i+1,j+1);
			copyregion->v->f(i,j)=region->v->f(i+1,j+1);
			copyregion->wmb->f(i,j)=region->wmb->f(i+1,j+1);
			copyregion->wl->f(i,j)=region->wl->f(i+1,j+1);
			copyregion->wmbl->f(i,j)=region->wmbl->f(i+1,j+1);
			copyregion->wcoax->f(i,j)=region->wcoax->f(i+1,j+1);
			
			

		}
	}
}
//copy every arrays to be reused when the folded region did not move right yet 
inline void scancopyend(OligoPclass *region, OligoPclass *copyregion) {

	int j,i;
	for (i=1;i<=copyregion->number;i++) {
		for (j=i;j<=copyregion->number;j++) {
						
			if(i==1)	copyregion->w5[j]=region->w5[j];
			copyregion->wca[i][j]=region->wca[i][j];
			copyregion->w->f(i,j)=region->w->f(i,j);
			copyregion->v->f(i,j)=region->v->f(i,j);
			copyregion->wmb->f(i,j)=region->wmb->f(i,j);
			copyregion->wl->f(i,j)=region->wl->f(i,j);
			copyregion->wmbl->f(i,j)=region->wmbl->f(i,j);
			copyregion->wcoax->f(i,j)=region->wcoax->f(i,j);
		
		}
	}
}


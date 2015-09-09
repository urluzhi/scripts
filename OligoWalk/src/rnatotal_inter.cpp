/************************************************************************
* rnatotal_inter.cpp:							*
* Usage: foo.o  [lisfile] [report file] [extension]			*
* Input: .lis .seq .dat .dh .ct						*
* Calculate free energy with N3 instead of N4 for  unlimited internal	*
* loop; using algorithm_inter.cpp file instead of algorithm.cpp  	*
* based on rnatotal_tem.cpp:						*
* Calculate Free energy at certain temperature  with 			*
* dg = dh - ds*T= dh - (dh-dg37)*T/310.15				*   	
* temperatures are read from lis files                                  *   
* dh are read from .dh files and stored in Class dhdata			*
* dg37 are read from .dat files and stored in Class data		*
* dg are calculated and are stored in Class dg				*
*************************************************************************
* Created:	2003;			Modified:	Dec.2005	*
* Copyright: Zhi John Lu & David Mathews				*	
*************************************************************************/

#include "stdafx.h"
#include "algorithm_inter.cpp"
#include "rna_library_inter.cpp"
#include "structure.cpp"
#include "dG_tem.cpp"
#define maxsequences 2000 //maximum number of sequences that can be folded
#define cntrl8 20
#define cntrl6 750
#define cntrl9 0
#define writefiles 1 
#define path "../../archive/" //path to archive of known structures
//#define path ""

void scorer(structure* correct, structure* test, int *score, int *basepairs);

void bpscorer(structure *correctct,structure *ct,int *bpscore);

void getdat(char *loop, char *stackf, char *tstackh, char *tstacki,
		char *tloop, char *miscloop, char *danglef, char *int22,
      char *int21,char *coax, char *tstackcoax,
      char *coaxstack, char *tstack,char *tstackm,char *triloop, char *int11, char *hexaloop,
	  char *tstacki23, char *tstacki1n);
void getdh(char *loop, char *stackf, char *tstackh, char *tstacki,
                char *tloop, char *miscloop, char *danglef, char *int22,
		      char *int21,char *coax, char *tstackcoax,
		            char *coaxstack, char *tstack,char *tstackm,char *triloop, char *int11, char *hexaloop,
			              char *tstacki23, char *tstacki1n);

void getoutfile (char *outfile);//gets name of a line printer output file
											//to which a structure will be written
//void getinfo (int *cntrl6,int *cntrl8, int *cntrl9);//get info for sub-otimal structure


int main(int argc, char *argv[]) {


	char loop[maxfil],stackf[maxfil],tstackh[maxfil],tstacki[maxfil],
		tloop[maxfil],miscloop[maxfil],danglef[maxfil],int22[maxfil],
      int21[maxfil],coax[maxfil],tstackcoax[maxfil],
      coaxstack[maxfil],tstack[maxfil],tstackm[maxfil],triloop[maxfil],int11[maxfil],hexaloop[maxfil],
	  tstacki23[maxfil], tstacki1n[maxfil];
   //int cntrl6,cntrl8,cntrl9,ir;
   datatable *data,*dhdata,*dg;
   structure *ct,*correctct;
   //TProgressDialog PD;
   
   
   char list[maxfil],outfile[maxfil],label[maxfil],prefix[maxsequences][maxfil];
   char sequence[maxfil],ctin[maxfil],ctoutdyn[maxfil],ctoutefn[maxfil];
   int score[maxstructures],basepairs,percent,mine;
   int totalbp,totaldyn,totalefn,totalbest,totalworst,totalbpwise;
   int dynenergy[maxstructures],dynnumber[maxstructures];
   register int c,i,j,ir;
   int cur,tempnumofstructures,bpscore;
   int efnenergy[maxstructures],efnnumber[maxstructures];
   TProgressDialog PD;
   float temperature[maxsequences],T[maxsequences];// degree for temperature and K for T
   double totaldynp, totalefnp, totalbestp, lowdyn, highdyn, lowefn, highefn,
   	lowbest, highbest, sxd, sxsd, sd, per, sxe, sxse, sxb, sxsb, sxbp, sxsbp,
      totalbpp,totaldiff,sxdiff,sxsdiff;
   //dynamic alloc datatalbes
   data=new datatable();
   dhdata=new datatable();
   dg=new datatable();
   
   //initialize the variables used for scoring the accuracy:
   totalbp = 0;
   totaldyn = 0;
   totalefn = 0;
   totalbest = 0;
   totalworst = 0;
   totaldynp = 0;
   totalefnp = 0;
   totalbestp = 0;
   lowdyn = 1;
   highdyn = 0;
   lowefn = 1;
   highefn = 0;
   lowbest = 1;
   highbest = 0;
   sxe = 0;
   sxse = 0;
   sxd = 0;
   sxsd = 0;
   sxb = 0;
   sxsb = 0;
   totalbpwise = 0;
   sxb = 0;
   sxsb = 0;
   totalbpp = 0;
   sxbp = 0;
   sxsbp = 0;
   totaldiff = 0;
   sxdiff = 0;
   sxsdiff = 0;

   //Get required information:
   if (argc!=4) {
   	cout << "Usage: rnatotal [lisfile] [report file] [extension]\n";



   	cout << "Enter the name of a .lis file: \n";
   	cin >> list;

   	cout << "Enter the name for a report file to be written:  \n";
   	cin >> outfile;

   	cout << "Enter an extension for this set of foldings:  \n";
   	cin >> label;
   }
   else {
   	  strcpy(list,argv[1]);
      strcpy(outfile,argv[2]);
	  strcat(outfile,".out");
      strcpy(label,argv[3]);
   }



   ifstream inf;
	inf.open(list);

	//out will contain a report about accuracy for this set of structure predictions
	ofstream out;
	out.open(outfile);

   out <<"RNAstructure Folding Report\n";
   out <<"List file:  "<<list<<"\n";
   out <<"For set of structures with the suffix: "<<label<<"\n";
   out <<"% sort: "<<cntrl8<<"\n";
   out <<"# structures: "<<cntrl6<<"\n";
   out <<"Window: "<<cntrl9<<"\n";



   //read the list of structures to be folded from a list file:
	for (ir=1;ir<=maxsequences;ir++) {
		
		inf >> prefix[ir];
		inf >> temperature[ir];    //read the temperatures of RNA at degree
		if(inf.eof()) break;
		T[ir]=temperature[ir]+273.15;// transfer degree to K
		cout << "ir = "<<ir<<" ct = "<<prefix[ir]<<" temp="<<temperature[ir]<<" T="<<T[ir]<<"\n"<<flush;

}

	ir--;


	//open the data files -- must reside with executable
   getdat (loop,stackf,tstackh,tstacki,tloop,miscloop,danglef,
   	int22,int21,coax,tstackcoax,coaxstack,tstack,tstackm,triloop,int11,hexaloop,tstacki23, tstacki1n);
	if (opendat (loop,stackf,tstackh,tstacki,tloop,miscloop,danglef,int22,int21,
   	coax,tstackcoax,coaxstack,tstack,tstackm,triloop,int11,hexaloop,tstacki23, tstacki1n,data)==0) {
      	cout << "A data file was lost";
         cin >> i;
   }
   //open the enthalpy files
 getdh (loop,stackf,tstackh,tstacki,tloop,miscloop,danglef,
           int22,int21,coax,tstackcoax,coaxstack,tstack,tstackm,triloop,int11,hexaloop,tstacki23, tstacki1n);
	if (opendat (loop,stackf,tstackh,tstacki,tloop,miscloop,danglef,int22,int21,
        coax,tstackcoax,coaxstack,tstack,tstackm,triloop,int11,hexaloop,tstacki23, tstacki1n,dhdata)==0) {
	cout << "A dhdata file was lost";
	 cin >> i;
   }

			    
   for (i=1;i<=ir;i++) {//fold and score every structure


   	ct = new structure;
      correctct = new structure;

      strcpy(sequence,path);
      strcat(sequence,prefix[i]);
      strcat(sequence,".seq");

    cout << i << "  "<<sequence<<"\n"<<flush;

      strcpy(ctin,path);
      strcat(ctin,prefix[i]);
      strcat(ctin,".ct");
      strcpy(ctoutdyn,prefix[i]);
      strcat(ctoutdyn,"_dyn");
      strcat(ctoutdyn,label);
      strcat(ctoutdyn,".ct");
      strcpy(ctoutefn,prefix[i]);
      strcat(ctoutefn,"_efn");
      strcat(ctoutefn,label);
      strcat(ctoutefn,".ct");


		//open the current sequence
      openseq (ct,sequence);
//calculated the free energy at temperature T[i]

dG_T(T[i],*data,*dhdata,*dg);

  
//predict the secondary structure
dynamic (ct,dg,cntrl6,cntrl8,cntrl9,&PD);

//open the ct with known structure
      openct (correctct,ctin);

     //score the preeiction 
      scorer(correctct,ct,score,&basepairs);

      
      totalbp = totalbp+basepairs;
      totaldyn = totaldyn + score[1];

      //output the dynamic algorithm information to the report:
      out << "\f\n\nStructure:  "<<prefix[i]<<"\n\n\n";
      out << "Temperature:   "<<temperature[i]<<"\n\n\n";
	  out << "# of basepairs:  "<<basepairs<<"\n\n";
      out << "Dynamic algorithm information: \n\n";
      out <<"#\tenergy\tscore\t%\n";
      percent = 100*score[1]/basepairs;
      out << "1\t"<<((double (ct->energy[1]))/10)<<"\t"<<score[1]<<"\t"<<percent<<"\n";
      percent = 100*score[2]/basepairs;
      out << "2\t"<<((double (ct->energy[2]))/10)<<"\t"<<score[2]<<"\t"<<percent<<"\n";
      percent = 100*score[3]/basepairs;
      out << "3\t"<<((double (ct->energy[3]))/10)<<"\t"<<score[3]<<"\t"<<percent<<"\n";
      percent = 100*score[4]/basepairs;
      out << "4\t"<<((double (ct->energy[4]))/10)<<"\t"<<score[4]<<"\t"<<percent<<"\n\n";

      per = double(score[1]) / double (basepairs);
      if (per<lowdyn) lowdyn = per;
      if (per>highdyn) highdyn = per;
      totaldynp = totaldynp + per;
      sxd = sxd + per;
      sxsd = sxsd + per*per;


      tempnumofstructures = ct->numofstructures;
      ct->numofstructures=1;
      if (writefiles) ctout(ct,ctoutdyn);
      ct->numofstructures=tempnumofstructures;

      //we're now going to sort by efn2 energies, but first, make an array
      //	to store dynamic number and dynamic energy:

      for (j=1;j<=ct->numofstructures;j++) {
       	dynnumber[j]=j;
         dynenergy[j]=ct->energy[j];

      }

    //predict efn2 energy -- note that this is archaic because coaxial
	  //stacking was added to the dynamic programming algorithm in RNAstructure 4
   	efn2 (dg,ct);
	
    
		//sort by efn2 energy
		for (c = 2; c<=(ct->numofstructures);c++){

			cur = c;

			while (cur>1) {
				if ((ct->energy[cur])<(ct->energy[cur-1])) {
      			swap(&ct->energy[cur],&ct->energy[cur-1]);
               swap(&(dynnumber[cur]),&(dynnumber[cur-1]));
               swap(&(dynenergy[cur]),&(dynenergy[cur-1]));
               swap(&(score[cur]),&(score[cur-1]));
         		for (j=1;j<=(ct->numofbases);j++) {
         			swap(&ct->basepr[cur][j],&ct->basepr[cur-1][j]);
         		}
         		cur--;
   			}
				else {
					break;
				}
			}
		}

      tempnumofstructures = ct->numofstructures;
      ct->numofstructures=1;
      if (writefiles) ctout(ct,ctoutefn);
      ct->numofstructures=tempnumofstructures;

      totalefn = totalefn + score[1];

      //output the efn2 info to the report
      out << "efn2 information: \n\n";
      out <<"#\tenergy\tscore\t%\tdyn #\tdyn E\n";
      percent = 100*score[1]/basepairs;
      out << "1\t"<<((double (ct->energy[1]))/10)<<"\t"<<score[1]
      	<<"\t"<<percent<<"\t"<<dynnumber[1]<<"\t"<<
         (double (dynenergy[1])/10)<<"\n";
      percent = 100*score[2]/basepairs;
      out << "2\t"<<((double (ct->energy[2]))/10)<<"\t"<<score[2]
      	<<"\t"<<percent<<"\t"<<dynnumber[2]<<"\t"<<
         (double (dynenergy[2])/10)<<"\n";
      percent = 100*score[3]/basepairs;
      out << "3\t"<<((double (ct->energy[3]))/10)<<"\t"<<score[3]
      	<<"\t"<<percent<<"\t"<<dynnumber[3]<<"\t"<<
         (double (dynenergy[3])/10)<<"\n";
      percent = 100*score[4]/basepairs;
      out << "4\t"<<((double (ct->energy[4]))/10)<<"\t"<<score[4]
      	<<"\t"<<percent<<"\t"<<dynnumber[4]<<"\t"<<
         (double (dynenergy[4])/10)<<"\n\n";


      //save lowest free energy in mine:
      mine = ct->energy[1];


      per = double(score[1]) / double (basepairs);
      if (per<lowefn) lowefn = per;
      if (per>highefn) highefn = per;
      totalefnp = totalefnp + per;
      sxe = sxe + per;
      sxse = sxse + per*per;

      //now sort by score for the best suboptimal structure:
      for (j=1;j<=ct->numofstructures;j++) {
       	efnnumber[j]=j;
         efnenergy[j]=ct->energy[j];

      }

      for (c = 2; c<=(ct->numofstructures);c++){

			cur = c;

			while (cur>1) {
				if ((score[cur])>(score[cur-1])) {
      			swap(&(efnnumber[cur]),&(efnnumber[cur-1]));
               swap(&(dynnumber[cur]),&(dynnumber[cur-1]));
               swap(&(dynenergy[cur]),&(dynenergy[cur-1]));
               swap(&(efnenergy[cur]),&(efnenergy[cur-1]));
               swap(&(score[cur]),&(score[cur-1]));
         		//for (j=1;j<=(ct.numofbases);j++) {
         		//	swap(&ct.basepr[cur][j],&ct.basepr[cur-1][j]);
         		//}
         		cur--;
   			}
				else {
					break;
				}
			}
		}

      totalbest = totalbest + score[1];

      //now output this info to a report:
      out << "best structure information: \n\n";
      out <<"#\tscore\t%\tefn #\tefn E\tdyn #\tdyn E\n";
      percent = 100*score[1]/basepairs;
      out << "1\t"<<score[1]<<"\t"<<percent<<"\t"<<efnnumber[1]<<"\t"<<
      	(((double) (efnenergy[1]))/10)<<"\t"<<dynnumber[1]<<"\t"<<
         (((double) (dynenergy[1]))/10)<<"\n";
      percent = 100*score[2]/basepairs;
      out << "2\t"<<score[2]<<"\t"<<percent<<"\t"<<efnnumber[2]<<"\t"<<
      	(((double) (efnenergy[2]))/10)<<"\t"<<dynnumber[2]<<"\t"<<
         (((double) (dynenergy[2]))/10)<<"\n";
      percent = 100*score[3]/basepairs;
      out << "3\t"<<score[3]<<"\t"<<percent<<"\t"<<efnnumber[3]<<"\t"<<
      	(((double) (efnenergy[3]))/10)<<"\t"<<dynnumber[3]<<"\t"<<
         (((double) (dynenergy[3]))/10)<<"\n";
      percent = 100*score[4]/basepairs;
      out << "4\t"<<score[4]<<"\t"<<percent<<"\t"<<efnnumber[4]<<"\t"<<
      	(((double) (efnenergy[4]))/10)<<"\t"<<dynnumber[4]<<"\t"<<
         (((double) (dynenergy[4]))/10)<<"\n";

      per = (double)(score[1]) / (double) (basepairs);
      if (per<lowbest) lowbest = per;
      if (per>highbest) highbest = per;
      totalbestp = totalbestp + per;
      sxb = sxb + per;
      sxsb = sxsb + per*per;

      //now calculate the % diff in energy between the best and the lfe structure
      per = ((double) (mine - efnenergy[1]))/(((double)(mine+efnenergy[1]))/2);
      totaldiff = totaldiff + per;
      sxdiff = sxdiff + per;
      sxsdiff = sxsdiff + per*per;


      //now sort for the worst structure
      for (c = 2; c<=(ct->numofstructures);c++){

			cur = c;

			while (cur>1) {
				if ((score[cur])<(score[cur-1])) {
      			swap(&(efnnumber[cur]),&(efnnumber[cur-1]));
               swap(&(dynnumber[cur]),&(dynnumber[cur-1]));
               swap(&(dynenergy[cur]),&(dynenergy[cur-1]));
               swap(&(efnenergy[cur]),&(efnenergy[cur-1]));
               swap(&(score[cur]),&(score[cur-1]));
         		//for (j=1;j<=(ct.numofbases);j++) {
         		//	swap(&ct.basepr[cur][j],&ct.basepr[cur-1][j]);
         		//}
         		cur--;
   			}
				else {
					break;
				}
			}
		}


      totalworst = totalworst + score[1];

      //now output this info to a report:
      out << "worst structure information: \n\n";
      out <<"#\tscore\t%\tefn #\tefn E\tdyn #\tdyn E\n";
      percent = 100*score[1]/basepairs;
      out << "1\t"<<score[1]<<"\t"<<percent<<"\t"<<efnnumber[1]<<"\t"<<
      	(((double) (efnenergy[1]))/10)<<"\t"<<dynnumber[1]<<"\t"<<
         (((double) (dynenergy[1]))/10)<<"\n";
      percent = 100*score[2]/basepairs;
      out << "2\t"<<score[2]<<"\t"<<percent<<"\t"<<efnnumber[2]<<"\t"<<
      	(((double) (efnenergy[2]))/10)<<"\t"<<dynnumber[2]<<"\t"<<
         (((double) (dynenergy[2]))/10)<<"\n";
      percent = 100*score[3]/basepairs;
      out << "3\t"<<score[3]<<"\t"<<percent<<"\t"<<efnnumber[3]<<"\t"<<
      	(((double) (efnenergy[3]))/10)<<"\t"<<dynnumber[3]<<"\t"<<
         (((double) (dynenergy[3]))/10)<<"\n";
      percent = 100*score[4]/basepairs;
      out << "4\t"<<score[4]<<"\t"<<percent<<"\t"<<efnnumber[4]<<"\t"<<
      	(((double) (efnenergy[4]))/10)<<"\t"<<dynnumber[4]<<"\t"<<
         (((double) (dynenergy[4]))/10)<<"\n";


      bpscorer(correctct,ct,&bpscore);
      totalbpwise = totalbpwise + bpscore;

      per = (double)(bpscore) / (double) (basepairs);
      totalbpp = totalbpp + per;
      sxbp = sxbp + per;
      sxsbp = sxsbp + per*per;

      out << "\nBP wise score: "<<per<<"\n\n";

   		delete ct;
      delete correctct;
	}

   //now ouput the summary page:

   out <<"\fSummary:\n\n\n";
   out <<"List file:  "<<list<<"\n\n\n";
   out <<"total basepairs in database: "<<totalbp<<"\n\n";

   percent = (int) (100.0*(((double) (totaldyn))/((double) (totalbp))));
   out <<"dyn: "<<"\t"<<totaldyn<<"\t"<<percent<<"\n\n";

   percent = (int) (100.0*(((double) (totalefn))/((double) (totalbp))));
   out <<"efn2: "<<"\t"<<totalefn<<"\t"<<percent<<"\n\n";


   percent = (int) (100.0*(((double) (totalbest))/((double) (totalbp))));
   out <<"best: "<<"\t"<<totalbest<<"\t"<<percent<<"\n\n";

   percent = (int) (100.0*(((double) (totalworst))/((double) (totalbp))));
   out <<"worst: "<<"\t"<<totalworst<<"\t"<<percent<<"\n\n";

   percent = (int) (100.0*(((double) (totalbpwise))/((double) (totalbp))));
   out <<"bp wise: "<<"\t"<<totalbpwise<<"\t"<<percent<<"\n\n";

   //ir = number of structures
   out << "Calculation by averaging each structure: \n\n";
   per = totaldynp/ir;
   sd = sqrt((ir*sxsd - sxd*sxd)/(ir*(ir-1)));
   out << "dyn ->: "<<"\t"<<per<<"  +/-  "<<sd<<"\n\n";

   per = totalefnp/ir;
   sd = sqrt((ir*sxse - sxe*sxe)/(ir*(ir-1)));
   out << "efn ->: "<<"\t"<<per<<"  +/-  "<<sd<<"\n\n";

   per = totalbestp/ir;
   sd = sqrt((ir*sxsb - sxb*sxb)/(ir*(ir-1)));
   out << "best ->: "<<"\t"<<per<<"  +/-  "<<sd<<"\n\n";

   per = totalbpp / ir;
   sd = sqrt((ir*sxsbp - sxbp*sxbp)/(ir*(ir-1)));
   out << "bp wise score ->: "<<"\t"<<per<<"  +/-  "<<sd<<"\n\n";

   per = totaldiff / ir;
   sd = sqrt((ir*sxsdiff - sxdiff*sxdiff)/(ir*(ir-1)));
   out << "Average energy diff between mfe and best: "<<"\t"<<per<<"  +/-  "<<sd<<"\n\n";

   out << "\nRanges of percent:  \n\n";

   out << "dynamic from "<<lowdyn<< " to "<< highdyn<<"\n\n";
   out << "efn2 from "<<lowefn<< " to "<< highefn<<"\n\n";
   out << "best from "<<lowbest<< " to "<< highbest<<"\n\n";
   delete data;
   delete dhdata;
   delete dg;
   return 0;

}


void bpscorer(structure *correctct,structure *test,int *bpscore) {
//This routine will look to see if each bp in the phylogenetic structure
//	occurs in at least one suboptimal structure

int i,j;
*bpscore = 0;

if (correctct->numofstructures!=1) {
	cout << "There is more than one structure in the phylogenetic ct file.";
   cout << "\nThis is not allowed.\n";
   //cin >> correctct;
}
if (correctct->numofbases!=test->numofbases) {
	cout << "The two ct files have different lengths.\n";
   cout << "This is not allowed.";
   //cin >> correctct;
}


for (i=1;i<=correctct->numofbases;i++) {//look for each bp
	if (correctct->basepr[1][i]>i) {//bp found
   	j = 1;
      while (j<=test->numofstructures) {

       			if (test->basepr[j][i]==correctct->basepr[1][i]) {
               	(*bpscore)++;
                  j = test->numofstructures+10;
               }

         		else if ((test->basepr[j][i]+1)==correctct->basepr[1][i]) {
               	(*bpscore)++;
                  j = test->numofstructures+10;
               }

         		else if ((test->basepr[j][i]-1)==correctct->basepr[1][i]){
               	(*bpscore)++;
                  j = test->numofstructures+10;
               }

					else if ((test->basepr[j][i+1])==correctct->basepr[1][i]){
               	(*bpscore)++;
                  j = test->numofstructures+10;
               }

					 else if (i-1)
						{
						if (test->basepr[j][i-1]==correctct->basepr[1][i]) 
							{
               				(*bpscore)++;
							j = test->numofstructures+10;
               				}
						}
               j++;



      }


   }
}





}

void scorer(structure* correct, structure* test, int *score, int *basepairs){
int i,j;


if (correct->numofstructures!=1) {
	cout << "There is more than one structure in the phylogenetic ct file.";
   cout << "\nThis is not allowed.\n";
   //cin >> correctct;
}
if (correct->numofbases!=test->numofbases) {
	cout << "The two ct files have different lengths.\n";
   cout << "This is not allowed.";
   //cin >> correctct;
}


(*basepairs)=0;
//count the bases in the phylogenetic ct:
for (i=1;i<=correct->numofbases;i++) {
	if (correct->basepr[1][i]>i) (*basepairs)++;
}

//Now check every pair
for (j=1;j<=test->numofstructures;j++) {
	//cout << "j= "<<j<<"\n";
	//cin >> basepairs;
	score[j]=0;
	for (i=1;i<=test->numofbases;i++) {
		//cout<<"score = "<<score[j]<<"\n";
		//cout<<"pairs = "<<test.basepr[j][i]<<"  "<<correct.basepr[j][i]<<"\n";
    		if (correct->basepr[1][i]>i) {
      			if (test->basepr[j][i]==correct->basepr[1][i])
         			(score[j])++;
         		else if ((test->basepr[j][i]+1)==correct->basepr[1][i])
         			(score[j])++;
         		else if ((test->basepr[j][i]-1)==correct->basepr[1][i])
         			(score[j])++;
			else if ((test->basepr[j][i+1])==correct->basepr[1][i])
				(score[j])++;
			else if (i-1)
				{
				if((test->basepr[j][i-1])==correct->basepr[1][i])
				(score[j])++;
				}	
      		}
   	}
}

}


/*	Function getdat

	Function gets the names of data files to open

*/

void getdat(char *loop, char *stackf, char *tstackh, char *tstacki,
		char *tloop, char *miscloop, char *danglef, char *int22,
      char *int21,char *coax, char *tstackcoax,
      char *coaxstack, char *tstack, char *tstackm,char *triloop, char *int11,
	  char *hexaloop,char *tstacki23, char *tstacki1n )

{
strcpy (loop,"/home/john/data/loop.dat");
strcpy (stackf,"/home/john/data/stack.dat");
strcpy (tstackh,"/home/john/data/tstackh.dat");
strcpy (tstacki,"/home/john/data/tstacki.dat");
strcpy (tloop,"/home/john/data/tloop.dat");
strcpy (miscloop,"/home/john/data/miscloop.dat");
strcpy (danglef,"/home/john/data/dangle.dat");
strcpy (int22,"/home/john/data/int22.dat");
strcpy (int21,"/home/john/data/int21.dat");
strcpy (coax,"/home/john/data/coaxial.dat");
strcpy (tstackcoax,"/home/john/data/tstackcoax.dat");
strcpy (coaxstack,"/home/john/data/coaxstack.dat");
strcpy (tstack,"/home/john/data/tstack.dat");
strcpy (tstackm,"/home/john/data/tstackm.dat");
strcpy (triloop,"/home/john/data/triloop.dat");
strcpy (int11,"/home/john/data/int11.dat");
strcpy (hexaloop,"/home/john/data/hexaloop.dat");
strcpy (tstacki23,"/home/john/data/tstacki23.dat");
strcpy (tstacki1n,"/home/john/data/tstacki1n.dat");
}

/*      Function getdh      
Function gets the names of dhdata files to open */

void getdh(char *loop, char *stackf, char *tstackh, char *tstacki,char *tloop, char *miscloop, char *danglef, char *int22,char *int21,char *coax, char *tstackcoax,char *coaxstack, char *tstack, char *tstackm,char *triloop, char *int11,char *hexaloop,char *tstacki23, char *tstacki1n )
{
strcpy (loop,"/home/john/data/loop.dh");
strcpy (stackf,"/home/john/data/stack.dh");
strcpy (tstackh,"/home/john/data/tstackh.dh");
strcpy (tstacki,"/home/john/data/tstacki.dh");
strcpy (tloop,"/home/john/data/tloop.dh");
strcpy (miscloop,"/home/john/data/miscloop.dh");
strcpy (danglef,"/home/john/data/dangle.dh");
strcpy (int22,"/home/john/data/int22.dh");
strcpy (int21,"/home/john/data/int21.dh");
strcpy (coax,"/home/john/data/coaxial.dh");
strcpy (tstackcoax,"/home/john/data/tstackcoax.dh");
strcpy (coaxstack,"/home/john/data/coaxstack.dh");
strcpy (tstack,"/home/john/data/tstack.dh");
strcpy (tstackm,"/home/john/data/tstackm.dh");
strcpy (triloop,"/home/john/data/triloop.dh");
strcpy (int11,"/home/john/data/int11.dh");
strcpy (hexaloop,"/home/john/data/hexaloop.dh");
strcpy (tstacki23,"/home/john/data/tstacki23.dh");
strcpy (tstacki1n,"/home/john/data/tstacki1n.dh");
}




/*void getinfo (int *cntrl6,int *cntrl8, int *cntrl9) {

cout << "Enter % for sort (as a integer):    ";
cin >> *cntrl8;
cout << "Enter the maximum number of structures:     ";
cin >> *cntrl6;
cout << "Enter the window size:    ";
cin >> *cntrl9;
return;
} */

void errmsg(int err,int erri) {

if (err==30) {
	cout << "End Reached at traceback #"<<erri<<"\n";
   exit(1);
}
if (err==100) {
	cout << "error # "<<erri;
   exit(1);
}
switch (err) {
	case 1:
   	cout << "Could not allocate enough memory";
      break;
   case 2:
   	cout << "Too many possible base pairs";
      break;
   case 3:
   	cout << "Too many helixes in multibranch loop";
   case 4:
   	cout << "Too many structures in CT file";
   default:
   	cout << "Unknown error";
}
cin >> err;
exit(1);
return;

}

void update(int i) {

	cout<< i<<"\n"<<flush;
}


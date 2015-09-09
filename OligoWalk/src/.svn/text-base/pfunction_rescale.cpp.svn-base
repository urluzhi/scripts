#include "stdafx.h"
#include "pfunction_rescale.h"
#include <math.h>

//#undef pfdebugmode  //flag to indicate debugging
#define pfdebugmode  //flag to indicate debugging

#undef equiout
//#define equiout //flag to indicate equilibria should be written to file

#undef timer
//#define timer //flag to indicate the code execution should be timed

#undef disablecoax
//#define disablecoax


inline void write(ofstream *out,double *i) {

	out->write((char *) i,sizeof(*i));
}

inline void read(ifstream *out,double *i) {

	out->read((char *) i,sizeof(*i));
}

//This function cacluates a partition function for the sequence in CT
//If quickQ == true, return the partition function value in Q
//	otherwise, save the partial partition functions in the datafile named save
//If updates on progress are unwanted, set update=NULL
void pfunction(structure* ct,pfdatatable* data, TProgressDialog* update, char* save, bool quickQ, PFPRECISION *Q)
{

		
int d,ip,jp,ii,jj,lowi;
register int i,j;
int k,l,m,n,o,p;
register int inc[6][6]={{0,0,0,0,0,0},{0,0,0,0,1,0},{0,0,0,1,0,0},{0,0,1,0,1,0},
	{0,1,0,1,0,0},{0,0,0,0,0,0}};
bool *lfce,*mod;//[maxbases+1][maxbases+1];
int before,after;
register PFPRECISION e;
PFPRECISION *w5,*w3,*wca;
//pfunctionclass *w2,*wmb2;
register int number,maxj;
register PFPRECISION twoscaling,rarray;

 
#ifdef equiout
	ofstream kout;
	kout.open("k.out");
	kout << "sequence length = "<<ct->numofbases<<"\n";
	for (i=1;i<=ct->numofbases;i++) {
		kout << tobase(ct->numseq[i]);
		if (i%20==0) kout << "\n";
	}
	kout << "\n";
	kout.close();
#endif

#ifdef timer
	#include <time.h>
	ofstream timeout;
	int seconds;
	char timerstring[100];
	char timelength[10];
	strcpy(timerstring,"time_pf_");
	sprintf(timelength,"%i",ct->numofbases);
	strcat(timerstring,timelength);
	strcat(timerstring,".out");

	timeout.open(timerstring);
	timeout<<time(NULL)<<"\n";
	seconds = time(NULL);
#endif



//array *v,*w;
//inc is an array that saves time by showing which bases can pair before
//	erg is called

//number is the number of bases being folded
//v[i][j] is the best energy for subsequence i to j when i and j are paired
//	for i<j<n, it is the interior framgment between nucleotides i and j
//	for i<n<j, it is the exterior structure fragmnet that is 5' to j-n and 3' to i
//w[i][j] is the best energy for subsequence i to j

number = (ct->numofbases);//place the number of bases in a registered integer

//scaling is a per nucleotide scale factor for which W and V are divided
//This is necessary to keep the doubles from overflowing:

//scaling = 0.2;//this factor assumes about 1 kcal/mol/base 
twoscaling = data->scaling*data->scaling;


//allocate space for the v and w arrays:
pfunctionclass w(number);
pfunctionclass v(number);
pfunctionclass wmb(number);
pfunctionclass wl(number);
pfunctionclass wmbl(number);
pfunctionclass wcoax(number);
forceclass fce(number);

if (ct->intermolecular) {
	//take advantage of templating to prevent intramolecular base pairs

	ct->allocatetem();
	for (i=1;i<ct->inter[0];i++) {
		for (j=i+1;j<=ct->inter[2];j++) {

			ct->tem[j][i]=false;

		}
	}
	for (i=ct->inter[2]+1;i<ct->numofbases;i++) {
		for (j=i+1;j<=ct->numofbases;j++) {

			ct->tem[j][i]=false;

		}
	}


}


//add a second array for intermolecular folding:

/*if (ct->intermolecular) {
	w2 = new pfunctionclass(number);
	wmb2 = new pfunctionclass(number);
	
   
	

}

else {

	wmb2=NULL;
	w2=NULL;

}*/


   	

lfce = new bool [2*number+1];
mod = new bool [2*number+1];

for (i=0;i<=2*number;i++) {
	lfce[i] = false;
	mod[i] = false;
}


for (i=1;i<=ct->nmod;i++) {

	if (ct->mod[i]!=1&&ct->mod[i]!=ct->numofbases) {
		mod[ct->mod[i]]=true;
		mod[ct->mod[i]+ct->numofbases]=true;
	}
}




w5 = new PFPRECISION [number+1];
w3 = new PFPRECISION [number+2];
wca = new PFPRECISION [number+1];
   
   

   
w5[0] = 1;
w3[number+1] = 1;
   
for (i=0;i<=number;i++) {
	wca[i] = 0;
}
   



force(ct,&fce,lfce);


//This is the fill routine:

if (quickQ) maxj = number;
else maxj = 2*number-1;

for (j=1;j<=maxj;j++) {

	if (((j%10)==0)&&update) update->update((100*j)/(maxj));
	if (j<=number) lowi = 1;
	else lowi = j-number+1+minloop;
	for (i=min(j,number);i>=lowi;i--) {

#ifdef pfdebugmode

		if (i==10&&j==22) {
			i=i;
		}
#endif

	


		rarray=0;
	
	

		if (ct->templated) {
			if (i>ct->numofbases) ii = i - ct->numofbases;
			else ii = i;
			if (j>ct->numofbases) jj = j - ct->numofbases;
			else jj = j;
			if (jj<ii) {
				p = jj;
      			jj = ii;
				ii = p;
			}
   			if (!ct->tem[jj][ii]) goto sub2;
		}

		//Compute v[i][j], the minimum energy of the substructure from i to j,
		//inclusive, where i and j are base paired
		if (fce.f(i,j)&SINGLE) {
			//i or j is forced single-stranded
			//v.f(i,j) = 0;
			goto sub2;
		}
		if (fce.f(i,j)&NOPAIR) {
			//i or j is forced into a pair elsewhere
   			//v.f(i,j)= 0;
		
			goto sub2;
		}

		if (j<=(number)) {
			if ((j-i)<=minloop) goto sub3;
		}

		if (inc[ct->numseq[i]][ct->numseq[j]]==0) {
			//v.f(i,j)= 0;
			goto sub2;
		}

   	
		//force u's into gu pairs
		for (ip=0;ip<ct->ngu;ip++) {
			if (ct->gu[ip]==i) {
				if (ct->numseq[j]!=3) {
         			rarray = 0;
         			goto sub2;
				}
			}
      
			else if (ct->gu[ip]==j) {
       			if (ct->numseq[i]!=3) {
         			rarray = 0;
         			goto sub2;
				}
			}
			else if ((ct->gu[ip]+number)==j) {
       			if (ct->numseq[i]!=3) {
          			rarray = 0;
					goto sub2;
				}
			}

		}


		//now check to make sure that this isn't an isolated pair:
		//	(consider a pair separated by a bulge as not! stacked)

		//before = 0 if a stacked pair cannot form 5' to i
		before =0;
		if ((i>1&&j<(2*number)&&j!=number)) {
			if ((j>number&&((i-j+number)>minloop+2))||j<number) {
				before = inc[ct->numseq[i-1]][ct->numseq[j+1]];
			}
		}
	

		//after = 0 if a stacked pair cannot form 3' to i
		if ((((j-i)>minloop+2)&&(j<=number)||(j>number+1))&&(i!=number)) {
			after = inc[ct->numseq[i+1]][ct->numseq[j-1]];

		}
		else after = 0;

		//if there are no stackable pairs to i.j then don't allow a pair i,j
		if ((before==0)&&(after==0)) {
			//v.f(i,j)= 0;
			goto sub2;
		}
	
		if (i==(number)||j==((number)+1)) goto sub1;


   		//Perhaps i and j close a hairpin:
		rarray=erg3(i,j,ct,data,fce.f(i,j));

		if ((j-i-1)>=(minloop+2)||j>(number)) {
      		//Perhaps i,j stacks over i+1,j-1
			if (!mod[i]&&!mod[j])  //make sure this is not a site of chemical modification
				rarray+=erg1(i,j,i+1,j-1,ct,data)*v.f(i+1,j-1);
			else {
				//allow G-U to be modified or a pair next to a G-U to be modified
				if ((ct->numseq[i]==3&&ct->numseq[j]==4)||(ct->numseq[i]==4&&ct->numseq[j]==3)) {
					rarray+=erg1(i,j,i+1,j-1,ct,data)*v.f(i+1,j-1);	

				}
				else if ((ct->numseq[i+1]==3&&ct->numseq[j-1]==4)||(ct->numseq[i+1]==4&&ct->numseq[j-1]==3)) {

					rarray+=erg1(i,j,i+1,j-1,ct,data)*v.f(i+1,j-1);

				}
				else if (i-1>0&&j+1<2*number) {
					if ((ct->numseq[i-1]==3&&ct->numseq[j+1]==4)||(ct->numseq[i-1]==4&&ct->numseq[j+1]==3)) {

						rarray+=erg1(i,j,i+1,j-1,ct,data)*v.f(i+1,j-1);
				
					}

				}

			}
		}
		
		//Perhaps i,j closes an interior or bulge loop, enumerate all possibilities
      		//possibility
		if (((j-i-1)>=(minloop+3))||(j>(number))) {
      		for (d=(j-i-3);d>=1;d--) {
         		for (ip=(i+1);ip<=(j-1-d);ip++) {
            		jp = d+ip;
					if ((j-i-2-d)>(data->maxintloopsize)) goto sub1;
					if (abs(ip-i+jp-j)<=(data->maxintloopsize)) {
               			if (ip>(number)) {
                  		//if (jp<=number) {

                  		//	v.f(i,j)=min(v.f(i,j),(erg2(i,j,ip,jp,ct,data,fce[i][ip-number],
						//		fce[jp][j])+
						//		v.f(ip-(number),jp)));

						//}
						//else {
                     		rarray+=erg2(i,j,ip,jp,ct,data,fce.f(i,ip-number),fce.f(jp-number,j-number))*
                     			v.f(ip-(number),jp-(number));

							if ((mod[ip]||mod[jp])) if (inc[ct->numseq[ip]][ct->numseq[jp]]&&notgu(ip-(number),jp-(number),ct)&&!(fce.f(ip,jp)&SINGLE)) {
								//ip or jp is modified

								rarray+=erg2(i,j,ip,jp,ct,data,fce.f(i,ip-number),fce.f(jp-number,j-number))*
                     				v.f(ip-(number)+1,jp-(number)-1)*erg1(ip-number,jp-number,ip+1-number,jp-1-number,ct,data);

							}
						//}
						}
						else {
							if (jp<=number) {


                  				rarray+=erg2(i,j,ip,jp,ct,data,fce.f(i,ip),fce.f(jp,j))*v.f(ip,jp);

								if ((mod[ip]||mod[jp])) if (inc[ct->numseq[ip]][ct->numseq[jp]]&&notgu(ip,jp,ct)&&!(fce.f(ip,jp)&SINGLE)) {
									//i or j is modified
									rarray+=erg2(i,j,ip,jp,ct,data,fce.f(i,ip),fce.f(jp,j))*
                  						v.f(ip+1,jp-1)*erg1(ip,jp,ip+1,jp-1,ct,data);

								}
						


						
							}
							else {

						 
                  				rarray+=erg2(i,j,ip,jp,ct,data,fce.f(i,ip),fce.f(jp-number,j-number))*
                  					v.f(ip,jp);


								if ((mod[ip]||mod[jp])) if (inc[ct->numseq[ip]][ct->numseq[jp]]&&notgu(ip,jp,ct)
									&&!(fce.f(ip,jp)&SINGLE)) {
									//i or j is modified
									rarray+=erg2(i,j,ip,jp,ct,data,fce.f(i,ip),fce.f(jp-number,j-number))*
                  						v.f(ip+1,jp-1)*erg1(ip,jp,ip+1,jp-1,ct,data);

								}

							}




						}
					}


			   
				}
			}
		}
		
		//Perhaps i,j closes a multibranch or exterior loop, enumerate all possibilities

		
		sub1:
		
	  
	
	  
		if (((j-i-1)>=(2*minloop+4))||(j>(number))) {
     		

			//consider the exterior loop closed by i,j
			if (j>number) {
         		rarray+= w3[i+1]*w5[j-number-1]*penalty(i,j,ct,data)*twoscaling;
            

				if (i!=number) rarray+= erg4(i,j,i+1,1,ct,data,lfce[i+1])*penalty(i,j,ct,data)*w3[i+2]*w5[j-number-1]*twoscaling;
				if (j!=(number+1)) rarray+= erg4(i,j,j-1,2,ct,data,lfce[j-1])*penalty(i,j,ct,data)*w3[i+1]*w5[j-number-2]*twoscaling;
				if ((i!=number)&&(j!=(number+1))) {
            		rarray+= data->tstack[ct->numseq[i]][ct->numseq[j]][ct->numseq[i+1]][ct->numseq[j-1]]*pfchecknp(lfce[i+1],lfce[j-1])*w3[i+2]*
						w5[j-number-2]*penalty(i,j,ct,data)*twoscaling;

				}
			
			
				//consider the coaxial stacking of a helix from i to ip onto helix ip+1 or ip+2 to j:
				#ifndef disablecoax //a flag that can turn of coaxial stacking
				//first consider a helix stacking from the 5' sequence fragment:
				for (ip=j-number-minloop-1;ip>0;ip--) {
					//first consider flush stacking
					rarray+=
						w3[i+1]*w5[ip-1]*penalty(i,j,ct,data)*penalty(j-number-1,ip,ct,data)*
						ergcoaxflushbases(ip,j-number-1,j-number,i,ct,data)*v.f(ip,j-number-1)*twoscaling;

					if ((mod[ip]||mod[j-number-1])) if (j-number-2>0&&notgu(ip,j-number-1,ct)&&!(fce.f(ip,j-number-1)&SINGLE)) {
						if (inc[ct->numseq[ip+1]][ct->numseq[j-number-2]]) {
							rarray+=
								w3[i+1]*w5[ip-1]*penalty(i,j,ct,data)*penalty(j-number-1,ip,ct,data)*
								ergcoaxflushbases(ip,j-number-1,j-number,i,ct,data)*v.f(ip+1,j-number-2)*
								erg1(ip,j-number-1,ip+1,j-number-2,ct,data)*twoscaling;
						}

					}


					if (j-number-2>0) {
						//now consider an intervening nuc
						if(i<number) {
							rarray+=
								w3[i+2]*w5[ip-1]*penalty(i,j,ct,data)*penalty(ip,j-number-2,ct,data)*
								ergcoaxinterbases2(ip,j-number-2,j-number,i,ct,data)*v.f(ip,j-number-2)*twoscaling
								*pfchecknp(lfce[j-number-1],lfce[i+1]);

				
							if ((mod[ip]||mod[j-number-2])) if (inc[ct->numseq[ip+1]][ct->numseq[j-number-3]]&&notgu(ip,j-number-2,ct)
								&&!(fce.f(ip,j-number-2)&SINGLE)) {
								rarray+=
									w3[i+2]*w5[ip-1]*penalty(i,j,ct,data)*penalty(ip,j-number-2,ct,data)*
									ergcoaxinterbases2(ip,j-number-2,j-number,i,ct,data)*v.f(ip+1,j-number-3)*
									erg1(ip,j-number-2,ip+1,j-number-3,ct,data)*twoscaling*pfchecknp(lfce[j-number-1],lfce[i+1]);


							}
						}


						//consider the other possibility for an intervening nuc
						rarray+=
							w3[i+1]*w5[ip-1]*penalty(i,j,ct,data)*penalty(ip+1,j-number-2,ct,data)*
							ergcoaxinterbases1(ip+1,j-number-2,j-number,i,ct,data)*v.f(ip+1,j-number-2)*twoscaling
							*pfchecknp(lfce[j-number-1],lfce[ip]);


						if ((mod[ip+1]||mod[j-number-2])) if (inc[ct->numseq[ip+2]][ct->numseq[j-number-3]]&&notgu(ip+1,j-number-2,ct)
							&&!(fce.f(ip+1,j-number-2)&SINGLE)) {
							rarray+=
								w3[i+1]*w5[ip-1]*penalty(i,j,ct,data)*penalty(ip+1,j-number-2,ct,data)*
								ergcoaxinterbases1(ip+1,j-number-2,j-number,i,ct,data)*v.f(ip+2,j-number-3)
								*erg1(ip+1,j-number-2,ip+2,j-number-3,ct,data)*twoscaling
								*pfchecknp(lfce[j-number-1],lfce[ip]);
						}


					}

		
				}

				//now consider a helix stacking from the 3' sequence fragment:
				for (ip=i+minloop+1;ip<=number;ip++) {
					//first consider flush stacking
				
					rarray+=
						w3[ip+1]*w5[j-number-1]*penalty(i,j,ct,data)*penalty(ip,i+1,ct,data)*
						ergcoaxflushbases(j-number,i,i+1,ip,ct,data)*v.f(i+1,ip)*twoscaling;


					if ((mod[i+1]||mod[ip])) if (inc[ct->numseq[i+2]][ct->numseq[ip-1]]&&notgu(i+1,ip,ct)
						&&!(fce.f(i+1,ip)&SINGLE)) {

						rarray+=
							w3[ip+1]*w5[j-number-1]*penalty(i,j,ct,data)*penalty(ip,i+1,ct,data)*
							ergcoaxflushbases(j-number,i,i+1,ip,ct,data)*v.f(i+2,ip-1)
							*erg1(i+1,ip,i+2,ip-1,ct,data)*twoscaling;

					}
				
					//now consider an intervening nuc
					if (j-number>1) {
						rarray+=
							w3[ip+1]*w5[j-number-2]*penalty(i,j,ct,data)*penalty(ip,i+2,ct,data)*
							ergcoaxinterbases1(j-number,i,i+2,ip,ct,data)*v.f(i+2,ip)*twoscaling
							*pfchecknp(lfce[i+1],lfce[j-number-1]);


						if ((mod[i+2]||mod[ip])) if (inc[ct->numseq[i+3]][ct->numseq[ip-1]]&&notgu(i+2,ip,ct)
							&&!(fce.f(i+2,ip)&SINGLE)) {

							rarray+=
								w3[ip+1]*w5[j-number-2]*penalty(i,j,ct,data)*penalty(ip,i+2,ct,data)*
								ergcoaxinterbases1(j-number,i,i+2,ip,ct,data)*v.f(i+3,ip-1)
								*erg1(i+2,ip,i+3,ip-1,ct,data)*twoscaling
								*pfchecknp(lfce[i+1],lfce[j-number-1]);

						}
					}


					//consider the other possibility for an intervening nuc
					rarray+=
						w3[ip+1]*w5[j-number-1]*penalty(i,j,ct,data)*penalty(ip-1,i+2,ct,data)*
						ergcoaxinterbases2(j-number,i,i+2,ip-1,ct,data)*v.f(i+2,ip-1)*twoscaling
						*pfchecknp(lfce[i+1],lfce[ip]);

					if ((mod[i+2]||mod[ip-1])) if(inc[ct->numseq[i+3]][ct->numseq[ip-2]]&&notgu(i+2,ip-1,ct)
						&&!(fce.f(i+2,ip-1)&SINGLE)) {

						rarray+=
							w3[ip+1]*w5[j-number-1]*penalty(i,j,ct,data)*penalty(ip-1,i+2,ct,data)*
							ergcoaxinterbases2(j-number,i,i+2,ip-1,ct,data)*v.f(i+3,ip-2)
							*erg1(i+2,ip-1,i+3,ip-2,ct,data)*twoscaling*pfchecknp(lfce[i+1],lfce[ip]);
	

					}
				


				}
				#endif //ifndef disablecoax


			}



		


			//consider the multiloop closed by i,j
			if ((j-i)>(2*minloop+4)&&i!=number) {
          		//no dangling ends on i-j pair:
				if (j-1!=number) {
					rarray+=wmb.f(i+1,j-1)*data->eparam[5]*data->eparam[10]
            			*penalty(i,j,ct,data)*twoscaling;


					//i+1 dangles on i-j pair:
			
					/*if (i!=number)*/ rarray+=erg4(i,j,i+1,1,ct,data,lfce[i+1])*penalty(i,j,ct,data)*
            			wmb.f(i+2,j-1)* data->eparam[5] * data->eparam[6] * data->eparam[10]*twoscaling;
				}
				if (j-2!=number) {
					//j-1 dangles
					if (j!=(number+1))rarray+=erg4(i,j,j-1,2,ct,data,lfce[j-1]) * penalty(i,j,ct,data) *
            			wmb.f(i+1,j-2) * data->eparam[5] * data->eparam[6] * data->eparam[10]*twoscaling;
					//both i+1 and j-1 dangle
					if (/*(i!=number)&&*/(j!=(number+1))) {
            			rarray+=data->tstkm[ct->numseq[i]][ct->numseq[j]][ct->numseq[i+1]][ct->numseq[j-1]]*
									pfchecknp(lfce[i+1],lfce[j-1])*
									wmb.f(i+2,j-2) * data->eparam[5] * data->eparam[6] * data->eparam[6]* data->eparam[10]
									*penalty(i,j,ct,data)*twoscaling;
					}
				}

		
			

		
				//consider the coaxial stacking of a helix from i to j onto helix i+1 or i+2 to ip:
				#ifndef disablecoax //a flag to turn off coaxial stacking
				for (ip=i+1;(ip<j);ip++) {
					//first consider flush stacking


				

					//conditions guarantee that the coaxial stacking isn't considering an exterior loop 
					//if ((i!=number)/*&&(i+1!=number)*//*&&((j>number)||(ip!=number)&&(ip+1!=number))&&(j-1!=number)*/) {
					if (i!=number&&ip!=number&&j-1!=number) {
						rarray+=penalty(i,j,ct,data)*v.f(i+1,ip)*
							penalty(i+1,ip,ct,data)*data->eparam[5]
							*data->eparam[10]*data->eparam[10]*(w.f(ip+1,j-1)+wmb.f(ip+1,j-1))*ergcoaxflushbases(j,i,i+1,ip,ct,data)
							*twoscaling;

			
						if((mod[i+1]||mod[ip])) if (inc[ct->numseq[i+2]][ct->numseq[ip-1]]&&notgu(i+1,ip,ct)&&!(fce.f(i+1,ip)&SINGLE)) {

							rarray+=penalty(i,j,ct,data)*v.f(i+2,ip-1)*
								penalty(i+1,ip,ct,data)*data->eparam[5]
								*data->eparam[10]*data->eparam[10]*(w.f(ip+1,j-1)+wmb.f(ip+1,j-1))*ergcoaxflushbases(j,i,i+1,ip,ct,data)
								*erg1(i+1,ip,i+2,ip-1,ct,data)*twoscaling;

						}
				


						//if ((ip<j-1)&&(i+2!=number)) {
						if (ip+2<j-1&&i+1!=number&&ip+1!=number) {
						//now consider an intervening nuc
							if ((ip+2<j-1)/*&&(j>number||ip+2!=number)*/)
							rarray+=penalty(i,j,ct,data)*v.f(i+2,ip)*
								penalty(i+2,ip,ct,data)*data->eparam[5]
								*data->eparam[6]*data->eparam[6]*data->eparam[10]*data->eparam[10]*
								(w.f(ip+2,j-1)+wmb.f(ip+2,j-1))
								*ergcoaxinterbases2(j,i,i+2,ip,ct,data)*twoscaling*pfchecknp(lfce[i+1],lfce[ip+1]);

							if((mod[i+2]||mod[ip])) if (inc[ct->numseq[i+3]][ct->numseq[ip-1]]&&notgu(i+2,ip,ct)
								&&!(fce.f(i+2,ip)&SINGLE)) {

								rarray+=penalty(i,j,ct,data)*v.f(i+3,ip-1)*
									penalty(i+2,ip,ct,data)*data->eparam[5]
									*data->eparam[6]*data->eparam[6]*data->eparam[10]*data->eparam[10]*
									(w.f(ip+2,j-1)+wmb.f(ip+2,j-1))
									*ergcoaxinterbases2(j,i,i+2,ip,ct,data)
									*erg1(i+2,ip,i+3,ip-1,ct,data)*twoscaling*pfchecknp(lfce[i+1],lfce[ip+1]);

							}
		


							if (ip+1<j-2&&j-2!=number)
							rarray+=penalty(i,j,ct,data)*v.f(i+2,ip)*
								penalty(i+2,ip,ct,data)*data->eparam[5]
								*data->eparam[6]*data->eparam[6]*data->eparam[10]*data->eparam[10]*
								(w.f(ip+1,j-2)+wmb.f(ip+1,j-2))
								*ergcoaxinterbases1(j,i,i+2,ip,ct,data)*twoscaling
								*pfchecknp(lfce[i+1],lfce[j-1]);

							if((mod[i+2]||mod[ip])) if (inc[ct->numseq[i+3]][ct->numseq[ip-1]]&&notgu(i+2,ip,ct)
								&&!(fce.f(i+2,ip)&SINGLE)) {

								rarray+=penalty(i,j,ct,data)*v.f(i+3,ip-1)*
									penalty(i+2,ip,ct,data)*data->eparam[5]
									*data->eparam[6]*data->eparam[6]*data->eparam[10]*data->eparam[10]*
									(w.f(ip+1,j-2)+wmb.f(ip+1,j-2))
									*ergcoaxinterbases1(j,i,i+2,ip,ct,data)
									*erg1(i+2,ip,i+3,ip-1,ct,data)*twoscaling*pfchecknp(lfce[i+1],lfce[j-1]);

							}

				


						}




			
					}


				}

			

				//consider the coaxial stacking of a helix from i to j onto helix ip to j-2 or j-1:
				for (ip=j-1;ip>i;ip--) {


					//conditions guarantee that the coaxial stacking isn't considering an exterior loop
					//if ((i!=number)&&(i+1!=number)&&((j>number)||(ip!=number)&&(ip-1!=number))&&(j-1!=number)) {
					if (j-1!=number&&ip-1!=number&&i!=number) {
						//first consider flush stacking
						rarray+=penalty(i,j,ct,data)*v.f(ip,j-1)*
							penalty(j-1,ip,ct,data)*data->eparam[5]
							*data->eparam[10]*data->eparam[10]*(w.f(i+1,ip-1)+wmb.f(i+1,ip-1))*ergcoaxflushbases(ip,j-1,j,i,ct,data)
							*twoscaling;

					
						if((mod[ip]||mod[j-1])) if(inc[ct->numseq[ip+1]][ct->numseq[j-2]]&&notgu(ip,j-1,ct)&&!(fce.f(ip,j-1)&SINGLE)) {
							rarray+=penalty(i,j,ct,data)*v.f(ip+1,j-2)*
								penalty(j-1,ip,ct,data)*data->eparam[5]
								*data->eparam[10]*data->eparam[10]*(w.f(i+1,ip-1)+wmb.f(i+1,ip-1))
								*ergcoaxflushbases(ip,j-1,j,i,ct,data)
								*erg1(ip,j-1,ip+1,j-2,ct,data)*twoscaling;

						}

				

		

						if (j-2!=number) {
							//now consider an intervening nuc
							//if ((ip>i+1)&&(j>number||ip-2!=number))
							if (ip-2>i+1&&ip-2!=number) {
								rarray+=penalty(i,j,ct,data)*v.f(ip,j-2)*
									penalty(j-2,ip,ct,data)*data->eparam[5]
									*data->eparam[6]*data->eparam[6]*data->eparam[10]*data->eparam[10]*
									(w.f(i+1,ip-2)+wmb.f(i+1,ip-2))
									*ergcoaxinterbases1(ip,j-2,j,i,ct,data)*twoscaling*pfchecknp(lfce[j-1],lfce[ip-1]);
				


								if((mod[ip]||mod[j-2])) if(inc[ct->numseq[ip+1]][ct->numseq[j-3]]&&notgu(ip,j-2,ct)&&!(fce.f(ip,j-2)&SINGLE)) {
									rarray+=penalty(i,j,ct,data)*v.f(ip+1,j-3)*
										penalty(j-2,ip,ct,data)*data->eparam[5]
										*data->eparam[6]*data->eparam[6]*data->eparam[10]*data->eparam[10]*
										(w.f(i+1,ip-2)+wmb.f(i+1,ip-2))
										*ergcoaxinterbases1(ip,j-2,j,i,ct,data)
										*erg1(ip,j-2,ip+1,j-3,ct,data)*twoscaling*pfchecknp(lfce[j-1],lfce[ip-1]);

								}
							}



							if ((ip-1>i+2)&&i+1!=number) {
								rarray+=penalty(i,j,ct,data)*v.f(ip,j-2)*
									penalty(j-2,ip,ct,data)*data->eparam[5]
									*data->eparam[6]*data->eparam[6]*data->eparam[10]*data->eparam[10]*
									(w.f(i+2,ip-1)+wmb.f(i+2,ip-1))
									*ergcoaxinterbases2(ip,j-2,j,i,ct,data)*twoscaling*pfchecknp(lfce[j-1],lfce[i+1]);

								if((mod[ip]||mod[j-2])) if(inc[ct->numseq[ip+1]][ct->numseq[j-3]]&&notgu(ip,j-2,ct)
									&&!(fce.f(ip,j-2)&SINGLE)) {
									rarray+=penalty(i,j,ct,data)*v.f(ip+1,j-3)*
										penalty(j-2,ip,ct,data)*data->eparam[5]
										*data->eparam[6]*data->eparam[6]*data->eparam[10]*data->eparam[10]*
										(w.f(i+2,ip-1)+wmb.f(i+2,ip-1))
										*ergcoaxinterbases2(ip,j-2,j,i,ct,data)
										*erg1(ip,j-2,ip+1,j-3,ct,data)*twoscaling*pfchecknp(lfce[j-1],lfce[i+1]);

								}
							}

						
						}
					

				
					}

				
				
				
				}
				#endif //ifndef disablecoax


				/*if (ct->intermolecular) {

            		//intermolecular, so consider wmb2,
					//don't add the multiloop penalties because this is a exterior loop

            		e[1] = min(e[1],wmb2->f(i+1,j-1) + penalty(i,j,ct,data)+infinity);


            		//i+1 dangles on i-j pair:
            		if (i!=number) e[2] = min(e[2],erg4(i,j,i+1,1,ct,data,lfce[i+1]) + penalty(i,j,ct,data) +
            			wmb2->f(i+2,j-1)+infinity);
            		//j-1 dangles
            		if (j!=(number+1)) e[3] = min(e[3],erg4(i,j,j-1,2,ct,data,lfce[j-1]) + penalty(i,j,ct,data) +
            			wmb2->f(i+1,j-2)+infinity);
            		//both i+1 and j-1 dangle
            		if ((i!=number)&&(j!=(number+1))) {
            			e[4] = min(e[4],
            			data->tstkm[ct->numseq[i]][ct->numseq[j]][ct->numseq[i+1]][ct->numseq[j-1]] +
							pfchecknp(lfce[i+1],lfce[j-1]) +
               				wmb2->f(i+2,j-2) + penalty(i,j,ct,data)+infinity);

					}

				


				}*/


			}


      
		}

      
		sub2:

		v.f(i,j) = rarray;

		if (fce.f(i,j)&PAIR)  {//force a pair between i and j
	  		w.f(i,j) = v.f(i,j)*data->eparam[10]*penalty(i,j,ct,data);
			wl.f(i,j) = v.f(i,j)*data->eparam[10]*penalty(i,j,ct,data);
	  		goto sub3;
		}
	  
		////fill wmb:
		rarray = 0;
		if (((j-i-1)>(2*minloop+2))||j>number) {
		         
		
			#ifdef pfdebugmode
				ofstream dump;
				if (i==10&&j==22) {	
						
					dump.open("dump.out");
						 
				}
			#endif
         


			//also consider the coaxial stacking of two helixes in wv
			e = 0;
			#ifndef disablecoax //a flag to diable coaxial stacking
			for (ip=i+minloop+1;ip<j-minloop-1;ip++) {
				//first consider flush stacking
		
				if (ip!=number) {
					rarray+=v.f(i,ip)*v.f(ip+1,j)*penalty(i,ip,ct,data)
						*penalty(ip+1,j,ct,data)*ergcoaxflushbases(i,ip,ip+1,j,ct,data);
					
					#ifdef pfdebugmode
				
					if (i==10&&j==22) {	
						
						dump<< "i= "<<i<<" j= "<<j<<" ip= "<<ip<<" 11.d.v.1.d.ii.1 rarray+= "<<v.f(i,ip)*v.f(ip+1,j)*penalty(i,ip,ct,data)
						*penalty(ip+1,j,ct,data)*ergcoaxflushbases(i,ip,ip+1,j,ct,data)<<" rarray = "<<rarray<<"\n";

						 
					}
					#endif

					if ((mod[i]||mod[ip]||mod[ip+1]||mod[j])) {

						if ((mod[i]||mod[ip])&&(mod[ip+1]||mod[j])&&inc[ct->numseq[i+1]][ct->numseq[ip-1]]
							&&inc[ct->numseq[ip+2]][ct->numseq[j-1]]&&notgu(i,ip,ct)&&notgu(ip+1,j,ct)
								&&!(fce.f(ip+1,j)&SINGLE)&&!(fce.f(i,ip)&SINGLE)) {

							rarray+=v.f(i+1,ip-1)*v.f(ip+2,j-1)*penalty(i,ip,ct,data)
								*penalty(ip+1,j,ct,data)*ergcoaxflushbases(i,ip,ip+1,j,ct,data)
								*erg1(i,ip,i+1,ip-1,ct,data)*erg1(ip+1,j,ip+2,j-1,ct,data);


						}

						if ((mod[i]||mod[ip])&&inc[ct->numseq[i+1]][ct->numseq[ip-1]]&&notgu(i,ip,ct)&&!(fce.f(i,ip)&SINGLE)) {
							
							rarray+=v.f(i+1,ip-1)*v.f(ip+1,j)*penalty(i,ip,ct,data)
								*penalty(ip+1,j,ct,data)*ergcoaxflushbases(i,ip,ip+1,j,ct,data)
								*erg1(i,ip,i+1,ip-1,ct,data);


						}

						if ((mod[ip+1]||mod[j])&&inc[ct->numseq[ip+2]][ct->numseq[j-1]]&&notgu(ip+1,j,ct)&&!(fce.f(ip+1,j)&SINGLE)) {


							rarray+=v.f(i,ip)*v.f(ip+2,j-1)*penalty(i,ip,ct,data)
								*penalty(ip+1,j,ct,data)*ergcoaxflushbases(i,ip,ip+1,j,ct,data)
								*erg1(ip+1,j,ip+2,j-1,ct,data);

						}


					}
				


					if (ip+1!=number&&j!=number+1) {
						if (!lfce[ip+1]&&!lfce[j]) {
							//now consider an intervening mismatch
							e+=v.f(i,ip)*v.f(ip+2,j-1)*penalty(i,ip,ct,data)
								*penalty(ip+2,j-1,ct,data)*ergcoaxinterbases2(i,ip,ip+2,j-1,ct,data);

					#ifdef pfdebugmode
				
					if (i==10&&j==22) {	
						
						dump<< "i= "<<i<<" j= "<<j<<" ip= "<<ip<<" 11.d.v.1.d.ii.2 earray+= "<<v.f(i,ip)*v.f(ip+2,j-1)*penalty(i,ip,ct,data)
								*penalty(ip+2,j-1,ct,data)*ergcoaxinterbases2(i,ip,ip+2,j-1,ct,data)
								<<" earray = "<<e<<"\n";

						 
					}
					#endif

							if (mod[i]||mod[ip]||mod[ip+2]||mod[j-1]) {
								if ((mod[i]||mod[ip])&&(mod[ip+2]||mod[j-1])&&inc[ct->numseq[i+1]][ct->numseq[ip-1]]
									&&inc[ct->numseq[ip+3]][ct->numseq[j-2]]&&notgu(i,ip,ct)&&notgu(ip+2,j-1,ct)
										&&!(fce.f(i,ip)&SINGLE)&&!(fce.f(ip+2,j-1)&SINGLE)) {

									 e+=v.f(i+1,ip-1)*v.f(ip+3,j-2)*penalty(i,ip,ct,data)
										*penalty(ip+2,j-1,ct,data)*ergcoaxinterbases2(i,ip,ip+2,j-1,ct,data)
										*erg1(i,ip,i+1,ip-1,ct,data)*erg1(ip+2,j-1,ip+3,j-2,ct,data);


								}

								if ((mod[i]||mod[ip])&&inc[ct->numseq[i+1]][ct->numseq[ip-1]]&&notgu(i,ip,ct)&&!(fce.f(i,ip)&SINGLE)) {
							
									e+=v.f(i+1,ip-1)*v.f(ip+2,j-1)*penalty(i,ip,ct,data)
										*penalty(ip+2,j-1,ct,data)*ergcoaxinterbases2(i,ip,ip+2,j-1,ct,data)
										*erg1(i,ip,i+1,ip-1,ct,data);


								}

								if ((mod[ip+2]||mod[j-1])&&inc[ct->numseq[ip+3]][ct->numseq[j-2]]&&notgu(ip+2,j-1,ct)&&!(fce.f(ip+2,j-1)&SINGLE)) {


									e+=v.f(i,ip)*v.f(ip+3,j-2)*penalty(i,ip,ct,data)
										*penalty(ip+2,j-1,ct,data)*ergcoaxinterbases2(i,ip,ip+2,j-1,ct,data)
										*erg1(ip+2,j-1,ip+3,j-2,ct,data);

								}
							}
						}
		
						if(!lfce[i]&&!lfce[ip+1]&&i!=number) {
							e+=v.f(i+1,ip)*v.f(ip+2,j)*penalty(i+1,ip,ct,data)
								*penalty(ip+2,j,ct,data)*ergcoaxinterbases1(i+1,ip,ip+2,j,ct,data);

					#ifdef pfdebugmode
				
					if (i==10&&j==22) {	
						
						dump<< "i= "<<i<<" j= "<<j<<" ip= "<<ip<<" 11.d.v.1.d.ii.3 earray+= "<<v.f(i+1,ip)*v.f(ip+2,j)*penalty(i+1,ip,ct,data)
								*penalty(ip+2,j,ct,data)*ergcoaxinterbases1(i+1,ip,ip+2,j,ct,data)
								<<" earray = "<<e<<"\n";

						 
					}
					#endif
							if (mod[i+1]||mod[ip]||mod[ip+2]||mod[j]) {
								if ((mod[i+1]||mod[ip])&&(mod[ip+2]||mod[j])&&inc[ct->numseq[i+2]][ct->numseq[ip-1]]
									&&inc[ct->numseq[ip+3]][ct->numseq[j-1]]&&notgu(i+1,ip,ct)&&notgu(ip+2,j,ct)
									&&!(fce.f(i+1,ip)&SINGLE)&&!(fce.f(ip+2,j)&SINGLE)	) {

									e+=v.f(i+2,ip-1)*v.f(ip+3,j-1)*penalty(i+1,ip,ct,data)
										*penalty(ip+2,j,ct,data)*ergcoaxinterbases1(i+1,ip,ip+2,j,ct,data)
										*erg1(i+1,ip,i+2,ip-1,ct,data)*erg1(ip+2,j,ip+3,j-1,ct,data);	


							
								}
								if ((mod[i+1]||mod[ip])&&inc[ct->numseq[i+2]][ct->numseq[ip-1]]&&notgu(i+1,ip,ct)&&!(fce.f(i+1,ip)&SINGLE)) {
							
									e+=v.f(i+2,ip-1)*v.f(ip+2,j)*penalty(i+1,ip,ct,data)
										*penalty(ip+2,j,ct,data)*ergcoaxinterbases1(i+1,ip,ip+2,j,ct,data)
										*erg1(i+1,ip,i+2,ip-1,ct,data);


								}

								if ((mod[ip+2]||mod[j])&&inc[ct->numseq[ip+3]][ct->numseq[j-1]]&&notgu(ip+2,j,ct)&&!(fce.f(ip+2,j)&SINGLE)) {


									e+=v.f(i+1,ip)*v.f(ip+3,j-1)*penalty(i+1,ip,ct,data)
										*penalty(ip+2,j,ct,data)*ergcoaxinterbases1(i+1,ip,ip+2,j,ct,data)
										*erg1(ip+2,j,ip+3,j-1,ct,data);

								}
							}
						}
					}
				}

			



			}
			#endif //ifndef disablecoax
			#ifdef pfdebugmode
				
				if (i==10&&j==22) {	
						
					dump<< "i= "<<i<<" j= "<<j<<" ip= "<<ip<<" right before 11.d.v.1.d.iii execution rarray= "<<rarray<<" earray = "<<e<<"\n";
					dump.close();
						 
				}
			#endif
			
			wca[i] = rarray+e;
			rarray =(rarray+e*data->eparam[6]*data->eparam[6])*data->eparam[10]*data->eparam[10];
			wcoax.f(i,j) = rarray;

			//search for an open bifurcation:
			for (k=i+1;k<j;k++) {
				//e = 0;
				if (k!=number) {
					if (!lfce[i]&&i!=number)
						rarray+=(wl.f(i,k)-wl.f(i+1,k)*data->eparam[6]*data->scaling+wcoax.f(i,k))*(wl.f(k+1,j)+wmbl.f(k+1,j));

					else rarray+=(wl.f(i,k)+wcoax.f(i,k))*(wl.f(k+1,j)+wmbl.f(k+1,j));

				}
         	}

			
			

			if (i!=number)
				if (!lfce[i]) rarray+=wmbl.f(i+1,j)*data->eparam[6]*data->scaling;

			wmbl.f(i,j) = rarray;

			wmb.f(i,j) = rarray;
			if (j!=number+1)
				if (!lfce[j]) wmb.f(i,j)+=wmb.f(i,j-1)*data->eparam[6]*data->scaling;	
			

			/*if (ct->intermolecular) {
         		//intermolecular folding:

		
			

         	//search for an open bifurcation:
         	for (k=i;k<=j;k++) {

				if (k!=number) wmb2->f(i,j) = min(wmb2->f(i,j),w2->f(i,k)+work2[k+1][jmt]);

				
         	}

		
			
			if (i!=number)
				if (!(fce.f(i,i)&INTER)) wmb2->f(i,j) = min(wmb2->f(i,j) ,wmb2->f(i+1,j) );
				else  wmb2->f(i,j) = min(wmb2->f(i,j) ,wmb2->f(i+1,j) + data->init - infinity);

			if (j!=number+1)
				if (!(fce.f(j,j)&INTER)) wmb2->f(i,j)  = min(wmb2->f(i,j) ,wmb2->f(i,j-1));
				else wmb2->f(i,j)  = min(wmb2->f(i,j) ,wmb2->f(i,j-1) +data->init-infinity);


         	
			w2->f(i,j) = min(w2->f(i,j),wmb2->f(i,j) );

		}*/
		

		}

		//Compute w[i][j]
		if (j>number||(j-i>minloop)) {
			wl.f(i,j)= data->eparam[10]*v.f(i,j)*penalty(j,i,ct,data);

			if ((mod[i]||mod[j])) if (inc[ct->numseq[i+1]][ct->numseq[j-1]]&&notgu(i,j,ct)) {

				wl.f(i,j)+= data->eparam[10]*v.f(i+1,j-1)*penalty(j,i,ct,data)*erg1(i,j,i+1,j-1,ct,data);			

			} 



			if (i!=number) {
      			//calculate the energy of i stacked onto the pair of i+1,j

				wl.f(i,j)+= v.f(i+1,j)*data->eparam[10]*data->eparam[6]*
         			erg4(j,i+1,i,2,ct,data,lfce[i])*penalty(i+1,j,ct,data);

				if ((mod[i+1]||mod[j])) if(inc[ct->numseq[i+2]][ct->numseq[j-1]]&&notgu(i+1,j,ct)&&!(fce.f(i+1,j)&SINGLE)) {

					wl.f(i,j)+= v.f(i+2,j-1) * data->eparam[10] *data->eparam[6] *
         				erg4(j,i+1,i,2,ct,data,lfce[i])*penalty(i+1,j,ct,data)
						*erg1(i+1,j,i+2,j-1,ct,data);


				}
         
			}
			if (j!=((number)+1)) {
      			//calculate the energy of j stacked onto the pair of i,j-1
				if (j!=1) {
         			wl.f(i,j)+= v.f(i,j-1)* data->eparam[10] * data->eparam[6] *
         				erg4(j-1,i,j,1,ct,data,lfce[j])*penalty(i,j-1,ct,data);

					if ((mod[i]||mod[j-1])) if(inc[ct->numseq[i+1]][ct->numseq[j-2]]&&notgu(i,j-1,ct)&&!(fce.f(i,j-1)&SINGLE)) {

						wl.f(i,j)+= v.f(i+1,j-2) * data->eparam[10] * data->eparam[6] *
         					erg4(j-1,i,j,1,ct,data,lfce[j])*penalty(i,j-1,ct,data)
							*erg1(i,j-1,i+1,j-2,ct,data);
				


					}

         	
				}
			}
			if ((i!=(number))&&(j!=((number)+1))) {
      			//calculate i and j stacked onto the pair of i+1,j-1
				if (j!=1&&!lfce[i]&&!lfce[j]) {
         			wl.f(i,j)+= v.f(i+1,j-1) *data->eparam[10] * (data->eparam[6]*data->eparam[6]) *
         				data->tstkm[ct->numseq[j-1]][ct->numseq[i+1]][ct->numseq[j]][ct->numseq[i]]
						*penalty(j-1,i+1,ct,data);



					if ((mod[i+1]||mod[j-1])) if((j-2>0)&&!(fce.f(i+1,j-1)&SINGLE)) {
						if(inc[ct->numseq[i+2]][ct->numseq[j-2]]&&notgu(i+1,j-1,ct)) {

							wl.f(i,j)+= v.f(i+2,j-2) * data->eparam[10] * (data->eparam[6]*data->eparam[6]) *
         						data->tstkm[ct->numseq[j-1]][ct->numseq[i+1]][ct->numseq[j]][ct->numseq[i]]
								*penalty(j-1,i+1,ct,data)*erg1(i+1,j-1,i+2,j-2,ct,data);

						}
					}
				}
			}

	  
	  

			if (i!=number&&!lfce[i]) {
         		//if (!(fce.f(i,i)&INTER))
               //add a nuc to an existing loop:
         		wl.f(i,j)+=  wl.f(i+1,j)*data->eparam[6]*data->scaling;
            	//this is for when i represents the center of an intermolecular linker:
				// else e[4] = w.f(i+1,j) + data->eparam[6] + infinity;
			}

			w.f(i,j) = wl.f(i,j);
			if (j!=number+1&&!lfce[j]) {
             	//if (!(fce.f(j,j)&INTER)) {
               	//add a nuc to an existing loop:
               	w.f(i,j)+= w.f(i,j-1) * data->eparam[6]*data->scaling;
               //}
               //else e[5] = w.f(i,j-1) + data->eparam[6] + infinity;

			}
		}
		
		/* if (ct->intermolecular) {

      		//wmb2[i][j%3] = infinity;
      		//keep track of w2:
			for (ii=1;ii<=5;ii++) e[ii] = 2*infinity;


			if (i!=number) {
      			//calculate the energy of i stacked onto the pair of i+1,j

         		e[1] = v.f(i+1,j) +
         			erg4(j,i+1,i,2,ct,data,lfce[i])+penalty(i+1,j,ct,data);

				if ((mod[i+1]||mod[j])) if(inc[ct->numseq[i+2]][ct->numseq[j-1]]) {

					e[1] = min(e[1],v.f(i+2,j-1) +
         				erg4(j,i+1,i,2,ct,data,lfce[i])+penalty(i+1,j,ct,data)
						+erg1(i+1,j,i+2,j-1,ct,data));

				}


         		
         		e[4] = w2->f(i+1,j);
            	
			}
      		if (j!=((number)+1)) {
      		//calculate the energy of j stacked onto the pair of i,j-1
         		if (j!=1) {
         			e[2] = v.f(i,j-1)   +
         				erg4(j-1,i,j,1,ct,data,lfce[j])+penalty(i,j-1,ct,data);

					if ((mod[i]||mod[j-1])&&inc[ct->numseq[i+1]][ct->numseq[j-2]]) {

						e[2] = min(e[2],v.f(i+1,j-2) +
         					erg4(j-1,i,j,1,ct,data,lfce[j])+penalty(i,j-1,ct,data)
							+erg1(i,j-1,i+1,j-2,ct,data));

					}

         			
               		e[5] = w2->f(i,j-1);
                
         		}
      		}
      		if ((i!=(number))&&(j!=((number)+1))) {
      			//calculate i and j stacked onto the pair of i+1,j-1
         		if (j!=1) {
         			e[3] = v.f(i+1,j-1)   +
         				data->tstkm[ct->numseq[j-1]][ct->numseq[i+1]]
									[ct->numseq[j]][ct->numseq[i]]
					+pfchecknp(lfce[i+1],lfce[j-1])
               		+penalty(j-1,i+1,ct,data);



					if ((mod[i+1]||mod[j-1])&&inc[ct->numseq[i+2]][ct->numseq[j-2]]) {

						e[3] = min(e[3],v.f(i+2,j-2) +
         					data->tstkm[ct->numseq[j-1]][ct->numseq[i+1]][ct->numseq[j]][ct->numseq[i]]
							+pfchecknp(lfce[i+1],lfce[j-1])
							+penalty(j-1,i+1,ct,data)+erg1(i+1,j-1,i+2,j-2,ct,data));

					}
         		}
      		}

			e[1] = min(e[1],(v.f(i,j)+penalty(j,i,ct,data)));

			if (mod[i]||mod[j]&&inc[ct->numseq[i+1]][ct->numseq[j-1]]) {

				e[1] = min((v.f(i+1,j-1)+penalty(j,i,ct,data)+erg1(i,j,i+1,j-1,ct,data)),e[1]);			

			}


      		w2->f(i,j) = min(e[1],e[2]);
      		w2->f(i,j) = min(w2->f(i,j),e[3]);
      		w2->f(i,j) = min(w2->f(i,j),e[4]);
      		w2->f(i,j) = min(w2->f(i,j),e[5]);




		}*/

		        

      

	


      
		sub3:

      

      
   


		//Compute w5[i], the energy of the best folding from 1->i, and
      		//w3[i], the energy of the best folding from i-->numofbases
		if (i==(1)) {


			if (j<=minloop+1) {
				if (lfce[j]) w5[j]= 0;
				else  w5[j] = w5[j-1]*data->scaling;
			}
				   
			else {
      			if (lfce[j]) rarray = 0;
			
				else rarray = w5[j-1]*data->scaling;
				//if(j==2804)	cout<<"\nw5 rararay(w5[j-1]): "<<rarray<<"\n"; 
      		
			
      			for (k=0;k<=(j-4);k++) {



      				rarray+=w5[k]*v.f(k+1,j)*penalty(j,k+1,ct,data);
					

					if ((mod[k+1]||mod[j])) if(inc[ct->numseq[k+2]][ct->numseq[j-1]]&&notgu(k+1,j,ct)&&!(fce.f(k+1,j)&SINGLE)) {

						rarray+=w5[k]*v.f(k+2,j-1)*penalty(j,k+1,ct,data)
							*erg1(k+1,j,k+2,j-1,ct,data);
					}

			

					rarray+=w5[k]*erg4(j,k+2,k+1,2,ct,data,lfce[k+1])*v.f(k+2,j)*penalty(j,k+2,ct,data);

				
					if((mod[k+2]||mod[j])) if(inc[ct->numseq[k+3]][ct->numseq[j-1]]&&notgu(k+2,j,ct)
						&&!(fce.f(k+2,j)&SINGLE)) {
						rarray+=w5[k]*erg4(j,k+2,k+1,2,ct,data,lfce[k+1])*v.f(k+3,j-1)
							*penalty(j,k+2,ct,data)*erg1(k+2,j,k+3,j-1,ct,data);

					}


         			rarray+=w5[k]*erg4(j-1,k+1,j,1,ct,data,lfce[j])*v.f(k+1,j-1)*penalty(j-1,k+1,ct,data);

			
					if ((mod[k+1]||mod[j-1])) if(inc[ct->numseq[k+2]][ct->numseq[j-2]]&&notgu(k+1,j-1,ct)&&!(fce.f(k+1,j-1)&SINGLE)) {

						rarray+=w5[k]*erg4(j-1,k+1,j,1,ct,data,lfce[j])*v.f(k+2,j-2)
							*penalty(j-1,k+1,ct,data)*erg1(k+1,j-1,k+2,j-2,ct,data);
					}



					rarray+=w5[k]*data->tstack[ct->numseq[j-1]][ct->numseq[k+2]][ct->numseq[j]][ct->numseq[k+1]] 
									*pfchecknp(lfce[j],lfce[k+1]) * v.f(k+2,j-1)*
									penalty(j-1,k+2,ct,data);
				
					if ((mod[k+2]||mod[j-1])) if(inc[ct->numseq[k+3]][ct->numseq[j-2]]&&notgu(k+2,j-1,ct)&&!(fce.f(k+2,j-1)&SINGLE)) {

						rarray+=w5[k]*data->tstack[ct->numseq[j-1]][ct->numseq[k+2]][ct->numseq[j]][ct->numseq[k+1]] 
									*pfchecknp(lfce[j],lfce[k+1]) * v.f(k+3,j-2)*
									penalty(j-1,k+2,ct,data)*erg1(k+2,j-1,k+3,j-2,ct,data);

					}



					
				
		


					rarray+=w5[k]*wca[k+1];
					if(j==2804){
					//	cout<<k <<"w5 rararay: "<<rarray<<"\t"<<"wca[k+1]: "<<wca[k+1]<<"\t";
					//	cout<<" w5[k]: "<<w5[k]<<"\n";
					}
		
	
				}
			
				 
				w5[j] = rarray;
		
				//check to see if w5 is about to go out of bounds:
				if (w5[j]>PFMAX) {
					rescale(i,j,ct,data,&v,&w,&wl,&wcoax,&wmb,&wmbl,w5,w3,wca,SCALEDOWN);
					twoscaling = twoscaling*SCALEDOWN*SCALEDOWN;
				}
				else if (w5[j]<PFMIN&&w5[j]>0) {
					rescale(i,j,ct,data,&v,&w,&wl,&wcoax,&wmb,&wmbl,w5,w3,wca,SCALEUP);
					twoscaling=twoscaling*SCALEUP*SCALEUP;

				}					
			}


			if (j==number) {

				//w3[0] = 0;
				//w3[number+1] = 0;
				for (ii=(number);ii>=(number-minloop);ii--) {    //number+1 ... number-minloop
      				if (lfce[ii]) w3[ii] = 0;
					else w3[ii]=w3[ii+1]*data->scaling;
				}
				//w3[i]=0;
   				for (ii=((number)-minloop-1);ii>=1;ii--) {

      				if (lfce[ii]) rarray = 0;
         
   					else rarray = w3[ii+1]*data->scaling;
         
      	
							
					for (k=((number)+1);k>=(ii+4);k--) {
      					rarray+=v.f(ii,k-1)*w3[k]*penalty(k-1,ii,ct,data);

						if((mod[ii]||mod[k-1])) if(inc[ct->numseq[ii+1]][ct->numseq[k-2]]&&notgu(ii,k-1,ct)&&!(fce.f(ii,k-1)&SINGLE)) {
							rarray+=v.f(ii+1,k-2)*w3[k]*penalty(k-1,ii,ct,data)*erg1(ii,k-1,ii+1,k-2,ct,data);

						}


						rarray+=v.f(ii+1,k-1)*erg4(k-1,ii+1,ii,2,ct,data,lfce[ii])*penalty(k-1,ii+1,ct,data) * w3[k];

						if((mod[ii+1]||mod[k-1])) if(inc[ct->numseq[ii+2]][ct->numseq[k-2]]&&notgu(ii+1,k-1,ct)&&!(fce.f(ii+1,k-1)&SINGLE)) {

							rarray+=v.f(ii+2,k-2)*erg4(k-1,ii+1,ii,2,ct,data,lfce[ii])*
								penalty(k-1,ii+1,ct,data) *w3[k]*erg1(ii+1,k-1,ii+2,k-2,ct,data);

						}


						rarray+=v.f(ii,k-2)*erg4(k-2,ii,k-1,1,ct,data,lfce[k-1])*penalty(k-2,ii,ct,data)*w3[k];

						if((mod[ii]||mod[k-2]))if(inc[ct->numseq[ii+1]][ct->numseq[k-3]]&&notgu(ii,k-2,ct)&&!(fce.f(ii,k-2)&SINGLE)) {
							rarray+=v.f(ii+1,k-3)*erg4(k-2,ii,k-1,1,ct,data,lfce[k-1])* 
								penalty(k-2,ii,ct,data)*w3[k]*erg1(ii,k-2,ii+1,k-3,ct,data);

						}

						if (!lfce[ii]&&!lfce[k-1]) {
							rarray+=v.f(ii+1,k-2)*data->tstack[ct->numseq[k-2]][ct->numseq[ii+1]]
								[ct->numseq[k-1]][ct->numseq[ii]]
								*w3[k]*
								penalty(k-2,ii+1,ct,data);



							if((mod[ii+1]||mod[k-2]))if(inc[ct->numseq[ii+2]][ct->numseq[k-3]]&&notgu(ii+1,k-2,ct)&&!(fce.f(ii+1,k-2)&SINGLE)) {
								rarray+=v.f(ii+2,k-3)*data->tstack[ct->numseq[k-2]][ct->numseq[ii+1]]
									[ct->numseq[k-1]][ct->numseq[ii]]
									*pfchecknp(lfce[k-1],lfce[ii])*w3[k]*
									penalty(k-2,ii+1,ct,data)*erg1(ii+1,k-2,ii+2,k-3,ct,data);


							}
						}

						//also consider coaxial stacking:
						#ifndef disablecoax //a flag to disable coaxial stacking
						for (ip=k+minloop+1;ip<=number+1;ip++) {

				
							//first consider flush stacking:
							rarray+=v.f(ii,k-1)*v.f(k,ip-1)*w3[ip]*
								penalty(ii,k-1,ct,data)*penalty(k,ip-1,ct,data)*
								ergcoaxflushbases(ii,k-1,k,ip-1,ct,data);

							if(mod[ii]||mod[k-1]||mod[k]||mod[ip-1]) {

								if ((mod[ii]||mod[k-1])&&(mod[k]||mod[ip-1])&&inc[ct->numseq[ii+1]][ct->numseq[k-2]]
									&&inc[ct->numseq[k+1]][ct->numseq[ip-2]]&&notgu(ii,k-1,ct)&&notgu(k,ip-1,ct)
									&&!(fce.f(ii,k-1)&SINGLE)&&!(fce.f(k,ip-1)&SINGLE)) {

									rarray+=v.f(ii+1,k-2)*v.f(k+1,ip-2)*w3[ip]*
										penalty(ii,k-1,ct,data)*penalty(k,ip-1,ct,data)*
										ergcoaxflushbases(ii,k-1,k,ip-1,ct,data)
										*erg1(ii,k-1,ii+1,k-2,ct,data)*erg1(k,ip-1,k+1,ip-2,ct,data);
								}
								if((mod[ii]||mod[k-1])&&inc[ct->numseq[ii+1]][ct->numseq[k-2]]&&notgu(ii,k-1,ct)&&!(fce.f(ii,k-1)&SINGLE)) {
						
									rarray+=v.f(ii+1,k-2)*v.f(k,ip-1)*w3[ip]*
										penalty(ii,k-1,ct,data)*penalty(k,ip-1,ct,data)*
										ergcoaxflushbases(ii,k-1,k,ip-1,ct,data)
										*erg1(ii,k-1,ii+1,k-2,ct,data);

								}

								if((mod[k]||mod[ip-1])&&inc[ct->numseq[k+1]][ct->numseq[ip-2]]&&notgu(k,ip-1,ct)&&!(fce.f(k,ip-1)&SINGLE)) {

									rarray+=v.f(ii,k-1)*v.f(k+1,ip-2)*w3[ip]*
										penalty(ii,k-1,ct,data)*penalty(k,ip-1,ct,data)*
										ergcoaxflushbases(ii,k-1,k,ip-1,ct,data)
										*erg1(k,ip-1,k+1,ip-2,ct,data);
								}

							}


							//now consider an intervening mismatch:
							if (ii>0&&!lfce[ii]&&!lfce[k-1]) {
								rarray+=v.f(ii+1,k-2)*v.f(k,ip-1)*w3[ip]*
									penalty(ii+1,k-2,ct,data)*penalty(k,ip-1,ct,data)*
									ergcoaxinterbases1(ii+1,k-2,k,ip-1,ct,data);

								if(mod[ii+1]||mod[k-2]||mod[k]||mod[ip-1]){
						
									if((mod[ii+1]||mod[k-2])&&(mod[k]||mod[ip-1])&&inc[ct->numseq[ii+2]][ct->numseq[k-3]]
										&&inc[ct->numseq[k+1]][ct->numseq[ip-2]]&&notgu(ii+1,k-2,ct)&&notgu(k,ip-1,ct)
										&&!(fce.f(k,ip-1)&SINGLE)&&!(fce.f(ii+1,k-2)&SINGLE)){
										rarray+=v.f(ii+2,k-3)*v.f(k+1,ip-2)*w3[ip]*
											penalty(ii+1,k-2,ct,data)*penalty(k,ip-1,ct,data)*
											ergcoaxinterbases1(ii+1,k-2,k,ip-1,ct,data)
											*erg1(ii+1,k-2,ii+2,k-3,ct,data)*erg1(k,ip-1,k+1,ip-2,ct,data);

									}

									if((mod[ii+1]||mod[k-2])&&inc[ct->numseq[ii+2]][ct->numseq[k-3]]&&notgu(ii+1,k-2,ct)&&!(fce.f(ii+1,k-2)&SINGLE)) {
										rarray+=v.f(ii+2,k-3)*v.f(k,ip-1)*w3[ip]*
										penalty(ii+1,k-2,ct,data)*penalty(k,ip-1,ct,data)*
										ergcoaxinterbases1(ii+1,k-2,k,ip-1,ct,data)
										*erg1(ii+1,k-2,ii+2,k-3,ct,data);

									}
									if((mod[k]||mod[ip-1])&&inc[ct->numseq[k+1]][ct->numseq[ip-2]]&&notgu(k,ip-1,ct)&&!(fce.f(k,ip-1)&SINGLE)) {
										rarray+=v.f(ii+1,k-2)*v.f(k+1,ip-2)*w3[ip]*
											penalty(ii+1,k-2,ct,data)*penalty(k,ip-1,ct,data)*
											ergcoaxinterbases1(ii+1,k-2,k,ip-1,ct,data)
											*erg1(k,ip-1,k+1,ip-2,ct,data);	


									}


								}

							}
							if (!lfce[k-1]&&!lfce[ip-1]) {
					
								rarray+=v.f(ii,k-2)*v.f(k,ip-2)*w3[ip]*
									penalty(ii,k-2,ct,data)*penalty(k,ip-2,ct,data)*
									ergcoaxinterbases2(ii,k-2,k,ip-2,ct,data);

								if (mod[ii]||mod[k-2]||mod[k]||mod[ip-2]) {

									if ((mod[ii]||mod[k-2])&&(mod[k]||mod[ip-2])&&inc[ct->numseq[ii+1]][ct->numseq[k-3]]
										&&inc[ct->numseq[k+1]][ct->numseq[ip-3]]&&notgu(ii,k-2,ct)&&notgu(k,ip-2,ct)
										&&!(fce.f(ii,k-2)&SINGLE)&&!(fce.f(k,ip-2)&SINGLE)) {

										rarray+=v.f(ii+1,k-3)*v.f(k+1,ip-3)*w3[ip]*
											penalty(ii,k-2,ct,data)*penalty(k,ip-2,ct,data)*
											ergcoaxinterbases2(ii,k-2,k,ip-2,ct,data)
											*erg1(ii,k-2,ii+1,k-3,ct,data)*erg1(k,ip-2,k+1,ip-3,ct,data);
									}

									if ((mod[ii]||mod[k-2])&&inc[ct->numseq[ii+1]][ct->numseq[k-3]]&&notgu(ii,k-2,ct)&&!(fce.f(ii,k-2)&SINGLE)) {

										rarray+=v.f(ii+1,k-3)*v.f(k,ip-2)*w3[ip]*
											penalty(ii,k-2,ct,data)*penalty(k,ip-2,ct,data)*
											ergcoaxinterbases2(ii,k-2,k,ip-2,ct,data)
											*erg1(ii,k-2,ii+1,k-3,ct,data);
									}

									if ((mod[k]||mod[ip-2])&&inc[ct->numseq[k+1]][ct->numseq[ip-3]]&&notgu(k,ip-2,ct)&&!(fce.f(k,ip-2)&SINGLE)) {

										rarray+=v.f(ii,k-2)*v.f(k+1,ip-3)*w3[ip]*
											penalty(ii,k-2,ct,data)*penalty(k,ip-2,ct,data)*
											ergcoaxinterbases2(ii,k-2,k,ip-2,ct,data)
											*erg1(k,ip-2,k+1,ip-3,ct,data);
									}

								}
							}
					
						}
						#endif //ifndef disablecoax

					
					}
      	
					w3[ii] = rarray;
					//check to see if w5 is about to go out of bounds:
					if (w3[ii]>PFMAX) {
						rescale(1,number,ct,data,&v,&w,&wl,&wcoax,&wmb,&wmbl,w5,w3,wca,SCALEDOWN);
						twoscaling=twoscaling*SCALEDOWN*SCALEDOWN;
					}
					else if (w3[ii]<PFMIN&&w3[ii]>0) {
						rescale(1,number,ct,data,&v,&w,&wl,&wcoax,&wmb,&wmbl,w5,w3,wca,SCALEUP);
						twoscaling = twoscaling*SCALEUP*SCALEUP;
					}
				}
			}
  
   


			
		}
		//check to see if any of the 2-D arrays are about to go out oif bounds
		if (v.f(i,j)>PFMAX) {
			rescale(i,j,ct,data,&v,&w,&wl,&wcoax,&wmb,&wmbl,w5,w3,wca,SCALEDOWN);
			twoscaling = twoscaling*SCALEDOWN*SCALEDOWN;
		}
		else if (w.f(i,j)>PFMAX) {
			rescale(i,j,ct,data,&v,&w,&wl,&wcoax,&wmb,&wmbl,w5,w3,wca,SCALEDOWN);
			twoscaling = twoscaling*SCALEDOWN*SCALEDOWN;
		}
		else if (wl.f(i,j)>PFMAX) {
			rescale(i,j,ct,data,&v,&w,&wl,&wcoax,&wmb,&wmbl,w5,w3,wca,SCALEDOWN);
			twoscaling = twoscaling*SCALEDOWN*SCALEDOWN;
		}
		else if (wcoax.f(i,j)>PFMAX) {
			rescale(i,j,ct,data,&v,&w,&wl,&wcoax,&wmb,&wmbl,w5,w3,wca,SCALEDOWN);
			twoscaling = twoscaling*SCALEDOWN*SCALEDOWN;
		}
		else if (wmb.f(i,j)>PFMAX) {
			rescale(i,j,ct,data,&v,&w,&wl,&wcoax,&wmb,&wmbl,w5,w3,wca,SCALEDOWN);
			twoscaling = twoscaling*SCALEDOWN*SCALEDOWN;
		}
		else if (wmbl.f(i,j)>PFMAX) {
			rescale(i,j,ct,data,&v,&w,&wl,&wcoax,&wmb,&wmbl,w5,w3,wca,SCALEDOWN);
			twoscaling = twoscaling*SCALEDOWN*SCALEDOWN;
		}
		else if (v.f(i,j)<PFMIN&&v.f(i,j)>0) {
			rescale(i,j,ct,data,&v,&w,&wl,&wcoax,&wmb,&wmbl,w5,w3,wca,SCALEUP);
			twoscaling = twoscaling*SCALEUP*SCALEUP;
		}
		else if (w.f(i,j)<PFMIN&&w.f(i,j)>0) {
			rescale(i,j,ct,data,&v,&w,&wl,&wcoax,&wmb,&wmbl,w5,w3,wca,SCALEUP);
			twoscaling = twoscaling*SCALEUP*SCALEUP;
		}
		else if (wl.f(i,j)<PFMIN&&wl.f(i,j)>0) {
			rescale(i,j,ct,data,&v,&w,&wl,&wcoax,&wmb,&wmbl,w5,w3,wca,SCALEUP);
			twoscaling = twoscaling*SCALEUP*SCALEUP;
		}
		else if (wcoax.f(i,j)<PFMIN&&wcoax.f(i,j)>0) {
			rescale(i,j,ct,data,&v,&w,&wl,&wcoax,&wmb,&wmbl,w5,w3,wca,SCALEUP);
			twoscaling = twoscaling*SCALEUP*SCALEUP;
		}
		else if (wmb.f(i,j)<PFMIN&&wmb.f(i,j)>0) {
			rescale(i,j,ct,data,&v,&w,&wl,&wcoax,&wmb,&wmbl,w5,w3,wca,SCALEUP);
			twoscaling = twoscaling*SCALEUP*SCALEUP;
		}
		else if (wmbl.f(i,j)<PFMIN&&wmbl.f(i,j)>0) {
			rescale(i,j,ct,data,&v,&w,&wl,&wcoax,&wmb,&wmbl,w5,w3,wca,SCALEUP);
			twoscaling = twoscaling*SCALEUP*SCALEUP;
		}
   }





   
}






delete[] wca;

#ifdef timer
timeout << time(NULL)<<"\n";
timeout << time(NULL) - seconds;
timeout.close();
#endif

//////////////////////////
//output V, W, WMB, and W2V:
#if defined (pfdebugmode)
	ofstream foo;
	foo.open("arrays.out");
	foo << "i" << "\t"<<"j"<<"\t"<<"v.f(i,j)"<<"\t"<<"w.f(i,j)"<<"\t"<<"wmb.f(i,j)\twmbl.f(i,j)\twcoax.f(i,j)"<<"\t"<<"wl.f(i,j)"<<"\t"<<"v.f(j,i+number)"<<"\t"<<"w.f(j,i+number)"<<"\t"<<"wmb.f(j,i+number)"<<"\t"<<"wl.f(j,i+number)"<<"\t"<<"wmbl.f(j,i+numer)\twcoax.f(j,i+number)"<<"\n";
	for (j=1;j<=number;j++) {
		for (i=1;i<=j;i++) {

			foo << i << "\t"<<j<<"\t"<<v.f(i,j)<<"\t"<<w.f(i,j)<<"\t"<<wmb.f(i,j)<<"\t"<<wmbl.f(i,j)<<"\t"<<wcoax.f(i,j)<<"\t"<<wl.f(i,j)<<"\t"<<v.f(j,i+number)<<"\t"<<w.f(j,i+number)<<"\t"<<wmb.f(j,i+number)<<"\t"<<wl.f(j,i+number)<<"\t"<<wmbl.f(j,i+number)<<"\t"<<wcoax.f(j,i+number)<<"\n";

		}
	}

	foo <<"\n\n\n";
	foo << "i" << "\t" << "w5[i]" << "\t" << "w3[i]" << "\n";
	for (i=0;i<=number;i++) foo << i << "\t" << w5[i] << "\t" << w3[i] << "\n";

	foo.close();

#endif


if (save!=0) {
	ofstream sav(save,ios::binary);
	
	//write the save file information so that the sequence can be re-folded,
		//include thermodynamic data to prevent traceback errors

	//start with structure information
	write(&sav,&(ct->numofbases));
	write(&sav,&(ct->intermolecular));
	write(&sav,&data->scaling);
	
	write(&sav,&(ct->npair));
	for (i=0;i<=ct->npair;i++) {
		write(&sav,&(ct->pair[i][0]));
		write(&sav,&(ct->pair[i][1]));
	}
	for (i=0;i<=ct->numofbases;i++) {
		
		write(&sav,&(ct->hnumber[i]));
		sav.write(&(ct->nucs[i]),1);

	}

	for (i=0;i<=2*ct->numofbases;i++) write(&sav,&(ct->numseq[i]));
	
	write(&sav,&(ct->ndbl));
	for (i=0;i<=ct->ndbl;i++) write(&sav,&(ct->dbl[i]));

	
	if (ct->intermolecular) {
		for (i=0;i<3;i++) write(&sav,&(ct->inter[i]));
		
	}

	write(&sav,&(ct->nnopair));
	for (i=0;i<=ct->nnopair;i++) write(&sav,&(ct->nopair[i]));

	write(&sav,&(ct->nmod));
	for (i=0;i<=ct->nmod;i++) write(&sav,&(ct->mod[i]));
	
	write(&sav,&(ct->ngu));
	for (i=0;i<=ct->ngu;i++) write(&sav,&(ct->gu[i]));  
	
	write(&sav,ct->ctlabel[1]);

	write(&sav,&(ct->templated));
	if (ct->templated) {
		for (i=0;i<=ct->numofbases;i++) {
			for (j=0;j<=i;j++) write(&sav,&(ct->tem[i][j]));	

		}

	}


	//now write the array class data for v, w, and wmb:
	for (i=0;i<=ct->numofbases;i++) {
		write(&sav,&(w3[i]));
		write(&sav,&(w5[i]));
		for (j=0;j<=ct->numofbases;j++) {
			write(&sav,&(v.dg[i][j]));
			write(&sav,&(w.dg[i][j]));
			write(&sav,&(wmb.dg[i][j]));
			write(&sav,&(wmbl.dg[i][j]));
			write(&sav,&(wl.dg[i][j]));
			write(&sav,&(wcoax.dg[i][j]));
			writesinglechar(&sav,&(fce.dg[i][j]));
			//if (ct->intermolecular) {
			//	write(&sav,&(w2->dg[i][j]));
			//	write(&sav,&(wmb2->dg[i][j]));

			//}


		}	
		

	}

	write(&sav,&(w3[ct->numofbases+1]));
	for (i=0;i<=2*ct->numofbases;i++) {
		write(&sav,&(lfce[i]));
		write(&sav,&(mod[i]));

	}
	

	


	//now write the thermodynamic data:
	write(&sav,&(data->temp));
	for (i=0;i<5;i++) write(&sav,&(data->poppen[i]));
	write(&sav,&(data->maxpen));
	for (i=0;i<11;i++) write(&sav,&(data->eparam[i]));
	for (i=0;i<31;i++) {
		write(&sav,&(data->inter[i]));
		write(&sav,&(data->bulge[i]));
		write(&sav,&(data->hairpin[i]));

	}
	for (i=0;i<6;i++) {
		for (j=0;j<6;j++) {
			for (k=0;k<6;k++) {
				for (l=0;l<3;l++) {
					write(&sav,&(data->dangle[i][j][k][l]));	
				}
				for (l=0;l<6;l++) {
					write(&sav,&(data->stack[i][j][k][l]));
					write(&sav,&(data->tstkh[i][j][k][l]));
					write(&sav,&(data->tstki[i][j][k][l]));
					write(&sav,&(data->coax[i][j][k][l]));
					write(&sav,&(data->tstackcoax[i][j][k][l]));
					write(&sav,&(data->coaxstack[i][j][k][l]));
					write(&sav,&(data->tstack[i][j][k][l]));
					write(&sav,&(data->tstkm[i][j][k][l]));
					write(&sav,&(data->tstki23[i][j][k][l]));
					write(&sav,&(data->tstki1n[i][j][k][l]));
					for (m=0;m<6;m++) {
						for (n=0;n<6;n++) {
							write(&sav,&(data->iloop11[i][j][k][l][m][n]));
							for (o=0;o<6;o++) {
								if (inc[i][j]&&inc[n][o]) write(&sav,&(data->iloop21[i][j][k][l][m][n][o]));
								for (p=0;p<6;p++) {
									if (inc[i][k]&&inc[j][l])
										write(&sav,&(data->iloop22[i][j][k][l][m][n][o][p]));
								}
							}
							

						}
					}
				}
			}
		}
	}
	write(&sav,&(data->numoftloops));
	for (i=0;i<=data->numoftloops;i++) {
		write(&sav,&(data->itloop[i]));
		write(&sav,&(data->tloop[i]));

	}
	write(&sav,&(data->numoftriloops));
	for (i=0;i<=data->numoftriloops;i++) {
		write(&sav,&(data->itriloop[i]));
		write(&sav,&(data->triloop[i]));

	}
	write(&sav,&(data->numofhexaloops));
	for (i=0;i<=data->numofhexaloops;i++) {
		write(&sav,&(data->ihexaloop[i]));
		write(&sav,&(data->hexaloop[i]));

	}
	write(&sav,&(data->auend));
	write(&sav,&(data->gubonus));
	write(&sav,&(data->cint));
	write(&sav,&(data->cslope));
	write(&sav,&(data->c3));
	write(&sav,&(data->efn2a));
	write(&sav,&(data->efn2b));
	write(&sav,&(data->efn2c));
	write(&sav,&(data->init));
	write(&sav,&(data->mlasym));
	write(&sav,&(data->strain));
	write(&sav,&(data->prelog));
	write(&sav,&(data->singlecbulge));
	write(&sav,&(data->maxintloopsize));
	

	
	sav.close();
}

if (quickQ) *Q = w5[ct->numofbases];


delete[] lfce;
delete[] mod;




delete[] w5;
delete[] w3;


/*if (ct->intermolecular) {
	delete w2;
	delete wmb2;
}*/



return;
}


////////////////////////////////////////////////////////////////////////
//pfunctionclass encapsulates the large 2-d arrays of w and v, used by the
//	partition function

      //the constructor allocates the space needed by the arrays
pfunctionclass::pfunctionclass(int size) {
	//zero indicates whether the array should be set to zero as opposed
		//to being set to infinity, it is false by default
      	

      	
	infinite = 0;

    Size = size;
    register int i,j;
    dg = new PFPRECISION *[size+1];

	for (i=0;i<=(size);i++)  {
   		dg[i] = new PFPRECISION [size+1];
   	}
    for (i=0;i<=size;i++) {
         for (j=0;j<size+1;j++) {
			 
             dg[i][j] = 0;
         }
    }

}

//the destructor deallocates the space used
pfunctionclass::~pfunctionclass() {

      	
	int i;
       	
    for (i=0;i<=Size;i++) {
        delete[] dg[i];
    }
     delete[] dg;
}

      //f is an integer function that references the correct element of the array
inline PFPRECISION &pfunctionclass::f(int i, int j) {
   


   if (i>j) {
        return infinite;
    }
   else if (i>Size) return f(i-Size,j-Size);
   else return dg[i][j-i];
         
}

//When considering mismatch at the end of a helix, consult this function to check
//	whether the nucs are required to pair
PFPRECISION pfchecknp(bool lfce1,bool lfce2) {

	if (lfce1||lfce2) return 0;
	else return 1;
}

inline PFPRECISION boltzman(short i, PFPRECISION temp) {

	if (i==infinity) return 0;
	else return exp((-((PFPRECISION) i)/((PFPRECISION)conversionfactor))/(RKC*temp));

}


pfdatatable::pfdatatable() {
}

pfdatatable::pfdatatable(datatable *data, PFPRECISION Scaling, PFPRECISION  Temp) {
	//the partition function datatable needs to be initialized from the datatable
	short i,j,k,l,m,n,o,p;

	scaling = Scaling;

	//derive the temperature from the datatable
	temp = Temp;

	for (i=0;i<5;i++) poppen[i]=boltzman(data->poppen[i],temp);
	maxpen = boltzman(data->maxpen,temp);
	for (i=0;i<11;i++) eparam[i] = boltzman(data->eparam[i],temp);
	maxintloopsize = data->eparam[7];
	for (i=0;i<31;i++) {
		inter[i] = boltzman(data->inter[i],temp)*pow(scaling,i+2);
		bulge[i] = boltzman(data->bulge[i],temp)*pow(scaling,i+2);
		hairpin[i] = boltzman(data->hairpin[i],temp)*pow(scaling,i+2);

	}
	for (i=0;i<6;i++) {
		for (j=0;j<6;j++) {
			for (k=0;k<6;k++) {
				for (l=0;l<3;l++) {
					dangle[i][j][k][l] = boltzman(data->dangle[i][j][k][l],temp)*scaling;	
				}
				for (l=0;l<6;l++) {
					stack[i][j][k][l]=boltzman(data->stack[i][j][k][l],temp)*pow(scaling,2);
					tstkh[i][j][k][l]=boltzman(data->tstkh[i][j][k][l],temp);
					tstki[i][j][k][l]=boltzman(data->tstki[i][j][k][l],temp);
					coax[i][j][k][l]=boltzman(data->coax[i][j][k][l],temp);
					tstackcoax[i][j][k][l]=boltzman(data->tstackcoax[i][j][k][l],temp)*pow(scaling,2);
					coaxstack[i][j][k][l] = boltzman(data->coaxstack[i][j][k][l],temp);
					tstack[i][j][k][l]=boltzman(data->tstack[i][j][k][l],temp)*pow(scaling,2);
					tstkm[i][j][k][l]=boltzman(data->tstkm[i][j][k][l],temp)*pow(scaling,2);
					tstki23[i][j][k][l]=boltzman(data->tstki23[i][j][k][l],temp);
					tstki1n[i][j][k][l]=boltzman(data->tstki1n[i][j][k][l],temp);
					for (m=0;m<6;m++) {
						for (n=0;n<6;n++) {
							iloop11[i][j][k][l][m][n]=boltzman(data->iloop11[i][j][k][l][m][n],temp)*pow(scaling,4);
							for (o=0;o<6;o++) {
								iloop21[i][j][k][l][m][n][o]=boltzman(data->iloop21[i][j][k][l][m][n][o],temp)*pow(scaling,5);
								for (p=0;p<6;p++) {
									iloop22[i][j][k][l][m][n][o][p]=boltzman(data->iloop22[i][j][k][l][m][n][o][p],temp)*pow(scaling,6);
								}
							}
							

						}
					}
				}
			}
		}
	}
	numoftloops = data->numoftloops;
	for (i=0;i<=data->numoftloops;i++) {
		itloop[i]=data->tloop[i][0];
		tloop[i] = boltzman(data->tloop[i][1],temp)*pow(scaling,6);

	}
	numoftriloops=data->numoftriloops;
	for (i=0;i<=data->numoftriloops;i++) {
		itriloop[i] = data->triloop[i][0];
		triloop[i] = boltzman(data->triloop[i][1],temp)*pow(scaling,5);
	}
	numofhexaloops=data->numofhexaloops;
	for (i=0;i<=data->numofhexaloops;i++) {
		ihexaloop[i]=data->hexaloop[i][0];
		hexaloop[i] = boltzman(data->hexaloop[i][1],temp)*pow(scaling,8);

	}
	auend = boltzman(data->auend,temp);
	gubonus = boltzman(data->gubonus,temp);
	cint = boltzman(data->cint,temp);
	cslope = boltzman(data->cslope,temp);
	c3 = boltzman(data->c3,temp);
	efn2a = boltzman(data->efn2a,temp);
	efn2b = boltzman(data->efn2b,temp);
	efn2c = boltzman(data->efn2c,temp);
	init = boltzman(data->init,temp);
	mlasym = boltzman(data->mlasym,temp);
	strain = boltzman(data->strain,temp);
	prelog = data->prelog/conversionfactor;
	singlecbulge = boltzman(data->singlecbulge,temp);





}



//calculate the energy of stacked base pairs

//oriented by:
//5' i ip 3'
//   |  |
//   j jp

PFPRECISION erg1(int i,int j,int ip,int jp,structure *ct, pfdatatable *data)
{

		PFPRECISION energy;

		 if ((i==(ct->numofbases))||(j==((ct->numofbases)+1))) {
      		//this is not allowed because n and n+1 are not cavalently attached
			energy = 0;
		}
		else {
      		energy = data->stack[(ct->numseq[i])][(ct->numseq[j])]
				[(ct->numseq[ip])][(ct->numseq[jp])]*data->eparam[1];
		}

		#ifdef equiout
		ofstream kout;
		kout.open("k.out",ofstream::app);
		kout << "k base pair stack ("<<i<<","<<j<<","<<ip<<","<<jp<<") ="<<energy<< "\n";
		kout.close();
		#endif

		return energy;
}


//calculate the energy of a bulge/internal loop
//where i is paired to j; ip is paired to jp; ip > i; j > jp
PFPRECISION erg2(int i,int j,int ip,int jp,structure *ct, pfdatatable *data,
	char a, char b)
{

	int size,size1,size2, lopsid, count;
	PFPRECISION energy,loginc;
	/* size,size1,size2 = size of a loop
		energy = energy calculated
		loginc = the value of a log used in large hairpin loops
	*/
   	if (((i<=(ct->numofbases))&&(ip>(ct->numofbases)))||((
      	jp<=(ct->numofbases))&&(j>(ct->numofbases)))) {
         //A loop cannot contain the ends of the sequence
         
         return 0;
      }


	
      size1 = ip-i-1;
		size2 = j - jp - 1;

	if ((a>0)||(b>0)) {
      	if ((a&DOUBLE)||(b&DOUBLE)) return 0;//the loop contains a nuc that
      		//should be double stranded
/*      	else if ((a&INTER)) {
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
	  }*/
	  /*else if (b&INTER) {
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

         }*/
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
						*data->bulge[size]*data->eparam[2]/(data->scaling*data->scaling);
				if (size1==1)  {
					
					//count the number of alternative bulges that exist:
					
					//k = i;
					//while (ct->numseq[k]==ct->numseq[i+1]) {
					//	count++;
					//	k--;
					//}
					//k=ip;
					//while(ct->numseq[k]==ct->numseq[i+1]) {
					//	count++;
					//	k++;
					//}
					//give bonus to C adjacent to single C bulge
					if (ct->numseq[i+1]==2&&(ct->numseq[i+2]==2||ct->numseq[i]==2)) energy= energy*data->singlecbulge;
					
				}
				
				else {
					//size2 == 1
					
					//count the number of alternative bulges that exist:
					
					//k = jp;
					//while (ct->numseq[k]==ct->numseq[jp+1]) {
					//	count++;
					//	k--;
					//}
					//k=j;
					//while(ct->numseq[k]==ct->numseq[jp+1]) {
					//	count++;
					//	k++;
					//}
					//give bonus to C adjacent to single C bulge
					if (ct->numseq[j-1]==2&&(ct->numseq[j-2]==2||ct->numseq[j]==2)) energy=energy*data->singlecbulge;
					
				}
				//do not apply a correction for the number of equivalent states because 
					//the bulge can move to adjacent sites in the partition function calc
				//energy-= (int) (rt*conversionfactor* log ((double) count));
			}
			else if (size>30) {

				loginc = ((data->prelog)*log(PFPRECISION ((size)/30.0)));
				energy = data->bulge[30]*exp(-loginc/(RKC*data->temp))*data->eparam[2];
				energy = energy*penalty(i,j,ct,data)*penalty(jp,ip,ct,data)*pow(data->scaling,size-30);

			}
			else {
         		energy = data->bulge[size]*data->eparam[2];
				energy = energy*penalty(i,j,ct,data)*penalty(jp,ip,ct,data);
			}
		}
		else {//internal loop
			size = size1 + size2;
			lopsid = abs(size1-size2);

			if (size>30) {

				loginc = ((data->prelog)*log((PFPRECISION ((size))/30.0)));
				if (size1==1||size2==1) {
            		energy = data->tstki1n[ct->numseq[i]][ct->numseq[j]]
						[ct->numseq[i+1]][ct->numseq[j-1]]*
						data->tstki1n[ct->numseq[jp]][ct->numseq[ip]]
						[ct->numseq[jp+1]][ct->numseq[ip-1]]*
						data->inter[30]* exp(-loginc/(RKC*data->temp)) *data->eparam[3]*
						max(data->maxpen,
						pow(data->poppen[min(2,min(size1,size2))],lopsid))
						*pow(data->scaling,size-30);

				}

				else {
					energy = data->tstki[ct->numseq[i]][ct->numseq[j]]
						[ct->numseq[i+1]][ct->numseq[j-1]]*
						data->tstki[ct->numseq[jp]][ct->numseq[ip]]
						[ct->numseq[jp+1]][ct->numseq[ip-1]]*
						data->inter[30]* exp(-loginc/(RKC*data->temp))*data->eparam[3] *
						max(data->maxpen,
						pow(data->poppen[min(2,min(size1,size2))],lopsid))
						*pow(data->scaling,size-30);
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
				energy = data->tstki1n[ct->numseq[i]][ct->numseq[j]][ct->numseq[i+1]][ct->numseq[j-1]] *
						data->tstki1n[ct->numseq[jp]][ct->numseq[ip]][ct->numseq[jp+1]][ct->numseq[ip-1]] *
						data->inter[size] * data->eparam[3] *
					max(data->maxpen,pow(data->poppen[min(2,min(size1,size2))],lopsid));
			}
        

			else if ((size1==2&&size2==3)||(size1==3&&size2==2)) {
			//this is a 2x3 loop
				energy = data->tstki23[ct->numseq[i]][ct->numseq[j]][ct->numseq[i+1]][ct->numseq[j-1]] *
					data->tstki23[ct->numseq[jp]][ct->numseq[ip]][ct->numseq[jp+1]][ct->numseq[ip-1]] *
					data->inter[size] * data->eparam[3] *
					max(data->maxpen,pow(data->poppen[min(2,min(size1,size2))],lopsid));


			}
			else {
         	


         		energy = data->tstki[ct->numseq[i]][ct->numseq[j]][ct->numseq[i+1]][ct->numseq[j-1]] *
					data->tstki[ct->numseq[jp]][ct->numseq[ip]][ct->numseq[jp+1]][ct->numseq[ip-1]] *
					data->inter[size] * data->eparam[3] *
					max(data->maxpen,pow(data->poppen[min(2,min(size1,size2))],lopsid));
			}
		}
		#ifdef equiout
		ofstream kout;
		kout.open("k.out",ofstream::app);
		kout << "k internal/bulge ("<<i<<","<<j<<","<<ip<<","<<jp<<") ="<<energy<< "\n";
		kout.close();
		#endif

		return energy;
}



//calculate the energy of a hairpin loop:
#ifndef equiout
PFPRECISION erg3(int i,int j,structure *ct, pfdatatable *data,char dbl) {
#else
PFPRECISION erg3indirect(int i,int j,structure *ct, pfdatatable *data,char dbl) {
#endif
int size,count,key,k;
PFPRECISION energy,loginc;
	/* size,size1,size2 = size of a loop
		energy = energy calculated
		loginc = the value of a log used in large hairpin loops
	*/


		if (dbl&DOUBLE) return 0;//the loop contains a base that should be
      										//double stranded

      else if (dbl&INTER) {//intermolecular interaction
      	//intermolecular "hairpin" free energy is that of intermolecular
         //	initiation plus the stacked mismatch

         energy = data->init * data->tstack[ct->numseq[i]][ct->numseq[j]]
         	[ct->numseq[i+1]][ct->numseq[j-1]]*penalty(i,j,ct,data);

         return energy;
      }


   	if ((i<=(ct->numofbases))&&(j>(ct->numofbases))) {
      	//A hairpin cannot contain the ends of the sequence
         
         return 0;
      }

		size = j-i-1;



		if (size>30) {

			loginc = ((data->prelog)*log((PFPRECISION ((size))/30.0)));



			energy = data->tstkh[ct->numseq[i]][ct->numseq[j]]
				[ct->numseq[i+1]][ct->numseq[j-1]]
				* data->hairpin[30]*exp(-loginc/(RKC*data->temp))*data->eparam[4]*pow(data->scaling,size-30);
		}
		else if (size<3) {
      		energy = data->hairpin[size]*data->eparam[4];
				if (ct->numseq[i]==4||ct->numseq[j]==4) energy = energy*exp(-.6/(RKC*data->temp));
		}
		else if (size==4) {
			
			key = (ct->numseq[j])*3125 + (ct->numseq[i+4])*625 +
				(ct->numseq[i+3])*125 + (ct->numseq[i+2])*25+(ct->numseq[i+1])*5+(ct->numseq[i]);
			for (count=1;count<=data->numoftloops;count++) {
					if (key==data->itloop[count]) return data->tloop[count];
			}
			energy = data->tstkh[ct->numseq[i]][ct->numseq[j]]
				[ct->numseq[i+1]][ct->numseq[j-1]]
				* data->hairpin[size] * data->eparam[4];
		}
		else if (size==3) {
			
			key = (ct->numseq[j])*625 +
				(ct->numseq[i+3])*125 + (ct->numseq[i+2])*25+(ct->numseq[i+1])*5+(ct->numseq[i]);
			for (count=1;count<=data->numoftriloops;count++) {
				if (key==data->itriloop[count]) return data->triloop[count];
			}
			
			energy =	data->hairpin[size] * data->eparam[4]
         	*penalty(i,j,ct,data);
		}
		else if (size==6) {
			key = (ct->numseq[j])*78125 + (ct->numseq[i+6])*15625 + (ct->numseq[i+5])*3125 
				+ (ct->numseq[i+4])*625 +
				(ct->numseq[i+3])*125 + (ct->numseq[i+2])*25+(ct->numseq[i+1])*5+(ct->numseq[i]);
			for (count=1;count<=data->numofhexaloops;count++) {
				if (key==data->ihexaloop[count]) return data->hexaloop[count];
			}

			energy = data->tstkh[ct->numseq[i]][ct->numseq[j]]
				[ct->numseq[i+1]][ct->numseq[j-1]]
				* data->hairpin[size] * data->eparam[4];
		}

		else {
			energy = data->tstkh[ct->numseq[i]][ct->numseq[j]]
				[ct->numseq[i+1]][ct->numseq[j-1]]
				* data->hairpin[size] * data->eparam[4];
		}




		//check for GU closeure preceded by GG
      if (ct->numseq[i]==3&&ct->numseq[j]==4) {
      	if ((i>2&&i<ct->numofbases)||(i>ct->numofbases+2))
       		if (ct->numseq[i-1]==3&&ct->numseq[i-2]==3) {

         		energy = energy * data->gubonus;
            	


         	}
      }

      //check for an oligo-c loop
	 
      for (k=1;(k<=size);k++) {
       	if (ct->numseq[i+k] != 2) return energy;//this is not an oligo-C loop
      }
      //this is a poly c loop so penalize
      if (size==3) return (energy *data->c3);
      else return (energy * data->cint * pow(data->cslope,size));
      

}

#ifdef equiout
PFPRECISION erg3(int i,int j,structure *ct, pfdatatable *data,char dbl) {
	PFPRECISION energy;
	
	energy = erg3indirect(i,j,ct, data,dbl,temp);
	ofstream kout;
	kout.open("k.out",ofstream::app);
	kout << "k hairpin ("<<i<<","<<j<<") ="<<energy<< "\n";
	kout.close();
	return energy;
}
#endif

//calculate the energy of a dangling end:
PFPRECISION erg4(int i,int j,int ip,int jp,structure *ct, pfdatatable *data, bool lfce)
{

//dangling base
		// jp = 1 => 3' dangle
		// jp = 2 => 5' dangle



      if (lfce) return 0;//stacked nuc should be double stranded

	  //commented out 11/8/99
      //if (ip==5) return 0;//dangling nuc is an intermolecular linker
	    #ifdef equiout
			ofstream kout;
			kout.open("k.out",ofstream::app);
			kout << "k dangle ("<<i<<","<<j<<","<<ip<<","<<jp<<") ="<<data->dangle[ct->numseq[i]][ct->numseq[j]][ct->numseq[ip]][jp]<< "\n";
			kout.close();
		#endif

		return data->dangle[ct->numseq[i]][ct->numseq[j]][ct->numseq[ip]][jp];
		
}

//this function calculates whether a terminal pair i,j requires the end penalty
PFPRECISION penalty(int i,int j,structure* ct, pfdatatable *data) {
	
	#ifdef equiout
		ofstream kout;
		kout.open("k.out",ofstream::app);
		if (ct->numseq[i]==4||ct->numseq[j]==4) kout << "k penalty ("<<i<<","<<j<<") ="<<data->auend<< "\n";
		else kout << "k penalty ("<<i<<","<<j<<") ="<<1.0<< "\n";
		kout.close();
	#endif
   
	if (ct->numseq[i]==4||ct->numseq[j]==4)
   	return data->auend;
	else return 1;//no end penalty
	



}

//this function calculates whether a terminal pair i,j requires the end penalty
PFPRECISION penalty2(int i,int j, pfdatatable *data) {


   if (i==4||j==4)
   	return data->auend;
   else return 1;//no end penalty


}



//this function calculates the free energy for coaxial stacking of pair i-j onto ip-jp
//	require that the backbone continues directly between j and ip without a nucleotide in
//	a canonical pair
//k indicates the intervening nuc in intervening mismatch that is not sandwiched between
//	j and ip
//#ifdef equiout
//PFPRECISION ergcoaxindirect(int i, int j, int ip, int jp, int k, structure *ct, pfdatatable *data) {
//#else
PFPRECISION ergcoaxflushbases(int i, int j, int ip, int jp, structure *ct, pfdatatable *data) {
//#endif
	
		//flush stacking
		
		return data->coax[ct->numseq[i]][ct->numseq[j]][ct->numseq[ip]][ct->numseq[jp]];

	

	

}

PFPRECISION ergcoaxinterbases1(int i, int j, int ip, int jp, structure *ct, pfdatatable *data) {
//#endif
		//coaxial stacking with an intervening mismatch
		//(k==i-1) 
			
			return data->tstackcoax[ct->numseq[j]][ct->numseq[i]]
				[ct->numseq[j+1]][ct->numseq[i-1]] * 
				data->coaxstack[ct->numseq[j+1]][ct->numseq[i-1]]
				[ct->numseq[ip]][ct->numseq[jp]];	
			

}

PFPRECISION ergcoaxinterbases2(int i, int j, int ip, int jp, structure *ct, pfdatatable *data) {
//#endif
	//coaxial stacking with an intervening mismatch
	/*(k==jp+1) {*/
			
			return data->tstackcoax[ct->numseq[jp]][ct->numseq[ip]]
				[ct->numseq[jp+1]][ct->numseq[ip-1]] *
				data->coaxstack[ct->numseq[j]][ct->numseq[i]]
				[ct->numseq[j+1]][ct->numseq[jp+1]];	

}

//#ifdef equiout
//PFPRECISION ergcoax(int i, int j, int ip, int jp, int k, structure *ct, pfdatatable *data) {
//PFPRECISION energy;

//	energy = ergcoaxindirect(i,j,ip,jp,k,ct,data);

//	ofstream kout;

//	if (((i<j&&j<ip&&ip<jp)||(i>j&&j>ip&&ip>jp))&&k!=i&&k!=j&&k!=ip&&k!=jp) {
//		kout.open("k.out",ofstream::app);
//		kout << "k coax ("<<i<<","<<j<<","<<ip<<","<<jp<<","<<k<<") ="<<energy<< "\n";
//		kout.close();
//	}	
//	return energy;
//}
//#endif


void readpfsave(char *filename, structure *ct, 
			 PFPRECISION *w5, PFPRECISION *w3, 
			 pfunctionclass *v, pfunctionclass *w, pfunctionclass *wmb, pfunctionclass *wl, pfunctionclass *wmbl, pfunctionclass *wcoax,
			 forceclass *fce, PFPRECISION *scaling, bool *mod, bool *lfce, pfdatatable *data) {
	 int i,j,k,l,m,n,o,p;
	ifstream sav(filename,ios::binary);
	int inc[6][6]={{0,0,0,0,0,0},{0,0,0,0,1,0},{0,0,0,1,0,0},{0,0,1,0,1,0},
	{0,1,0,1,0,0},{0,0,0,0,0,0}};
	
	//read the save file

	//start with structure information
	read(&sav,&(ct->numofbases));
	read(&sav,&(ct->intermolecular));
	read(&sav,scaling);
	read(&sav,&(ct->npair));
	for (i=0;i<=ct->npair;i++) {
		read(&sav,&(ct->pair[i][0]));
		read(&sav,&(ct->pair[i][1]));
	}
	for (i=0;i<=ct->numofbases;i++) {
		
		read(&sav,&(ct->hnumber[i]));
		sav.read(&(ct->nucs[i]),1);

	}

	for (i=0;i<=2*ct->numofbases;i++) read(&sav,&(ct->numseq[i]));
	
	read(&sav,&(ct->ndbl));
	for (i=0;i<=ct->ndbl;i++) read(&sav,&(ct->dbl[i]));

	
	if (ct->intermolecular) {
		for (i=0;i<3;i++) read(&sav,&(ct->inter[i]));
		
	}

	read(&sav,&(ct->nnopair));
	for (i=0;i<=ct->nnopair;i++) read(&sav,&(ct->nopair[i]));

	read(&sav,&(ct->nmod));
	for (i=0;i<=ct->nmod;i++) read(&sav,&(ct->mod[i]));
	
	read(&sav,&(ct->ngu));
	for (i=0;i<=ct->ngu;i++) read(&sav,&(ct->gu[i]));  
	
	read(&sav,ct->ctlabel[1]);

	read(&sav,&(ct->templated));
	if (ct->templated) {
		ct->allocatetem();
		for (i=0;i<=ct->numofbases;i++) {
			for (j=0;j<=i;j++) read(&sav,&(ct->tem[i][j]));	

		}

	}


	//now read the array class data for v, w, and wmb:
	for (i=0;i<=ct->numofbases;i++) {
		read(&sav,&(w3[i]));
		read(&sav,&(w5[i]));
		for (j=0;j<=ct->numofbases;j++) {
			read(&sav,&(v->dg[i][j]));
			read(&sav,&(w->dg[i][j]));
			read(&sav,&(wmb->dg[i][j]));
			read(&sav,&(wmbl->dg[i][j]));
			read(&sav,&(wl->dg[i][j]));
			read(&sav,&(wcoax->dg[i][j]));
			readsinglechar(&sav,&(fce->dg[i][j]));
			//if (ct->intermolecular) {
			//	read(&sav,&(w2->dg[i][j]));
			//	read(&sav,&(wmb2->dg[i][j]));

			//}


		}	
		

	}

	read(&sav,&(w3[ct->numofbases+1]));
	for (i=0;i<=2*ct->numofbases;i++) {
		read(&sav,&(lfce[i]));
		read(&sav,&(mod[i]));

	}
	

	


	//now read the thermodynamic data:
	read(&sav,&(data->temp));
	for (i=0;i<5;i++) read(&sav,&(data->poppen[i]));
	read(&sav,&(data->maxpen));
	for (i=0;i<11;i++) read(&sav,&(data->eparam[i]));
	for (i=0;i<31;i++) {
		read(&sav,&(data->inter[i]));
		read(&sav,&(data->bulge[i]));
		read(&sav,&(data->hairpin[i]));

	}
	for (i=0;i<6;i++) {
		for (j=0;j<6;j++) {
			for (k=0;k<6;k++) {
				for (l=0;l<3;l++) {
					read(&sav,&(data->dangle[i][j][k][l]));	
				}
				for (l=0;l<6;l++) {
					read(&sav,&(data->stack[i][j][k][l]));
					read(&sav,&(data->tstkh[i][j][k][l]));
					read(&sav,&(data->tstki[i][j][k][l]));
					read(&sav,&(data->coax[i][j][k][l]));
					read(&sav,&(data->tstackcoax[i][j][k][l]));
					read(&sav,&(data->coaxstack[i][j][k][l]));
					read(&sav,&(data->tstack[i][j][k][l]));
					read(&sav,&(data->tstkm[i][j][k][l]));
					read(&sav,&(data->tstki23[i][j][k][l]));
					read(&sav,&(data->tstki1n[i][j][k][l]));
					for (m=0;m<6;m++) {
						for (n=0;n<6;n++) {
							read(&sav,&(data->iloop11[i][j][k][l][m][n]));
							for (o=0;o<6;o++) {
								if (inc[i][j]&&inc[n][o]) read(&sav,&(data->iloop21[i][j][k][l][m][n][o]));
								for (p=0;p<6;p++) {
									if (inc[i][k]&&inc[j][l])
										read(&sav,&(data->iloop22[i][j][k][l][m][n][o][p]));
								}
							}
							

						}
					}
				}
			}
		}
	}
	read(&sav,&(data->numoftloops));
	for (i=0;i<=data->numoftloops;i++) {
		read(&sav,&(data->itloop[i]));
		read(&sav,&(data->tloop[i]));

	}
	read(&sav,&(data->numoftriloops));
	for (i=0;i<=data->numoftriloops;i++) {
		read(&sav,&(data->itriloop[i]));
		read(&sav,&(data->triloop[i]));

	}
	read(&sav,&(data->numofhexaloops));
	for (i=0;i<=data->numofhexaloops;i++) {
		read(&sav,&(data->ihexaloop[i]));
		read(&sav,&(data->hexaloop[i]));

	}
	read(&sav,&(data->auend));
	read(&sav,&(data->gubonus));
	read(&sav,&(data->cint));
	read(&sav,&(data->cslope));
	read(&sav,&(data->c3));
	read(&sav,&(data->efn2a));
	read(&sav,&(data->efn2b));
	read(&sav,&(data->efn2c));
	read(&sav,&(data->init));
	read(&sav,&(data->mlasym));
	read(&sav,&(data->strain));
	read(&sav,&(data->prelog));
	read(&sav,&(data->singlecbulge));
	read(&sav,&(data->maxintloopsize));
	

	
	sav.close();
}


PFPRECISION calculateprobability(int i, int j, pfunctionclass *v, PFPRECISION *w5, structure *ct, pfdatatable *data, bool *lfce, bool *mod, PFPRECISION scaling, forceclass *fce) {
	PFPRECISION interior, exterior;
	short before,after;
	int inc[6][6]={{0,0,0,0,0,0},{0,0,0,0,1,0},{0,0,0,1,0,0},{0,0,1,0,1,0},
	{0,1,0,1,0,0},{0,0,0,0,0,0}};
	bool adjacentgu;

	if (!mod[i]&&!mod[j]) return (v->f(i,j)*v->f(j,i+ct->numofbases))/(w5[ct->numofbases]*scaling*scaling);
	else {
		if (!(fce->f(i,j)&SINGLE)) {
			before =0;
			if ((i>1&&j<(2*ct->numofbases)&&j!=ct->numofbases)) {
				if ((j>ct->numofbases&&((i-j+ct->numofbases)>minloop+2))||j<ct->numofbases) {
					before = inc[ct->numseq[i-1]][ct->numseq[j+1]];
				}
			}
	

			//after = 0 if a stacked pair cannot form 3' to i
			if ((((j-i)>minloop+2)&&(j<=ct->numofbases)||(j>ct->numofbases+1))&&(i!=ct->numofbases)) {
				after = inc[ct->numseq[i+1]][ct->numseq[j-1]];

			}
			else after = 0;

			adjacentgu = false;
			//check for preceding or following GU or whether the pair itself is gu
			if (ct->numseq[i+1]==3&&ct->numseq[j-1]==4) adjacentgu = true;
			else if (ct->numseq[i+1]==4&&ct->numseq[j-1]==3) adjacentgu = true;
			else if (ct->numseq[i]==3&&ct->numseq[j]==4) adjacentgu = true;
			else if (ct->numseq[i]==4&&ct->numseq[j]==3) adjacentgu = true;
			else if (i-1>0&&j+1<=ct->numofbases) {
				if (ct->numseq[i-1]==3&&ct->numseq[j+1]==4) adjacentgu = true;
				else if (ct->numseq[i-1]==4&&ct->numseq[j+1]==3) adjacentgu = true;

			}

			//if there are no stackable pairs to i.j then don't allow a pair i,j
			if ((before!=0)||(after!=0)) {
				
				if (i+1<j-1&&!adjacentgu) interior = erg1(i,j,i+1,j-1,ct,data)*v->f(i+1,j-1);
				else interior = 0;
				if (j+1<=ct->numofbases&&!adjacentgu) exterior = erg1(j,i+ct->numofbases,j+1,i+ct->numofbases-1,ct,data)*v->f(j+1,i+ct->numofbases-1);
				else exterior = 0;
				return ((v->f(i,j)+interior)*(v->f(j,i+ct->numofbases)+exterior)-interior*exterior)/(w5[ct->numofbases]*scaling*scaling);
			}
			else return 0;
				

		}
		else return 0;

	}

}

//function to rescale all arrays when partition function calculation is headed out of bounds
void rescale(int i, int j, structure *ct, pfdatatable *data, pfunctionclass *v, pfunctionclass *w, pfunctionclass *wl, pfunctionclass *wcoax,
			 pfunctionclass *wmb,pfunctionclass *wmbl, PFPRECISION *w5, PFPRECISION *w3, PFPRECISION *wca, PFPRECISION rescalefactor) {
	int ii,jj,lowi,nucs,index;
	double multiplier;

	for (jj=1;jj<=j;jj++) {

		if (jj==j) lowi=i;
		else if (jj<=ct->numofbases) lowi = 1;
		else lowi = jj-ct->numofbases+1+minloop;
		for (ii=min(jj,ct->numofbases);ii>=lowi;ii--) {
			//rescale v,w,wl,wcoax,wmb,wmbl
			
			//if (jj>ct->numofbases) nucs = ii + jj - ct->numofbases + 1;
			nucs = jj-ii+1;
			multiplier = pow(rescalefactor,(double) nucs);
				
				
			v->f(ii,jj)=v->f(ii,jj)*multiplier;
			w->f(ii,jj)=w->f(ii,jj)*multiplier;
			wl->f(ii,jj)=wl->f(ii,jj)*multiplier;
			wcoax->f(ii,jj)=wcoax->f(ii,jj)*multiplier;
			wmb->f(ii,jj)=wmb->f(ii,jj)*multiplier;
			wmbl->f(ii,jj)=wmbl->f(ii,jj)*multiplier;
			



			if (ii==(1)) {
				//rescale w5
				w5[jj]=w5[jj]*pow(rescalefactor,(double) jj);

				if (jj==ct->numofbases) {
					//rescale w3
					for (index=1;index<=ct->numofbases;index++) w3[index]=w3[index]*pow(rescalefactor,(double) (ct->numofbases-index+1));

				}
			}


		}
	}
	//for (ii=1; ii<=ct->numofbases;ii++)	cout <<"\nWCA: "<<wca[ii]<<"\n";
	for (ii=min(j,ct->numofbases);ii>=lowi;ii--) wca[ii]=wca[ii]*pow(rescalefactor,(double) (j-ii+1));
	//for (ii=1; ii<=ct->numofbases;ii++)	cout <<"\nWCA_rescaled: "<<wca[ii]<<"\n";
	data->rescaledatatable(rescalefactor);


}


//function to rescale all arrays when partition function calculation is headed out
	//of bounds when claculating w3
//void rescaleatw3(int ii, structure *ct, pfdatatable *data, pfunctionclass *v, pfunctionclass *w, pfunctionclass *wl, pfunctionclass *wcoax,
//				 pfunctionclass *wmb,pfunctionclass *wmbl, PFPRECISION *w5, PFPRECISION *w3, double rescalefactor) {

	 

//}

//function to rescale all arrays when partition function calculation is headed out
	//of bounds when claculating w5
//void rescaleatw5(int jj,structure *ct, pfdatatable *data, pfunctionclass *v, pfunctionclass *w, pfunctionclass *wl, pfunctionclass *wcoax,
//				 pfunctionclass *wmb,pfunctionclass *wmbl, PFPRECISION *w5, PFPRECISION *w3, double rescalefactor) { 

	//rescale the previously filled arrays
//	rescale(1, jj, ct, data, v, w, wl, wcoax, wmb, wmbl, w5, w3, rescalefactor);
	
//	w5[jj]=w5[jj]*pow(rescalefactor,(double) jj);

//}

//rescale the entries in datatable
void pfdatatable::rescaledatatable(PFPRECISION rescalefactor) {

	scaling=scaling*rescalefactor;
	int i,j,k,l,m,n,o,p;

//debug rescale
	cout <<"Rescaled now. Scaling: \t"<<scaling<<"\n";

	
	for (i=0;i<31;i++) {
		inter[i] = inter[i]*pow(rescalefactor,i+2);
		bulge[i] = bulge[i]*pow(rescalefactor,i+2);
		hairpin[i] = hairpin[i]*pow(rescalefactor,i+2);

	}
	for (i=0;i<6;i++) {
		for (j=0;j<6;j++) {
			for (k=0;k<6;k++) {
				for (l=0;l<3;l++) {
					dangle[i][j][k][l] = dangle[i][j][k][l]*rescalefactor;	
				}
				for (l=0;l<6;l++) {
					stack[i][j][k][l]=stack[i][j][k][l]*pow(rescalefactor,2);
					
					tstackcoax[i][j][k][l]=tstackcoax[i][j][k][l]*pow(rescalefactor,2);
					
					tstack[i][j][k][l]=tstack[i][j][k][l]*pow(rescalefactor,2);
					tstkm[i][j][k][l]=tstkm[i][j][k][l]*pow(rescalefactor,2);
					
					for (m=0;m<6;m++) {
						for (n=0;n<6;n++) {
							iloop11[i][j][k][l][m][n]=iloop11[i][j][k][l][m][n]*pow(rescalefactor,4);
							for (o=0;o<6;o++) {
								iloop21[i][j][k][l][m][n][o]=iloop21[i][j][k][l][m][n][o]*pow(rescalefactor,5);
								for (p=0;p<6;p++) {
									iloop22[i][j][k][l][m][n][o][p]=iloop22[i][j][k][l][m][n][o][p]*pow(rescalefactor,6);
								}
							}
							

						}
					}
				}
			}
		}
	}
	
	for (i=0;i<=numoftloops;i++) {
		
		tloop[i] = tloop[i]*pow(rescalefactor,6);

	}
	
	for (i=0;i<=numoftriloops;i++) {
		
		triloop[i] = triloop[i]*pow(rescalefactor,5);
	}
	
	for (i=0;i<=numofhexaloops;i++) {
		
		hexaloop[i] = hexaloop[i]*pow(rescalefactor,8);

	}
	


}


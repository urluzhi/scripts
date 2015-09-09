
#if !defined(STRUCTURE_H)
#define STRUCTURE_H




#include "defines.h"



//////////////////////////////////////////////////////////////////////
class structure //this structure contains all the info for a structure
{
public:
int numofstructures;//number of structures 
short int numofbases;//number of nucleotides in sequence
short int pair[maxforce][2],npair,nforbid,forbid[maxforce][2];
short int *numseq,*hnumber;
short int **basepr;
int *energy;//[maxstructures+1];
char **ctlabel;//[maxstructures+1][ctheaderlength];
short int ndbl, dbl[maxforce];
int inter[3],allocatedstructures;
short int nnopair,nopair[maxforce],nmod,mod[maxforce];
short int ngu,gu[maxgu];
char *nucs;
bool intermolecular,allocated,templated,stacking;
bool **tem;//tem stores template information as to whether a pair is allowed
structure(int structures = maxstructures+1);
~structure();
void allocate(int size = maxbases);
void allocatestructure();
void deletestructure();
void checknumberofstructures();//check to make sure there is room for one more structure
void allocatetem();//allocate space in **tem 
short int min_gu, min_g_or_u;//NMR-derived constraint variables
short int neighbors[maxforce][maxneighborlength],nneighbors;//also NMR-derived index this from zero in both dimensions
//regional NMR constraints:
short int nregion,rmin_gu[maxregions],rmin_g_or_u[maxregions];
short int rneighbors[maxregions][maxforce][maxneighborlength],rnneighbors[maxregions],start[maxregions],stop[maxregions];
//microarray type constraints:
short int nmicroarray,microstart[maxregions],microstop[maxregions],microunpair[maxregions];
bool *fcedbl;//pointer to a 2-D array used in Dynalign to track nucleotides that must be double-stranded
void sort();//sort the structures so that the lowest free energy structure is in position
void ReadSHAPE(char *filename, float SingleStrandThreshold, float ModificationThreshold);//Read SHAPE reactivity data from a file
void ReadSHAPE(char *filename);//Read SHAPE reactivity data from a file
bool limitdistance;//toggle to indicate that there is a limit on the maximum distance between nucs in base pairs
int maxdistance;//maximum distance between nucs in base pairs
double *SHAPE;//double array to contain SHAPE data -- values less than -500 are ignored
bool shaped;//keeps track of whether SHAPE data was loaded
double SHAPEslope,SHAPEintercept;//values of slope and intercept for SHAPE data modification of pairing stability

/*structure is set up to hold many possible structures of the same sequence
	numofbases = number of bases in sequence
	numofstructures = number of alternative structures of the sequence
				that is held by structure
	numseq[i] = a numeric that stands for the base in the ith position
				of the sequence,
			A = 1
			C = 2
			G = 3
			U = 4
	basepr[i][j] = base to which the jth base is paired in the ith structure
	force[i] = any information about the way the ith pair is treated when
				folding; eg: forced single, etc
	energy[i] = 10 * the Gibb's free energy of the ith structure, this
				is done so that integer math can be performed
				and the nearest tenth of a kcal/mol can be
				followed
	ctlabel = a string of information for each of the structures
	fce = an array that can keep track of how each i,j pair is being
				forced
   hnumber array stores the historical numbering of a sequence
   nucs is a character array to store the sequence information --
   	this will allow the program to keep T and U from getting confused

	stacking = a bool to keep track of whether stacking is being tracked in the ct.
		If stacking is to be tracked, stacking must be manually set to true before allocate
		is called.  The stacks are stored in basepr, by doubling the size of basepr, in basepr[index][i+N],
		where N is the number of nucleotides.  This lets the stacking info get sorted along with the
		pairing info in ::sort.
		When stacking is true, basepr[i][j]!=implies a stack for nucleotide j-numofbases.
		For j-numofbases unpaired and j-numofbases is not indicated in another stack, this nucleotide (j-numofbases) 
		is either half a terminal mismatch or a dangling end that is stacked on the indicated nucleotide.
		For j-numofbases paired, this is a coaxial stack in which j-numofbases is stacked onto the 
		nucleotide as indicated.  For flush stacks, both nucleotides in terminal basepairs indicate each other as
		stacks.  For stacks with an intervening mismatch, one helix will point to an unpaired nuc that mediates the stack
		and then that unpaired nucleotide will indicate the next helix in the coaxial stack.

	
*/



};

short int tonumi(char *base); //converts base to a numeric


char *tobase (int i);//convert a numeric value for a base to the familiar
								//character
void sortstructures (structure *ct);//this routine resorts the structures
					//predicted by dynamic according to energies calculated by efn2

void tonum(char *base,structure *ct,int count); //converts base to a numeric

void openct(structure *ct,char *ctfile);//reads ct file

int openseq (structure *ct,const char *seqfile);//inputs a sequence from file
															// seqfile
void ctout (structure *ct,const char *ctoutfile);//outputs a ct file


int ecompare(const void *i, const void *j);

void swap(int *a,int *b);//Swap two variables
void swap(short int *a,short int *b);//swap two variables
void swap(float *a,float *b);//swap two variables
#endif

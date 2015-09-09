/* This program finds the secondary RNA structure that has the maximum number of hydrogen bonds
 between pairing bases.
 
 The number of hydrogen bonds between allowed pairing bases.
 
 Allowed base pairing |  # of hydrogen bonds
 A --- U        |         2
 C --- G        |         3
 U --- G        |         2
 Author: Xinqiang Ding
 Version: 2012.8.9
 Bug report: dingxinqiang@gmail.com
 */

#include <iostream>
#include <cstring>
#include <cstdio>
#include <ctype.h>
using namespace std;

int NumHydrogenBonds(char baseOne, char baseTwo);
void BracketNotation(int **P, int i, int j);

int main(int argc, char *argv[]){
	int size; // the length of RNA.
	size = strlen(argv[1]);
	char *seq;
	int i; // the firt base index of a RNA segment.
	int j; // the end base index of a RNA segment.
	int l; // the base index between i and j, so i<=l<j.
	int k;  // the length of RNA between base i and base j, so k = j-i+1.
	int V[size][size]; /* V array: V(i,j) is the maximum hydrogen bond between base i and base j
                        on the condition that base i and base j are pairing.
                        */
	int W[size][size]; // W array: W(i,j) is the maximum hydrogen bonds between base i and base j.
	
    //int P[size][size];
    /* P(i,j) is the k so that W(i,j) = W(i,k) + W(k+1,j) and i <= k < j,
     but if W(i,j) = V(i,j), then P(i,j) = -2; if j-i+1<=4, then P(i,j)=-1.
     */
	int **P = new int * [size];
	for (i=0;i<size;i++){
		P[i] = new int [size];
	}
    seq = argv[1];
    // Fill the V, W and P array
	for (k = 1; k <= size; k++){
		if (k<=4){
			for (i = 0; i <= size-k; i++){
				V[i][i+k-1] = 0;
				W[i][i+k-1] = 0;
				P[i][i+k-1] = -1;
			}
		}
		else{
			for (i = 0; i <= size-k; i++){
				if (NumHydrogenBonds(seq[i],seq[i+k-1])==0)
                    V[i][i+k-1] = 0;
				else
                    V[i][i+k-1] = NumHydrogenBonds(seq[i],seq[i+k-1]) + W[i+1][i+k-2];
                W[i][i+k-1] = V[i][i+k-1];
                P[i][i+k-1] = -2;
				for (l=i;l<i+k-1;l++){
					if (W[i][i+k-1] < (W[i][l] + W[l+1][i+k-1])){
						W[i][i+k-1] = W[i][l] + W[l+1][i+k-1];
						P[i][i+k-1] = l;
					}
				}
                if (W[i][i+k-1] == 0) {
                    P[i][i+k-1] = -3;
                }
			}
		}
	}
    
    // Print the V, W and P array
  	cout<< "V array:\n";  
	for (j = size-1; j >= 0; j--){
		for (i = 0; i <=j ; i++){
			cout << V[i][j] << ' ';
		}
		cout << "\n";
	}
    
	cout << "\n\n";
  	cout<< "W array:\n";  
	for (j = size-1; j >= 0; j--){
		for (i = 0; i <=j ; i++){
			cout << W[i][j] << ' ';
		}
		cout << "\n";
	}
	cout << "\n\n";
   /* 
	for (j = size-1; j >= 0; j--){
		for (i = 0; i <=j ; i++){
			cout << P[i][j] << ' ';
		}
		cout << '\n';
	}
*/
	cout << "Predicted structure:\n";
	cout << seq << '\n';
	BracketNotation(P, 0, size-1);
	cout << '\n';

	return 0;
}


/*
 define the function NumHydrogenBonds to return the
 number of hydrogen bonds between baseOne and baseTwo
 */
int NumHydrogenBonds(char baseOne, char baseTwo){
    baseOne = toupper(baseOne);
	baseTwo = toupper(baseTwo);
	switch(baseOne){
		case 'A':{
			switch (baseTwo){
				case 'A': return 0;
				case 'U': return 2;
				case 'G': return 0;
				case 'C': return 0;
			}
		}
		case 'U':{
			switch (baseTwo){
				case 'A': return 2;
				case 'U': return 0;
				case 'G': return 2;
				case 'C': return 0;
			}
		}
		case 'C':{
			switch (baseTwo){
				case 'A': return 0;
				case 'U': return 0;
				case 'G': return 3;
				case 'C': return 0;
			}
		}
		case 'G':{
			switch (baseTwo){
				case 'A': return 0;
				case 'U': return 2;
				case 'G': return 0;
				case 'C': return 3;
			}
		}
    }
}

/*
 Print the bracket notation for RNA secondary structures
 */

void BracketNotation(int **P, int i, int j){
	if (i<j){
		if(P[i][j]==-2){
			cout << '(';
			BracketNotation(P,i+1,j-1);
			cout << ')';
		}
		else if(P[i][j]==-1){
			cout << '.';
			BracketNotation(P,i+1,j-1);
			cout << '.';
		}
        else if(P[i][j] == -3){
            cout << '.';
            BracketNotation(P, i+1, j-1);
            cout << '.';
        }
		else{
			BracketNotation(P,i,P[i][j]);
			BracketNotation(P,P[i][j]+1,j);
		}
	}
	else if (i==j){
		cout << '.';
		return;
	}
}



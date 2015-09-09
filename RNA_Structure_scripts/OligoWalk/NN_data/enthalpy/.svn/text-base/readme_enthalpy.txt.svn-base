Updated: Oct. 31, 2006

There are nineteen files of enthalpies for RNA-RNA interactions.
The derivation of these parameters is in
"Zhi John Lu, Douglas H. Turner, and David H. Mathews
A set of nearest neighbor parameters for predicting the enthalpy change
of RNA secondary structure formation
Nucleic Acids Res., October 2006; 34: 4912 - 4924."

Here is a brief summary of what the files contain:

Stack.dh contains the nearest neighbor enthalpy parameters for stacking in a helix.  

Dangle.dh presents the enthalpy of dangling ends (5 prime or 3 prime). 
These enthalpies are used in multibranch and exterior loops.  

The loop.dh table contains the initial penalty for closing hairpin, 
internal, and bulge loops. Each has its own penalty based on length.  

Tstackh.dh gives the enthalpy of stacking for a terminal mismatch in a 
hairpin loop.
Tstacki.dh gives the enthalpy of stacking for a terminal mismatch in an internal loop.  
Tstacki23.dh is used for terminal mismatch stacking in internal loops of size 2x3 unpaired nucleotides (and 3x2).
Tstacki1n.dh is used for terminal mismatch stacking in internal loops of size 1xn nucleotides where n>2. 
Tstackm.dh is used for terminal stacking in a multibranch loop.
Tstack.dh is for terminal mismatch stacking in exterior loops, i.e. loops that contain the ends of the sequence.  

Miscloop.dh contains miscellaneous loop parameters.

Int21.dh is a lookup table for the enthalpy of a 2x1 internal loop.  The values presented are  the total stability of the loop, including the initial penalty and the terminal stack bonus.

Int22.dh is a lookup table for the enthalpy of a 2x2 internal loop analogous to the 2x1 internal loop.  
Int11.dh is the lookup table for single mismatches.

Tloop.dh is the tetraloop bonus table. For hairpin loops with four unpaired nucleotides,
this table is consulted. If the sequence (starting with the 5 prime last paired nucleotide
and finishing with its 3 prime paired nucleotide) appears in the table, the 
bonus is applied to the hairpin's stability.  

Triloop.dh is a lookup table for hairpin loops of three.  It is applied analogously to the tetraloop 
table (above).  

Also, hexaloop.dh is a lookup table for hairpins of six 
nucleotides.  

Coaxial.dh is the enthalpy of coaxial stacking for two 
helixes with no intervening unpaired nucleotides.  When applied, the top 
nucleotides are the continuous sequence and the bottom nucleotides are 
interrupted by the rest of the multibranch loop or exterior loop.  

Coaxstack.dh is one of the two tables for coaxial stacking with an 
intervening mismatch.  This is the stack with the open backbone.  

Tstackcoax.dh is the second of the tables for coaxial stacking with an 
intervening mismatch.  This is the stack with the continuous backbone.

 

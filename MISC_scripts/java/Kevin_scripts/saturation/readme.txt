Author	Kevin Yuk-Lap Yip, Gerstein Lab, Yale University
	(project with Zhi John Lu)
Version	1.2 (March 5, 2010)

This program produces saturation plots from a set of files, each containing a
list of genomic regions. Each genomic region is specified using the following
format:

<ID><tab><start><tab><end>

where
<ID>	is the identifier of the region-at-large, such as the chromosome
<start>	is the starting position of the region
<end>	is the ending position of the region (this position is inside the
	region)

The y-axis could be the absolute number of nucleotides, or a fraction of an
input total number of nucleotides, such as the total number of nucleotides of
the coding transcripts in the example. To use the absolute number, input the
total as 0.

If the number of input files is no more than 31, the program can compute the
coverage from all combinations of the input datasets. If the number of input
files is more than 31, or if the number of combinations is more than a specified
threshold, a random sample of the combinations will be considered.

To run the program, use the example script file (run.bat or run.sh), or type
"java -classpath . SaturationPlotsCreator" for a description of the usage.

The program requires a Java virtual machine or a Java runtime environment
installed, preferrably with J2SE 1.6.0 or above. Visit
http://java.sun.com/javase/downloads/index.jsp for instructions.

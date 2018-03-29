#!/usr/bin/perl 
use strict;
use warnings;
use Bio::SeqIO;

#
#############################################################################
#print out the usage inforagion
my $usage = <<EOF;

USAGE:
	seqmaniuplator.pl [options] (PATTERN) [FILE1 FILE2 ...] 
DESCRIPTION:
	This is a RNA sequence manipulating package, it incudes some useful tools and subroutines to handle .seq file of RNA/DNA.
OPTIONS:
	-h,--help
		Usage information
	-translate [DNA.fasta]
		translate DNA to Protein (have to start with ATG)
	-rev_translate [Protein.fasta]
		reverse translate Protein to RNA by randomly picked genetic codon
	-fetch [FILE.fasta] chr start end strand(+/-)
		fetch sequence from start to end
	-mfetch [FILE.fasta] <FILE.loc>
		fetch sequences from start to end in FILE.loc
	-sp [FILE.fasta]
		split fasta files
	-fasta2seq [FILE.fasta] [FILE.seq]
		transform fasta file to .seq file
	-check [FILE.seq]
		Read and check the .seq file 
	-wseq [FILE.seq] NUCS_STRING; -wseq [FILE.seq] <STDIN>
		write nucleotides into a .seq file
	-plain2fasta [FILE.txt] [FILE.fasta]
		transform plain sequence file to .fasta file
	-seq2fasta [FILE.seq] [FILE.fasta]
		transform .seq file to .fasta file
	-gb2seq [FILE.gb] [FILE.seq]
		transform Genebank database file(.gb) to .seq file
	-align2maf
		transform strecher alignment to maf format
	-count [FILE.seq] ; -count <STDIN>
		Count the number of each nucleotide in a sequence
	-count_fasta [FILE.fasta]
		Count the Chromsome number and lengths
	-mylogo [SEQ_COUNT] [OUTPUT_1000SEQ] Length (nocount)
		generate logo seq for web logo (if "nocount": dup. not counted)
	-freq [FILE.seq] ; -freq <STDIN>
		Calculate the percentage frequence of each nucleotide in a sequence 
	-anti [FILE.seq] ; -anti <STDIN>    
		Generate antisense strand of a sequence
	-match_motif [ac]gaa[ct] [fasta]
		Match motif in fasta file

SUBROUTINES:
	&read_seq 
		my (\$seqname, \$nucs)= &read_seq("FILE.seq")
		read the .seq file into a string
	&write_seq
		write("FILE.seq",\$nucs)

Created: Sep. 28, 2006      Last Modified:  Sep. 29, 2007
by John(Zhi Lu)            zhi_lu\@urmc.rochester.edu
EOF

#outpust usage information
die $usage unless (@ARGV && $ARGV[0] ne "-h" && $ARGV[0] ne "--help");

#different usages
if ($ARGV[0] eq "-check" && $ARGV[1]) { #read .seq file
	my ($seqname,$nucs)= read_seq($ARGV[1]);	
#	print "$seqname :\t $nucs\n";	
	print "$seqname :\t Correct\n";	
}
if ($ARGV[0] eq "-wseq" && $ARGV[1]) { #read .seq file
	my $seqname=$ARGV[1];
	my $nucs;
	if ($ARGV[2]) 	{	$nucs=$ARGV[2];	}
	else			{	$nucs=<STDIN>;	}
	write_seq($seqname,$nucs);	
}
elsif ($ARGV[0] eq "-fasta2seq" && $ARGV[1] && $ARGV[2]) { #transform .fasta file
	my $nucs= read_fasta($ARGV[1]);
	write_seq($ARGV[2],$nucs);	
#print "\nDone: Transform sequence from $ARGV[1] to $ARGV[2]\n\n";
}
elsif ($ARGV[0] eq "-fetch" && $ARGV[1] && $ARGV[2]) { #transform .fasta file
	&fetch;	
}
elsif ($ARGV[0] eq "-mfetch" && $ARGV[1] ) { #transform .fasta file
	&mfetch;	
}
elsif ($ARGV[0] eq "-translate" && $ARGV[1] ) { #transform .fasta file
	&translate;	
}
elsif ($ARGV[0] eq "-rev_translate" && $ARGV[1] ) { #transform .fasta file
	&rev_translate;	
}

elsif ($ARGV[0] eq "-match_motif" && $ARGV[1] ) { #transform .fasta file
	&match_motif;	
}

elsif ($ARGV[0] eq "-seq2fasta" && $ARGV[1] && $ARGV[2]) { #transform .fasta file
	my ($seqname,$nucs)= read_seq($ARGV[1]);
	write_fasta($ARGV[2],$seqname,$nucs);	
	print ".\n";
}
elsif ($ARGV[0] eq "-plain2fasta" && $ARGV[1] && $ARGV[2]) { #transform .fasta file
	&plain2fasta;
}


elsif ($ARGV[0] eq "-gb2seq" && $ARGV[1] && $ARGV[2]) { #transform .gb file
	my $nucs= read_gb($ARGV[1]);
	write_seq($ARGV[2],$nucs);	
	print "\nDone: Transform sequence from $ARGV[1] to $ARGV[2]\n\n";
}
elsif ($ARGV[0] eq "-align2maf" && $ARGV[1] && $ARGV[2]) { #transform to.maf file
	my ($align1, $align2) = read_align($ARGV[1]);
	write_maf($ARGV[2],$align1,$align2);	
}
elsif ($ARGV[0] eq "-count") {#count the nucleotides' numbers
	my $nucs= &parse_input;
	count($nucs,0);
}
elsif ($ARGV[0] eq "-count_fasta") {#count the nucleotides' numbers
	&count_fasta;
}
elsif ($ARGV[0] eq "-mylogo") {#count the nucleotides' numbers
	&mylogo;
}

elsif ($ARGV[0] eq "-freq") {#count the nucleotides' frequences 
	my $nucs= &parse_input;
	count($nucs,1);
}
elsif ($ARGV[0] eq "-sp") {#count the nucleotides' frequences 
	&spfasta;
}
elsif ($ARGV[0] eq "-anti") {#antisense strand of a seqence
	my $nucs= &parse_input;
	my $compliment= anti_sense($nucs);
	print "$compliment\n";

}
else {
	print "\n\t\tUnrecognized options!!!\n\n\n";
	die $usage;
}
exit;
#subroutines:
##############################################################
#parse the input
sub parse_input {
	my $nucs;
	if ($ARGV[1])	{#read a .seq file	
		my $seqname;	
		($seqname,$nucs)= read_seq($ARGV[1]);	
	}
	else {	#read from a standard input
		$nucs=<STDIN>;		
		$nucs=~ s/[\s]//g;
	}
	return $nucs;
}

#read the nucleotides in a string from a .seq file
sub read_seq {
	my ($seq)=@_;
	my ($seqname,$nucs);
	open (SEQ,$seq) or die("Cannot open $seq file $!");
	my $linei=0;
	while (<SEQ>) {
		my $line= $_;
		#discard the ; line
	    if ($line =~ /;/ ) {
    	    $linei=1;
        	next;
    	}
	    #discard all the lines before ; line
    	if ($linei==0)  {
        	next;
    	}
    	#discard the first line after ; line
	    if ($linei==1) {
    	    $linei++;
        	$seqname=$line;		#read the sequence name from the line
        	$seqname=~ s/\n//;		#read the sequence name from the line
			next;
   		 }
    	$nucs .=$line;
	}
	#clean some white spaces
	$nucs=~ s/[\s1]//g;
	$nucs=~ s/\cM//g;
	#checking the format of nucleotides
	my $checknucs=$nucs;
	$checknucs=~ s/[AGTCUagtcuxXN]//g;
	&die("Sequence file has wrong format,codes other than AGTCUXNagtcux were found\n")	unless($checknucs eq '');
	close SEQ;
	return ($seqname,$nucs);
}

#calculate the frequencies of each nucleotides in a sequence
sub count {
	#readthe sequences
	my ($nucs,$isfreq) = @_;
	#transfomat to the standard nuc. format
	$nucs=~ tr/acgtuxn/ACGTUXN/;
	#Initiate the counts
	my $count_of_A = 0;
	my $count_of_G = 0;
	my $count_of_C = 0;
	my $count_of_U = 0;
	my $count_of_T = 0;
	my $count_of_X = 0;
	my $count_of_N = 0;
	my $count_of_err= 0;
	#count the bases
	for (my $position=0; $position < length $nucs; ++$position){
		my $base = substr($nucs, $position,1);
		if ($base eq 'A') {
			++$count_of_A;
		}	
		elsif ($base eq 'G') {
			++$count_of_G;
		}	
		elsif ($base eq 'C') {
			++$count_of_C;
		}	
		elsif ($base eq 'U') {
			++$count_of_U;
		}	
		elsif ($base eq 'T') {
			++$count_of_T;
		}	
		elsif ($base eq 'N') {
			++$count_of_N;
		}	
		elsif ($base eq 'X') {
			++$count_of_X;
		}	
		else {
			++$count_of_err;
		}	
	}
	#report the error counts
	if($count_of_err != 0)	{	die ("$count_of_err non-base(s) was/were found in  \n $nucs");	}
	#oupt the counts
	my	$all_counts=length $nucs;
	print "length= $all_counts\t";
	if ($isfreq) {	#calculate the frequency
		my	$freq_of_A = $count_of_A/$all_counts;
		my	$freq_of_G = $count_of_G/$all_counts;
		my	$freq_of_C = $count_of_C/$all_counts;
		my	$freq_of_U = $count_of_U/$all_counts;
		my	$freq_of_T = $count_of_T/$all_counts;
		my	$freq_of_N = $count_of_N/$all_counts;
		my	$freq_of_X = $count_of_X/$all_counts;
		printf "A= %.2f\t", $freq_of_A;
		printf "C= %.2f\t", $freq_of_C;
		printf "G= %.2f\t", $freq_of_G;
		printf "U= %.2f\t", $freq_of_U;
		printf "T= %.2f\t", $freq_of_T;
		if ($count_of_X != 0)	{	printf "X= %.2f\t", $freq_of_X;	}
		if ($count_of_N != 0)	{	printf "X= %.2f\t", $freq_of_N;	}
		print "\n";
	}
	else { #output the counts
		print "A= $count_of_A\t";
		print "C= $count_of_C\t";
		print "G= $count_of_G\t";
		print "U= $count_of_U\t";
		print "T= $count_of_T\t";
		if($count_of_X !=0)	{	print "X= $count_of_X\t";	}
		if($count_of_N !=0)	{	print "N= $count_of_N\t";	}
		print "\n";
	}
}
#calculate the frequencies of each nucleotides in a sequence
sub mylogo {
	my @pos; 
	#readthe sequences
	open (TABLE,"$ARGV[1]" ) or  &die("Can't open file: $ARGV[1] $!\n");
	while (<TABLE>) {
		my @t=split;
		# get nuc from each colum
		for (my $i=0;$i<=$ARGV[3];$i++) {
			my $count=1;
			unless ($ARGV[4])	{	$count=$t[1];	}
			for (my $j=0;$j<$count;$j++) {
				$pos[$i] .= substr($t[0],$i,1);	
			}
		}
	}
	close TABLE;
	#clulate freq. for each colum
	my $total= length($pos[0]);
	my @loc; my $j;
	print "Pos.\tA\tT\tG\tC\n";
	for (my $i=0;$i<$ARGV[3];$i++){
		#		print $pos[$i],"\n";
		my $freq_of_A= round (($pos[$i] =~ tr/A//)*1000/$total);
		my $freq_of_T= round(($pos[$i] =~ tr/T//)*1000/$total);
		my $freq_of_G= round(($pos[$i] =~ tr/G//)*1000/$total);
#		my $freq_of_C= ($pos[$i] =~ tr/C//);
		my $freq_of_C= 1000-$freq_of_A-$freq_of_T-$freq_of_G;
		print $i+1,"\t",$freq_of_A/10,"\t",$freq_of_T/10,"\t",$freq_of_G/10,"\t",$freq_of_C/10,"\n";
		for ($j=0;$j<$freq_of_A;$j++) {	$loc[$i] .= 'A'; 	}
		for ($j=0;$j<$freq_of_T;$j++) {	$loc[$i] .= 'T'; 	}
		for ($j=0;$j<$freq_of_G;$j++) {	$loc[$i] .= 'G'; 	}
		for ($j=0;$j<$freq_of_C;$j++) {	$loc[$i] .= 'C'; 	}
	}
	#output 1000 samples
	open (OUT, ">$ARGV[2]");	
	for (my $j=0;$j<1000;$j++) {
		for (my $i=0;$i<$ARGV[3];$i++){	print OUT  substr($loc[$i],$j,1);	}
		print OUT "\n";
	}	
	close OUT;
}
sub round {
	    my($number) = shift;
		#considering negative number as well	
		return int($number + .5 * ($number <=> 0));
}

#transcribe the sequence to the complementary strand
sub anti_sense {
		my ($seq) = @_;
		my $anti_seq="";
		for (my $pos=0;$pos<length($seq);$pos++) {
			$anti_seq= substr($seq,$pos,1).$anti_seq;
		}
		$anti_seq =~ tr/ACGTUacgtu/UGCAAugcaa/;
		return "$anti_seq\n";
}
sub anti_sense_dna {
		my ($seq) = @_;
		my $anti_seq="";
		for (my $pos=0;$pos<length($seq);$pos++) {
			$anti_seq= substr($seq,$pos,1).$anti_seq;
		}
		$anti_seq =~ tr/ACGTUacgtu/TGCAAtgcaa/;
		return "$anti_seq";
}


#read sequence from GenBank file to a string
sub read_gb {
	my ($gb)=@_;
	#open the GenBack file
	open (GB,"$gb" ) or  &die("Can't open file: $gb $!\n");
	my $save_input_separator=$/;
	#Set input separator and read in a record to a scalar
	$/ = "//\n"; 
	my $record = <GB>;
	#Set the separator back to default
	$/ = $save_input_separator;  
	#Now fetch the sequence data using regular expression
#	my ($sequence) = ($record =~ /^LOCUS.*ORIGIN.*\n\s*(1\s.*)\/\/\n/s);
	$record =~ /^ORIGIN.*\n((.*\n)*)\/\/\n/m;
	my $sequence = $1;	
	unless ($sequence) {
		print "$gb is not parsed correctly!!!!!!!!!!!!\n";
		exit;
	}
	#clean the sequence
	$sequence =~ s/[\d\s\/]//g;
#	$sequence =~ tr/T/U/;	#transcribe DNA to RNA nt 
	close GB;
	return $sequence;
}

#fetch nucs  from start to end on +/- strand
sub fetch {
	my ($fastafile)=$ARGV[1];
	#open the fasta file
	open (FASTA,"$fastafile" ) or  die("Can't open file: $fastafile $!\n");
	my $sequence='';
	my $foundit=0;
	while (<FASTA>){
		#discard blank line
		if ($_ =~ /^\s*$/)		{	next;	}
		#discard comment line
		elsif ($_ =~ /^\s*#/)	{	next;	}
		#discard the header line
		elsif ($_ =~ /^>$ARGV[2]/) 		{ $foundit=1;	next;	}
		#keep line
		if ($foundit) {
			if ( $_ =~ /^>/) { last;}
			if ($_ =~ /^\s*$/)	{	next;	}
			else {$sequence .= $_;}
		}
	}
	close FASTA;
	$sequence =~ s/\n//g;
	print ">$ARGV[2]_$ARGV[3]_$ARGV[4]_$ARGV[5]\n";
	my $out= substr($sequence, $ARGV[3]-1,$ARGV[4]-$ARGV[3]+1);
	if ($ARGV[5] eq '+') {
		print $out,"\n";
	}
	else {
		print anti_sense_dna($out),"\n";
	}
}
#mfetch nucs  from start to end on +/- strand
sub mfetch {
	#read the coordinates of each gene
	my %hash;	my %sequences;
	my @input=<STDIN>;
	foreach (@input) {
		my @t=split;
		my $chr=$t[0];
		$chr =~ s/CHROMOSOME_//i;
		$chr =~ s/chr//i;
		my $name=$chr."_".$t[1]."_".$t[2]."_".$t[3];
		$hash{$name}=$chr;
	}
	#read the genome (fasta)
    my $seq_in  = Bio::SeqIO->new( -format => 'Fasta', -file => $ARGV[1]);
	while(my $seq = $seq_in->next_seq() ) {
		my $id=$seq->id;
		$id =~ s/CHROMOSOME_//i;
		$id =~ s/chr//i;
		my $chr_l=$seq->length;
		foreach (keys %hash) {
			if ($hash{$_} eq $id ) {
				my $pos=$_;
				my @co=split /_/,$pos;
				my $start;my $end; my $strand;
				$start=$co[1];
				$end=$co[2];
				if($co[3]) {$strand=$co[3];}
				my $out=$seq->subseq($start,$end);
				if (!$co[3]||$strand eq '+') {
					$sequences{$pos}=$out;
				}
				else {
					$sequences{$pos}=anti_sense_dna($out);
				}
			}
		}
	}
	my $i=1;
	foreach (@input) {
		my @t=split;
		my $chr=$t[0];
		$chr =~ s/CHROMOSOME_//i;
		$chr =~ s/chr//i;
		my $name=$chr."_".$t[1]."_".$t[2]."_".$t[3];
		print ">seq",$i++,"\t";
		print "$name\n";
		print $sequences{$name};
		print "\n";
	}
}
#translate DNA to Protein
sub translate {
    my $seq_in  = Bio::SeqIO->new( -format => 'Fasta', -file => $ARGV[1]);
	my $seq_out = Bio::SeqIO->new( -format => 'Fasta');
	while(my $seq = $seq_in->next_seq() ) {
		if ($seq->subseq(1,3) eq 'ATG') {
			$seq_out->write_seq($seq->translate());
		}
	}
}

#reverse translate Protein to DNA
sub rev_translate {
	my %codons=&read_codon;
    my $seq_in  = Bio::SeqIO->new( -format => 'Fasta', -file => $ARGV[1]);
	my $seq_out = Bio::SeqIO->new( -format => 'Fasta');
	while(my $seq = $seq_in->next_seq() ) {
		my $sequences='';
		for (my $i=0;$i<$seq->length;$i++) {
			my $aa=substr($seq->seq,$i,1);		
			fisher_yates_shuffle($codons{$aa});
			$sequences.=$codons{$aa}[0];
		}
		my $myseq = Bio::Seq->new( -seq => $sequences, -id  => $seq->id);
		$seq_out->write_seq($myseq);
	}
}
#match motif
sub match_motif{
    my $seq_in  = Bio::SeqIO->new( -format => 'Fasta', -file => $ARGV[2]);
	while(my $seq = $seq_in->next_seq() ) {
		my $id=$seq->id;
		my $sequence=$seq->seq;
		while($sequence =~/$ARGV[1]/ig) {
			print $id,"\t-",length($sequence)-length($`),"\t+1\n";
		}
		my $anti=anti_sense_dna($sequence);
		while($anti =~/$ARGV[1]/ig) {
			print $id,"\t-",length($`)+1,"\t-1\n";
		}
	}
}

#read sequence from fasta file to a string
sub read_fasta {
	my ($fastafile)=@_;
	#open the GenBack file
	open (FASTA,"$fastafile" ) or  die("Can't open file: $fastafile $!\n");
	my $sequence='';
	while (<FASTA>){
		#discard blank line
		if ($_ =~ /^\s*$/)		{	next;	}
		#discard comment line
		elsif ($_ =~ /^\s*#/)	{	next;	}
		#discard the header line
		elsif ($_ =~ /^>/) 		{	next;	}
		#keep line
		else{
			$sequence .= $_;
		}
	}
	close FASTA;
	return $sequence;
}

#output .seq file
sub write_seq {
	my ($seqname, $sequence) = @_;
	#clean the sequence
	$sequence =~ s/[\d\s\/]//g;	
	$sequence =~ tr/acgtu/ACGTU/;
	$sequence =~ tr/T/U/;
	#write the .seq file
	open (SEQ, ">$seqname") or &die("Can't write file: $seqname $!\n");
	print SEQ "; ", length($sequence)," bases\n";
	$seqname=~ s/\.seq//;
	my @seq_name=split "\/",$seqname;
	print SEQ "$seq_name[scalar @seq_name -1]\n";
	print SEQ "$sequence";
	print SEQ "1\n";
	close SEQ;
}

#read sequence from fasta file to a string
sub write_fasta {
	my ($fastafile,$seqname,$nucs)=@_;
	#open the GenBack file
	open (FASTA,">$fastafile" ) or  die("Can't open file: $fastafile $!\n");
	print FASTA ">1\n";
	print FASTA "$nucs\n";

	close FASTA;
}

sub read_align {
	my ($file)=@_;
	open (FILE,"$file" ) or  die("Can't open file: $file $!\n");
	my $align1='';
	my $align2='';
	my $first=1;
	while (<FILE>){
		#discard blank line
		if ($_ =~ /^\s*$/)		{	next;	}
		#discard comment line
		elsif ($_ =~ /^#/)	{	next;	}
		#keep line
		elsif ($_ =~ /[AUGCTN]/) { 
			if ($first) {
				$align1 .= $_;
				$first= 0 ;
			}
			else {
				$align2 .=$_;
				$first=1;
			}
		}
	}
	$align1 =~ s/fasta\///g;
	$align1 =~ s/\s//g;
	$align1 =~ s/\d//g;

	$align2 =~ s/fasta\///g;
	$align2 =~ s/\s//g;
	$align2 =~ s/\d//g;
	close FILE;

	return ($align1,$align2);
}

#output .seq file
sub write_maf {
	my ($maf, $align1,$align2) = @_;
	#clean the sequence
	$align1 =~ tr/-acgtun/.ACGTUN/;
	$align2 =~ tr/-acgtun/.ACGTUN/;
	#change to RNA
	$align1=~ s/T/U/g;
	$align2=~ s/T/U/g;
	open(MAF,">$maf" ) or  &die("Can't write file: $maf $!\n");
	$maf=~ s/\.maf//;
	print MAF "##maf version=1 scoring=multiz.v7; $maf\n";
	print MAF "a score=99999\n";
    print MAF "s seq1\t 9999  99   +   99999\t $align1\n";
   	print MAF "s seq2\t 9999  99   +   99999\t $align2\n";
   	close MAF;
}

#plain 2 fasta file
sub plain2fasta {
	#read nucleotides from plain text
	open(TXT,"$ARGV[1]" ) or  &die("Can't write file: $ARGV[1] $!\n");
	#make  the FASTA file
	open (FASTA,">$ARGV[2]" ) or  die("Can't open file: $ARGV[2] $!\n");
	my $i=1;
	while (<TXT>){
		print FASTA ">$i\n";
		print FASTA "$_";
		$i++;	
	}
	close FASTA;
}
sub spfasta {
	my $seq_in  = Bio::SeqIO->new( -format => 'Fasta', -file => $ARGV[1]);
	my $i=0;
	while(my $seq = $seq_in->next_seq() ) {
		my $out = Bio::SeqIO->new(-file => ">s".$i , '-format' => 'Fasta');
		$out->write_seq($seq);
		$i++;
	}
}

sub count_fasta {
	my $seq_in = Bio::SeqIO->new( -format => 'Fasta', -file => $ARGV[1] );
	my $total=0;
	my $number=0;
	while (my $seq = $seq_in->next_seq() ) {
		my $id=$seq->id;
		my $length=$seq->length;
		print "$id\t$length\n";
		$total+=$length;	
		$number++;
	}
	print "All=$number\t$total\n";
}

sub read_codon  {
	my $CODONS= <<EOF;
Ala/A 	GCU, GCC, GCA, GCG 	
Leu/L 	UUA, UUG, CUU, CUC, CUA, CUG
Arg/R 	CGU, CGC, CGA, CGG, AGA, AGG 	
Lys/K 	AAA, AAG
Asn/N 	AAU, AAC 	
Met/M 	AUG
Asp/D 	GAU, GAC 	
Phe/F 	UUU, UUC
Cys/C 	UGU, UGC 	
Pro/P 	CCU, CCC, CCA, CCG
Gln/Q 	CAA, CAG 	
Ser/S 	UCU, UCC, UCA, UCG, AGU, AGC
Glu/E 	GAA, GAG 	
Thr/T 	ACU, ACC, ACA, ACG
Gly/G 	GGU, GGC, GGA, GGG 	
Trp/W 	UGG
His/H 	CAU, CAC 	
Tyr/Y 	UAU, UAC
Ile/I 	AUU, AUC, AUA 	
Val/V 	GUU, GUC, GUA, GUG
START/^ 	AUG 	
STOP/* 	UAA, UGA, UAG
EOF
	my %codons;
	my @lines=split "\n",$CODONS;
	foreach (@lines) {
		my @tab=split "[\t\/]", $_;
		$tab[1]=~s/\s//g;
		$tab[2]=~s/\s//g;
		$codons{$tab[1]}=[split ",",$tab[2]];
	}
	return %codons;
}
sub fisher_yates_shuffle {
	my $array = shift;
	my $i;
	for ($i = @$array; --$i; ) {
		my $j = int rand ($i+1);
		next if $i == $j;
		@$array[$i,$j] = @$array[$j,$i];
	}
}


#!/usr/bin/perl  
use strict;
use warnings;
use Statistics::Basic::Correlation;
use Statistics::LineFit;
use Math::Cephes qw(:all);

##############################################################################
#Created Dec. 10, 2005                  Last Modified  Aug. 22, 2006          #
#zhi_lu@urmc.rochester.edu													#
##############################################################################
#print out the usage informagion
&usage unless (@ARGV);
sub usage {
	print "USAGE:\n";
	print "\t filter.pl [options] (pattern) [FILE1 FILE2 ...] \n";
	print "DESCRIPTION:\n";
	print "\t read the database and filter the funtional siRNA or antisense oligos.\n";
	print "OPTIONS:\n";
	print "\t -rw  [DATABASE.db]  \n";
	print "\t\t Read the database file into arrays and output them to STDOUT.\n";
	print "\t -end 2 [DATABASE.db]  \n";
	print "\t\t Calculate the free energy difference between 5' and 3' ends, window size 2 nt.\n";
	print "\t -f <oligotype>  [DATABASE.db]/[oligo.log]  \n";
	print "\t\t <oligotype>: \n";
	print "\t\t 1		sperater scores (without hairpin Tm) of methods in Reynolds et.al. Nature Biotech. \n";
	print "\t\t 2		total score of methods in Reynolds et.al. Nature Biotech. \n";
	print "\t\t 3		siRNA features from Ladunga's paper, NAR, 2006. \n";
	print "\t\t 0		filter the functional antisense oligos. \n";



	exit;
}


######################################################################################
#main program, global variables
our (@oligonum,@targetname, @efficacy, @position, @positionstop,@targetseq,@oligoseq,@concentration);
our (@length,@lnA,@end_diff,@senseseq);
our (@dup_break,@duplex,@Tm,@breaking,@intraoligo,@interoligo);
our (@filter_score,@reads);
our $dna;


######################################################################################
#different usages
my ($reg_argv) = ($ARGV[0]=~ /(-\w*)\d/);
if ($ARGV[0] eq "-rw") {#Read data from database files
	&readdata($ARGV[1]); 
	&outputdata;
}
elsif ($ARGV[0] eq "-end") {#calculate the difference of 5' and 3' ends of oligo
#	&readdata($ARGV[2]); 
	&readtable($ARGV[2]);
	&end_diff($ARGV[1],0); #different windows size defined by $ARGV[1]
							#0 means it is a RNA not DNA oligo
}
elsif ($ARGV[0] eq "-f") {#filter the functional siRNA or antisense oligos
	
	if ($ARGV[1] == 1)	{
		&readdata ($ARGV[2]);
		&filter;
	}	#Reynold scores sperately	
	elsif ($ARGV[1]==2)	{	&readdata ($ARGV[2]); &total_filter;	}#total score of Reynold
	elsif ($ARGV[1]==3)	{	
		&readoligo($ARGV[2]);	
		&features;
	# &combine_features;
	}#siRNA features from Ladunga's paper, NAR, 2006.
	elsif ($ARGV[1]==0)	{	die ("not available");	}
		
}
else {
	print "\n\t\tUnrecognized options!!! $ARGV[0]\n\n\n";
	&usage;
}
exit;

sub filter {
	my $i=0;
	foreach (@oligoseq) {
		my @score=(0,0,0,0,0,0,0,0,0,0);
		my $seq=$_;
		#Criteria I: 30-52% G/C content 
		#+1 point
		my $numofGC= ($seq =~ tr/[GC]//);
		my	$contentofGC= $numofGC/(length $seq);
		if ($contentofGC >=0.3 && $contentofGC <=0.52)	{	$score[1]=1;	}

		#Criteria II: at least 3 'A/U' bases at positions 15-19(sense strand)( 1-5 for antisense)
		#+ 1xN points (1 point for each 'A/U', at most 5 points in total)
		for (my $j=1; $j<=5; $j++) {
			if (substr($seq,$j-1,1) =~ /[AU]/) {	$score[2]++;		}
		}

		#Criteria III: Absence of internal repeats (Tm of potential internal hairpin is< 57 degree)
		#+1 point, not using seperately
		#if ($Tm[$i] < 57) {	$score[3]=1;	}
		#Criteria IV: An 'A' base at position 19 (sense strand) (1 for antisense)
		#+1 point
		if (substr($seq,1-1,1) eq 'A')	{ 	$score[4]=1	};

		#Criteria V: An 'A' base at position 3(sense strand) (17 for antisense)
		#+1 point
		if (substr($seq,17-1,1) eq 'A')	{	$score[5]=1;	}

		#Criteria VI : A 'U' base at position 10 (sense strand) (10 for antisense)
		#+1 point
		if (substr($seq,10-1,1) eq 'U')		{	$score[6]=1;	}

		#Criteria VII: A base other than 'G' or 'C' at 19 (sense strand) (1 for antisense)
		#-1 point
		if (substr($seq,1-1,1) =~ /[GC]/)	{		$score[7]=-1;	}

		#Criteria VIII: A base other than 'G' at position 13 (sense strand) (7 for antisense)
		#-1 point
		if (substr($seq,7-1,1) eq 'G')	{	$score[8]=-1;	}

		for (my $j=1;$j<=8;$j++)	{	$score[0]+=$score[$j];	}
	
		print "$oligonum[$i]\t ";
		for (my $j=1;$j<=8;$j++)	{	
			if ($j!=3)	{				print "$score[$j]\t";		}
		}
		print "$score[0]\t $efficacy[$i]\n";
		$i++;
	}

#	print scalar @score,"\n";	
#	&correlation (\@score, \@lnA);
}

sub total_filter {
	
	my @score; #score of the filter
	my $seqfile="foo_filter.in";
	open (IN, ">$seqfile");
	foreach (@oligoseq) {
		#print sequences to a temperary file
		print IN "$_\n";
	}
	print IN "end\n"; #end to label the end of file
	close IN;
	#pipe the output of siRNAfilter to score array
	local *PIPE;
	open PIPE, "../cpp/siRNAfilter.o $seqfile |";
	while (<PIPE>) {		 
#		print $_;
		$_=~ s/\n//g;	
		push (@score,$_); 
	}			
	close PIPE;
	#remove the temperary file
  system("rm -f $seqfile");
	my $j=0;
	foreach (@score) {
		print "$oligonum[$j]\t  $targetname[$j]\t ";
		print "$score[$j]\t $efficacy[$j]\n";
		$j++;
	}
#	print scalar @score,"\n";	
	&correlation (\@score, \@lnA);
}


#siRNA features from Ladunga's paper, NAR, 2006.
sub features {
	my $i=0;
	my @score_array;
	my $tmp;
	foreach (@oligoseq) {
		my @score=(0);
		my $seq=$_;
		my @base=split('',$seq);
			
		#Criteria I: DG 1-2 
		push @score,&calc_DG(substr($seq,0,2),0);
		
		#Criteria II: 'U' base at position 1
		if ($base[1-1] eq 'U')	{ 	push @score, 1;	} 
		else 	{	push @score, 0;	}
		

		#Criteria III: 'G' base at position 1
		if ($base[1-1] eq 'G')	{ 	push @score, 1;	} 
		else 	{	push @score, 0;	}
		
		#Criteria IV: U all
		$tmp=0;
		for (my $j=1; $j<=19; $j++) {
			if ($base[$j-1] eq 'U' ) {	$tmp++;		}
		}
		push @score, $tmp; 

		#Criteria V: DG ALL
		my $DG_ALL = calc_DG(substr($seq,0,19),1);
		push @score,$DG_ALL;

		#Criteria VI : 'UU' base at position 1
		if (substr($seq,1-1,2) =~ /UU/)			{ 	push @score, 1;	} 
		else 	{	push @score, 0;	}

		#Criteria VII: 'G' all
		$tmp=0;
		for (my $j=1; $j<=19; $j++) {
			if ($base[$j-1] eq 'G' ) {	$tmp++;		}
		}
		push @score, $tmp; 

		#Criteria VIII: DDG 3-5 -- 19-21; ---- excluded
#		push @score, &calc_DG(substr($seq,3-1,3),0) - &calc_DG(substr($seq,19-1,3),0);

		#Criteria IX: DDG 1-3 -- 19-21  ----- excluded
#		push @score, &calc_DG(substr($seq,1-1,3),0) - &calc_DG(substr($seq,19-1,3),0);

		#Criteria X: 'GG' base at position 1
		if (substr($seq,1-1,2) =~ /GG/)	{ 	push @score, 1;	} 
		else 	{	push @score, 0;	}
		
		#Criteria XI: 'GC' base at position 1
		if (substr($seq,1-1,2) =~ /GC/)	{ 	push @score, 1;	} 
		else 	{	push @score, 0;	}
		
		#Criteria XII: 'UA' all
		$tmp=0;
		for (my $j=1; $j<=18; $j++) {
			if (substr($seq,$j-1,2) =~ /UA/ ) {	$tmp++;		}
		}
		push @score, $tmp; 

		#Criteria XIII: 'U' base at position 2
		if ($base[2-1] eq 'U')	{ 	push @score, 1;	} 
		else 	{	push @score, 0;	}

		#Criteria XIV:  'C' bases at position 1
		if ($base[1-1] eq 'C')	{ 	push @score, 1;	} 
		else 	{	push @score, 0;	}
		
		#Criteria XV: 'GG' all
		$tmp=0;
		for (my $j=1; $j<=18; $j++) {
			if (substr($seq,$j-1,2) =~ /GG/ ) {	$tmp++;		}
		}
		push @score, $tmp; 

		#Criteria XVI: DDG 1-5 -- 17-21  ---- excluded
#		push @score,&calc_DG(substr($seq,1-1,5),0) - &calc_DG(substr($seq,17-1,3),0);
		
		#Criteria XVII: DG 18
		push @score,&DG($base[18-1],$base[19-1]);
		 
		#Criteria XVIII :  DG 13
		push @score,&DG($base[13-1],$base[14-1]);
		
		#Criteria XIX:  DG 2
		push @score,&DG($base[2-1],$base[3-1]);
		
		#Criteria 20: 'GC' all
		$tmp=0;
		for (my $j=1; $j<=18; $j++) {
			if (substr($seq,$j-1,2) =~ /GC/ ) {	$tmp++;		}
		}
		push @score, $tmp; 

		#Criteria 21:  'CC' all
		$tmp=0;
		for (my $j=1; $j<=18; $j++) {
			if (substr($seq,$j-1,2) =~ /CC/ ) {	$tmp++;		}
		}
		push @score, $tmp; 

		#Criteria 22: 'UU' all
		$tmp=0;
		for (my $j=1; $j<=18; $j++) {
			if (substr($seq,$j-1,2) =~ /UU/ ) {	$tmp++;		}
		}
		push @score, $tmp; 

		#Criteria 23: 'CG' base at position 1
		if (substr($seq,1-1,2) =~ /CG/)		{ 	push @score, 1;	} 
		else 	{	push @score, 0;	}
		#Criteria 24: 'A' base at position 19
		if ($base[19-1] eq 'A')	{ 	push @score, 1;	} 
		else 	{	push @score, 0;	}
							
		#Criteria 25: 'CC' base at position 1
		if (substr($seq,1-1,2) =~ /CC/)		{ 	push @score, 1;	} 
		else 	{	push @score, 0;	}
		
		#Criteria 26: DH 1-2 
		push @score,&calc_DH(substr($seq,0,2),0);

		#Criteria 27: DS 1-2 
		push @score,(&calc_DH(substr($seq,0,2),0) - &calc_DG(substr($seq,0,2),0) )/0.31015;
	
		#Criteria 28: DH ALL 
		my $DH_ALL=&calc_DH(substr($seq,0,19),1);
		push @score,$DH_ALL;

		#Criteria 29: DS ALL 
		my $DS_ALL=($DH_ALL-$DG_ALL)/310.15;
		push @score,$DS_ALL*1000;

		#Criteria 30: DH/DS ALL 
		push @score,$DH_ALL/$DS_ALL;

		#Criteria 31: End of two ends this will replace criteria VIII,IX and XVI
		push @score,&DG($base[1-1],$base[2-1]) - &DG($base[18-1],$base[19-1]);
	
		#print "$oligonum[$i]\t ";
#		print $efficacy[$i]>90?1:0;
#		print "\t";
#		print "$efficacy[$i]\t ";
		for (my $j=1;$j<=28;$j++)	{	
#			print "$j:";			printf "%.1f\t",$score[$j];	
			push (@{$score_array[$j-1]},$score[$j]);

		}
#		print " \n";
		$i++;
	}

#	print scalar @score,"\n";	
	for my $row (@score_array) {		&correlation (\@$row, \@efficacy);	}
#	for my $row (@score_array) {		&correlation (\@$row, \@lnA);	}
}

#siRNA features from Ladunga's paper, NAR, 2006.
sub combine_features {
	my $i=0;
	foreach (@oligoseq) {
		my @score=(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
		my $seq=$_;
		my @base=split('',$seq);
		#Criteria 0: p3, sfold pairs possibility, not included
		
		#Criteria I: DG 
		$score[1]=&calc_DG(substr($seq,0,19));
#$score[1]=&calc_DG($seq);
		
		#Criteria II: CC
		for (my $j=1; $j<=20; $j++) {
			if ($base[$j-1] eq 'C' && $base[$j] eq 'C') {	$score[2]++;		}
		}

		#Criteria III: 'U' base at position 1
		if ($base[1-1] eq 'U')	{ 	$score[3]=1	};
		
		#Criteria IV: DH 18-19 bases
		
		$score[4]=&DH($base[18-1],$base[19-1]);

		#Criteria V: 'A' base at position 19
		if (substr($seq,19-1,1) eq 'A')	{	$score[5]=1;	}

		#Criteria VI : 'G' base at position 1
		if (substr($seq,1-1,1) eq 'G')		{	$score[6]=1;	}

		#Criteria VII: 'UU' at 18-19
		if (substr($seq,18-1,2) =~ /UU/)	{		$score[7]=1;	}

		#Criteria VIII: DH 20-21
		$score[8]=&DH($base[20-1],$base[21-1]);

		#Criteria IX: 'U' base at position 2
		if ($base[2-1] eq 'U')	{ 	$score[9]=1	};

		#Criteria X: 'A' base at position 2
		if ($base[2-1] eq 'A')	{ 	$score[10]=1	};
		
		#Criteria XI: 'AU' bases at position 6-7
		if (substr($seq,6-1,2) =~ /AU/)	{		$score[11]=1;	}
		
		#Criteria XII: 'AA' bases at position 17-18
		if (substr($seq,17-1,2) =~ /AA/)	{		$score[12]=1;	}

		#Criteria XIII:  'GG' bases at position 20-21
		if (substr($seq,20-1,2) =~ /GG/)	{		$score[13]=1;	}
		
		#Criteria XIV: 'AA' bases at position 18-19
		if (substr($seq,18-1,2) =~ /AA/)	{		$score[14]=1;	}
		
		#Criteria XV: 'AU' bases at position 9-10
		if (substr($seq,9-1,2) =~ /AU/)	{		$score[15]=1;	}
		
		#Criteria XVI: DG 3-4
		$score[16]=&DG($base[3-1],$base[4-1]);
		 
		#Criteria XVII : 'C' base at position 1
		if ($base[1-1] eq 'C')	{ 	$score[17]=1	};
		
		#Criteria XVIII: 'GG' bases at position 16-17
		if (substr($seq,16-1,2) =~ /GG/)	{		$score[18]=1;	}
		
		#Criteria XIX: 'CG' bases at position 1-2
		if (substr($seq,1-1,2) =~ /CG/)	{		$score[19]=1;	}
		
		#Criteria 20: 'AG' bases at position 20-21
		if (substr($seq,20-1,2) =~ /AG/)	{		$score[20]=1;	}
		
		#Criteria 21: 'G' base at position 14
		if ($base[14-1] eq 'G')	{ 	$score[21]=1	};

		#Criteria 22: 'UG' base at position 4-5
		if (substr($seq,4-1,2) =~ /UG/)	{		$score[22]=1;	}
		
		#Criteria 23: 'A' base at position 20
		if ($base[20-1] eq 'A')	{ 	$score[23]=1	};

		#Criteria 24: 'UG' base at position 20-21
		if (substr($seq,20-1,2) =~ /UG/)	{		$score[24]=1;	}
							
		#Criteria 25: 'CC' base at position 13-14
		if (substr($seq,13-1,2) =~ /CC/)	{		$score[25]=1;	}
								
		#Criteria 26: 'GU' base at position 5-6
		if (substr($seq,5-1,2) =~ /GU/)	{		$score[26]=1;	}
									
		#Criteria 27: 'A' base at position 1
		if ($base[1-1] eq 'A')	{ 	$score[27]=1	};

		#Criteria 28: 'CC' base at position 20-21
		if (substr($seq,20-1,2) =~ /CC/)	{		$score[28]=1;	}
		
		#Criteria 29: 'U' base at position 7
		if ($base[7-1] eq 'U')	{ 	$score[29]=1	};
	
		print "$oligonum[$i]\t ";
		for (my $j=1;$j<=29;$j++)	{	
			print "$j:$score[$j]\t";		
		}
		print " \n";
		$i++;
	}

#	print scalar @score,"\n";	
#	&correlation (\@score, \@lnA);
}


#subroutines:
###################################################################
#read the database table 
sub readdata {
	#open the database file as input information
	my($database) =@_;
#	print "Database:\t $database\n";
	open (OLIGO,"$database" ) or  &dienice("Can't open $database $!\n");
	my @filedata=<OLIGO>;
	close OLIGO; 
	#reading each site
	foreach my $line (@filedata) {
		#split the line to get each data
		my @oligo=split(" ", $line);
		push(@oligonum, $oligo[0]);
		push(@targetname,$oligo[1]);
		push(@efficacy,$oligo[2]);
		push(@position,$oligo[3]);
		push(@positionstop,$oligo[4]);
		push(@targetseq,$oligo[5]);
		push (@oligoseq,$oligo[6]);
		if ( scalar @oligo > 7 )	{	
			if($oligo[7] eq "0") 	{	$oligo[7]=100;	}
			#translate the unit from nM(?) to M
			push(@concentration,$oligo[7]/1e9);
		}
		#calculate some parameters	
		push(@length, $oligo[4]-$oligo[3]+1);
		if ($oligo[2] >= 100)	{ push (@lnA, log(1- 99.9/100));		 }
		elsif ($oligo[2] < 0)	{ push (@lnA, log(1));		 			 }
		else 			 		{ push (@lnA, log(1-$oligo[2]/100.0));	 }
	}		
}
#read a simple table
###########################
sub readtable {
	#open the database file as input information
	my($database) =@_;
#	print "Database:\t $database\n";
	open (OLIGO,"$database" ) or  &dienice("Can't open $database $!\n");
	my @filedata=<OLIGO>;
	close OLIGO; 
	#reading each site
	foreach my $line (@filedata) {
		#split the line to get each data
		my @oligo=split(" ", $line);
		push (@oligoseq,$oligo[0]);
		push (@reads,$oligo[1]);
		#calculate some parameters	
		push (@lnA, log($oligo[1]));
	}		
}


#read the database table 
sub readoligo {
	#open the database file as input information
	my($database) =@_;
#	print "Database:\t $database\n";
	open (OLIGO,"$database" ) or  &dienice("Can't open $database $!\n");
	my @filedata=<OLIGO>;
	close OLIGO; 
	#reading each site
	foreach my $line (@filedata) {
		#split the line to get each data
		my @oligo=split(" ", $line);
		push(@oligonum, $oligo[0]);
		push (@oligoseq,$oligo[1]);
		push (@efficacy,$oligo[2]);
		push (@lnA,$oligo[3]);

		
	}		
}


#output data to STDOUT
sub outputdata {
	my $i=0;
	print "oligonum  targetname  efficacy  position  positionstop  ";
	print "targetseq  oligoseq  length  lnA  ";
	if(@concentration)	{	print "concentration";	}
	print "\n";
	foreach (@oligonum) {
		print "$oligonum[$i]  $targetname[$i]  $efficacy[$i]  $position[$i]  $positionstop[$i]  ";
		print "$targetseq[$i]  $oligoseq[$i]  $length[$i]  ";
		printf "%.3f  ",$lnA[$i];
		if(@concentration) {	printf "%g ", $concentration[$i];	}
		print "\n";
		$i++;
	}
}

#write database with reversed oligo
sub writedata{
	#open the original database file
	my $database = $ARGV[1];
	&dienice("$ARGV[1] is not a text file.\n")	unless(-T $database);
	open (OLIGO,"$database" ) or  &dienice("Can't open .db file: $!\n");
	my @filedata=<OLIGO>;
	close OLIGO; 
	#open the output database file
	open(OUT,">$ARGV[2]" ) or  &dienice("Can't write file: $!\n");

	#read the sequence name from database
	my ($i,$reverse_oligo,@reverse_array);
	my (@oligo,@oligoarray);
	foreach my $line (@filedata) {
		@oligo=split(" ", $line);
		@oligoarray= split ("",$oligo[5]);
		@reverse_array="";
		foreach my $temp (@oligoarray) {
		unshift (@reverse_array,$temp);
		}
		$reverse_oligo = join("",@reverse_array);
		$reverse_oligo =~ tr/ACGTU/UGCAA/;
		$reverse_oligo =~ tr/U/T/;
		for ($i=0; $i<=4;$i++) 	{ print OUT "$oligo[$i]\t"; }
		print OUT "$reverse_oligo\t$oligo[5]\n";
	}
	close OUT;
}

#check the postion of oligo
sub checkseq {
	#open the database file
	my $database = $ARGV[1];
	open (OLIGO,"$database" ) or  &dienice("Can't open .db file: $!\n");
	my @filedata=<OLIGO>;
	close OLIGO; 
	my ($record, @oligo,@oligoid, $seqname, $sequence, $oligotest);
	my $save_input_separator=$/;
	#read the sequence name from database per line
	foreach my $line (@filedata) {
		@oligo=split(" ", $line);
		$seqname = "../sequences/".$oligo[1].".seq"; # fetch the sequence name
		#read the sequence from seq file
		open (SEQ,"$seqname" ) or  &dienice("Can't open sequence file:$seqname $!\n");
		#Set input separator and read in a record to a scalar
		$/ = "//\n"; 
		$record = <SEQ>;
		#Set the separator back to default
		$/ = $save_input_separator;  
		#fetch the sequence with regular expression
		$record =~ /(^;.*\n)*.*\n((.*\n)*)/m; 
		$sequence=$2;
		#clean the sequnce
		$sequence =~ s/\s//g;
		$sequence =~ s/\d//g;
		#check the position of the oligo
		$oligotest= substr($sequence,$oligo[3]-1,$oligo[4]-$oligo[3]+1);
		$oligotest =~ tr/acgtu/ACGTU/;
		$oligo[5] =~ tr/acgtu/ACGTU/;   
	#	$oligotest =~ s/U/T/g;		  #transcribe RNA to DNA 	
		if ($oligotest eq $oligo[5]) {
		#		print "$oligo[0] is checked to be Correct\n";
		}
		else {
			print "Target: $oligo[1] \t oligo_id: $oligo[0]\n";
			print "$oligotest\t	(target) \n";
			print "$oligo[5]\t (oligo_reverse)\n"; 
		}
		close SEQ;
	}
	print "Done! All positions are correct if no error output\n";
}

#count targets
sub counttargets{
	#open the total database file
	open (OLIGO,"$ARGV[1]" ) or  die("Can't open $ARGV[1] file: $!\n");
	my @filedata=<OLIGO>;
	close OLIGO;
	my @targets=();
	foreach my $line (@filedata) {
		my @oligo=split(" ", $line);
		#count targets
		my $findtarget=0;
		foreach (@targets) { 
			if ($_ eq $oligo[1]) {$findtarget=1;}
		}
		push (@targets, $oligo[1])	unless $findtarget;
	}
	foreach (@targets){
		my $nameoftarget=$_;
		my $numofbases;
		my $seqfile= "../sequences/$nameoftarget".".seq";
		open (SEQ,"$seqfile") or  die("Can't open seq file: $seqfile : $!\n");
			while (<SEQ>) {
				if ($_ =~ /^;\s*(\d*)\s*bases.*/){
					$numofbases=$1;	
				}
			}
		print "$nameoftarget\t\t$numofbases\n";
	}
}

#check the database from BMC bioinformatics 7:122 table 1
sub checkbmc {
	#open the total database file
	open (OLIGO,"$ARGV[1]" ) or  die("Can't open $ARGV[1] file: $!\n");
	my @filedata=<OLIGO>;
	close OLIGO;
 	open (BMC,"$ARGV[2]" ) or  die("Can't open $ARGV[2] file: $!\n");
	my @bmcdata=<BMC>;
	close BMC; 
	open (OUT,">$ARGV[3]" ) or  die("Can't open $ARGV[3] file: $!\n");
	my @targets=();
	#read the sequence name from database per line
	foreach (@bmcdata) {
		my $found=0;
		my @bmcline= split(" ", $_);
		foreach my $line (@filedata) {
			my @oligo=split(" ", $line);
			if ($bmcline[0] eq $oligo[0] && $bmcline[1] eq $oligo[6] ) {
			if ($bmcline[3] == $oligo[3] && ($bmcline[3] + $bmcline[4] -1) == $oligo[4]){
			if ( ($bmcline[5]*100-100+$oligo[2]) <1 ){
				print OUT $line;
				$found=1;
						}	
			}
			}
		}
		print "$bmcline[0]\t $bmcline[2]\n"	unless $found;
	}
	close OUT;
}

#Tranform the sequence in database from fasta to .seq file
sub transform {
#open the database file
my $database = $ARGV[1];
my $fasta = $ARGV[2];
&dienice("$ARGV[1] is not a text file.\n")	unless(-T $database);
&dienice("$ARGV[2] is not a text file.\n")	unless(-T $fasta);
open (OLIGO,"$database" ) or  &dienice("Can't open .db file: $!\n");
my @filedata=<OLIGO>;
close OLIGO; 
#open the fasta file
open (FASTA,"$fasta" ) or  &dienice("Can't open .fasta file: $!\n");
my @fastadata=<FASTA>;
close FASTA; 

#read the sequence name from database
my $linei=0;
my $i;
my @oligo;
my @oligoid;
my $exist;
my $fastaline;
my @seqname;
my $sequence;
my $findit;
my $pos;
foreach my $line (@filedata) {
	@oligo=split(" ", $line);
	#oligoid store the name of target sequnce from each line of database
	$oligoid[$linei]=$oligo[1];
	$exist=0;
	#test if the sequence is already stored in oligoid
	for ($i=0; $i<$linei;$i++) {
	if ($oligo[1] eq $oligoid[$i]) 	{ $exist=1; }
	}
	$linei++;
	#trastalte fasta to seq file if the sequence found for the first time	
	unless($exist) {
		$sequence='';
		$findit=0;
		foreach $fastaline (@fastadata) {
			if ($findit ==1) {
				if ($fastaline =~ /^>/)		{ last; }
				$sequence .= $fastaline;
			}
			if ($fastaline =~ /^>/) {
				@seqname=split(" ",$fastaline);
				substr($seqname[0],0,1)='';
				if ($seqname[0] eq $oligo[1]) {
				$findit=1;
				}
			}
		}
		if ($findit) {
			#output .seq file
			$seqname[0] .=".seq";
			$sequence =~ tr/acgtu/ACGTU/;
			$sequence =~ tr/T/U/;	#transcribe DNA to RNA nt 
			open (SEQ, ">$seqname[0]");
			print SEQ ";\n";
			print SEQ "$oligo[1]\n";
			chop $sequence;
			print SEQ "$sequence";
			print SEQ "1\n";
			close SEQ;
			print "Translated sequence from fasta: $oligo[1]\n";
		}
		else {	
		print "$oligo[1] is not founded in fasta\n";
		}
	}
}
}

#fetch the DNA sequence from GenBank files
sub fetch {
	#open the GenBack file
	my $gb = $ARGV[1];
	open (GB,"$gb" ) or  &dienice("Can't open file: $gb $!\n");
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
	$sequence =~ tr/acgtu/ACGTU/;
	$sequence =~ tr/T/U/;	#transcribe DNA to RNA nt 
	(my $seqname=$gb) =~ s/.gb/.seq/;
	#output .seq file
	open (SEQ, ">$seqname") or &dienice("Can't write file: $seqname $!\n");
	print SEQ "; ", length($sequence)," bases\n";
	$seqname=~ s/.seq//;
	print SEQ "$seqname\n";
	print SEQ "$sequence";
	print SEQ "1\n";
	close SEQ;
	print "\nDone: Transform sequence from GenBank file: $gb\n\n";
	close GB;
}

#output the errors
sub dienice {
	my($errmsg) = @_;
    print "Error\n";
    print "$errmsg\n";
    exit;
}

#initialize the free energy arrays:
sub stackenergy() {
	my ($isdna) =@_;
	my @stack;
	$stack[0][1]=0;
	$stack[0][2]=0;
	$stack[0][3]=0;
	$stack[0][4]=0;
	$stack[1][0]=0;
	$stack[2][0]=0;
	$stack[3][0]=0;
	$stack[4][0]=0;
	if ($isdna) {
	$stack[1][1]=-1.0;
	$stack[1][2]=-2.1;
	$stack[1][3]=-1.8;
	$stack[1][4]=-0.9;
	$stack[2][1]=-0.9;
	$stack[2][2]=-2.1;
	$stack[2][3]=-1.7;
	$stack[2][4]=-0.9;
	$stack[3][1]=-1.3;
	$stack[3][2]=-2.7;
	$stack[3][3]=-2.9;
	$stack[3][4]=-1.1;
	$stack[4][1]=-0.6;
	$stack[4][2]=-1.5;
	$stack[4][3]=-1.6;
	$stack[4][4]=-0.2;
	}
	else {
	$stack[1][1]=-.93;
	$stack[1][2]=-2.24;
	$stack[1][3]=-2.08;
	$stack[1][4]=-1.1;
	$stack[2][1]=-2.11;
	$stack[2][2]=-3.26;
	$stack[2][3]=-2.36;
	$stack[2][4]=-2.08;
	$stack[3][1]=-2.35;
	$stack[3][2]=-3.42;
	$stack[3][3]=-3.26;
	$stack[3][4]=-2.24;
	$stack[4][1]=-1.33;
	$stack[4][2]=-2.35;
	$stack[4][3]=-2.11;
	$stack[4][4]=-.93;
	}
	return @stack;
}

#initialize the free energy arrays:
sub stackenthalpy() {
	my ($isdna) =@_;
	my @stack;
	$stack[0][1]=0;
	$stack[0][2]=0;
	$stack[0][3]=0;
	$stack[0][4]=0;
	$stack[1][0]=0;
	$stack[2][0]=0;
	$stack[3][0]=0;
	$stack[4][0]=0;
	if ($isdna) {
	
	}
	else {
	$stack[1][1]=-6.82;
	$stack[1][2]=-11.40;
	$stack[1][3]=-10.48;
	$stack[1][4]=-9.38;
	$stack[2][1]=-10.44;
	$stack[2][2]=-13.39;
	$stack[2][3]=-10.64;
	$stack[2][4]=-10.48;
	$stack[3][1]=-12.44;
	$stack[3][2]=-14.88;
	$stack[3][3]=-13.39;
	$stack[3][4]=-11.40;
	$stack[4][1]=-7.69;
	$stack[4][2]=-12.44;
	$stack[4][3]=-10.44;
	$stack[4][4]=-6.82;
	}
	return @stack;
}




sub DG {
	my($base_1,$base_2)=@_;
	$base_1=~ tr/ACGUTIacguti/123445123445/;
	$base_2=~ tr/ACGUTIacguti/123445123445/;
	my @stack = &stackenergy(0); 
	return $stack[$base_1][$base_2];
}

sub DH {
	my($base_1,$base_2)=@_;
	$base_1=~ tr/ACGUTIacguti/123445123445/;
	$base_2=~ tr/ACGUTIacguti/123445123445/;
	my @enthalpy = &stackenthalpy(0); 
	return $enthalpy[$base_1][$base_2];
}


#energy different of two ends of oligo
sub end_diff {
	my ($windowsize,$isdna)=@_;
	#define the energy parameters
	my @end;
	if ($isdna)	{	@end=(0,0,0,0,0);	}
	else		{	@end=(0,0.45,0,0,0.45);	}
	my @stack = &stackenergy(0); 
	my $num=0;
	foreach (@oligoseq) {
		my $seq= $_;
		$seq=~ tr/ACGUTIacguti/123445123445/;
		my @numseq=split //,$seq;
		#5' end of the antisense strand
		my $DG5=0;
		my $i;
		for ($i=0;$i<$windowsize-1;$i++)	{	$DG5+=$stack[$numseq[$i]][$numseq[$i+1]];	}
		$DG5+=$end[$numseq[0]];
		$DG5-=$end[$numseq[$windowsize-1]];#better correlation with this term substrated
		#3' end of the antisense strand
		my $DG3=0;
		for ($i=length($seq)-$windowsize;$i<length($seq)-1;$i++)	{ $DG3+=$stack[$numseq[$i]][$numseq[$i+1]];}
		$DG3+=$end[$numseq[length($seq)-1]];
		$DG3-=$end[$numseq[length($seq)-$windowsize]];#better correlation with this term substrated
		#difference of two ends
		my $diff=$DG5-$DG3;
		print $oligoseq[$num],"\t",$diff,"\t",$lnA[$num],"\t",$reads[$num], "\n";
		push (@end_diff,$diff);
		$num++;
	}
	$num=0;
	print "End_diff-lnA:\t";
	print "Window size is $windowsize\n";
	#print  scalar @end_diff,"\t",scalar @lnA,"\n";
	&correlation(\@end_diff,\@lnA);
}


#DH of duplex
sub calc_DH {
	
	my ($seq,$init)=@_;
	$seq=~ tr/ACGUTIacguti/123445123445/;
	my $num=length $seq;
	my @base=split //,$seq;
	
	#define the energy parameters
	my @stack = &stackenthalpy(0); 
	my @end=(0,3.72,0,0,3.72);
	
	my $energy=0;
	#if considering duplex initiation
	if($init)	{	$energy=3.61;	}
	foreach (my $i=0;$i<=$num-2;$i++) {
		$energy+=$stack[$base[$i]][$base[$i+1]];	
	}
	if($init) {	
		$energy = $energy + $end[$base[0]] + $end[$base[$num-1]];
	}
	return $energy;
}

#DG of duplex
sub calc_DG {
	
	my ($seq,$init)=@_;
	$seq=~ tr/ACGUTIacguti/123445123445/;
	my $num=length $seq;
	my @base=split //,$seq;
	
	#define the energy parameters
	my @stack = &stackenergy(0); 
	my @end=(0,0.45,0,0,0.45);
	
	my $energy=0;
	#if considering duplex initiation
	if($init)	{	$energy=4.09;	}
	foreach (my $i=0;$i<=$num-2;$i++) {
		$energy+=$stack[$base[$i]][$base[$i+1]];	
	}
	if($init) {	
		$energy = $energy + $end[$base[0]] + $end[$base[$num-1]];
	}
	return $energy;
}


#statistic analysis:finding correlations
sub correlation {
	my ($X, $Y)=@_;
	my $co = new Statistics::Basic::Correlation ([@$X],[@$Y]);
	my $R= $co->query;
	my $RSquared = $R**2;
	my $linefit = Statistics::LineFit->new();
 	$linefit->setData (\@$X, \@$Y) or die "Invalid data";
 	(my $tStatIntercept, my $tStatSlope) = $linefit->tStatistics();
	my $num=@$Y;
	my $pvalue= stdtr($num-2,abs($tStatSlope)); #funciont in math module(cephes)
	my $significance=(1-$pvalue)*2; #two tail distribution
	printf "R:\t %.4f;\t R-Squared\t %.4f;\t",$R, $RSquared;
	printf  "Significance:\t %.3g\n",$significance;
}



#read the energies from output files of oligowalk, each output has multiple sites
sub mreadenergy {
	$dna=0;
	my($foldsize, $mode, $database) =@_;
	print "\nDatabase:\t $database\t Folding size:\t $foldsize \t Mode:\t $mode\n";
	#open the database file as input information
	open (OLIGO,"$database" ) or  &dienice("Can't open $database $!\n");
	my @filedata=<OLIGO>;
	close OLIGO; 
	#reading each site
	foreach my $line (@filedata) {
		my @oligo=split(" ", $line);
		push(@length, $oligo[4]-$oligo[3]+1);
		my $foldsize_tr= $foldsize-$oligo[4]+$oligo[3]-1;
		push(@oligonum, $oligo[0]);
		push(@targetname,$oligo[1]);
		push (@position,$oligo[3]);
		push (@efficacy,$oligo[2]);
		if ($oligo[2] >= 100)	{ push (@lnA, log(1- 99.9/100));		 }
		elsif ($oligo[2] < 0)	{ push (@lnA, log(1));		 			 }
		else 			 		{ push (@lnA, log(1-$oligo[2]/100.0));	 }
		push (@senseseq,$oligo[5]);
		push (@oligoseq,$oligo[6]);
		#reading the energies from file generated by oligowalk, named by target name	
		my $output="/home/john/oligowalk/output_siRNA_all/".$oligo[1]."_".$foldsize."_".$mode.".out";
		unless (-e $output) {
		$output="/home/john/oligowalk/output_siRNA_all/".$oligo[1]."_0_".$mode.".out";
		}
		unless (-e $output) {
		$output="/home/john/oligowalk/output_siRNA_all/".$oligo[1]."_4000_".$mode.".out";
		}
		open (OUT,"$output" ) or  &dienice("Can't open $output $!\n");
		my @outdata=<OUT>;
		close OUT;
		my $linei=0;#the flag for reading energy table
		my $outposition;
		my $outoligoseq;
		my $foundit = 0;
		foreach my $outline(@outdata) {
			if($outline=~/<\/table>/) {#check if line hit the end of energy table
				last;#break out of the loop
			} elsif ($linei <4 && $linei >=1) {#read from the 4th line of the table
				$linei++;
			} elsif($linei >=4)	{#read the energy talbe to find the oligo site
				$outline=~ s/<tr>//g;
				$outline=~ s/<td>//g;
				$outline=~ s/<\/tr>//g;
				$outline=~ s/<\/td>//g;
				my @energies=split(" ",$outline);
				$outposition= $energies[0];
				if($outposition == $oligo[3]) {#if the site in database was found in oligowalk output
					$outoligoseq= $energies[1];
					$foundit=1;
					push(@duplex,$energies[3]);
					push(@Tm, $energies[4]);
					push(@breaking,$energies[5]);
					push(@intraoligo, $energies[6]);
					push(@interoligo, $energies[7]);
					push(@end_diff, $energies[8]);
					
					last;#break out of the loop if the oligo site was found
				}
			} elsif($outline=~/<table>/) {#if line begins a energy table
				$linei=1;# set $linei flag and prepare to read next 
			}
		}
		#check if oligo site in database was found in oligowalk output file
		unless($foundit)	{print "error: pos. $oligo[3] in target $oligo[1] is not found\n"; exit;}
		#check if the oligo site's sequence is the same between database and oligowalk output
		if ($outoligoseq ne $oligo[6])	{print "error: Uncompatible oligo $oligo[6] in target $oligo[1]\n"; exit;}
	}		

}



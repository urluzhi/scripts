#!/usr/bin/perl 
#===============================================================================
#
#         FILE:  siRNAWalk.pl
#
#        USAGE:  ./siRNAWalk.pl 
#
#  DESCRIPTION:  
#
#      OPTIONS:  ---
# REQUIREMENTS:  ---
#         BUGS:  ---
#        NOTES:  ---
#       AUTHOR:  Zhi John Lu <urluzhi@gmail.com>
#      COMPANY:  University of Rochester, Medical Center
#      VERSION:  1.0
#      CREATED:  10/15/2007 15:33:11 EDT
#     REVISION:  ---
#===============================================================================

use strict;
use warnings;

#set global variables
#############################
my (%Form, $mainstring);
my (@position, @oligoseq,@length, @filter_score,@predict_prob);
my (@overall,@dup_break,@duplex,@Tm,@breaking,@intraoligo,@interoligo,@end_diff);
my $svm_scaler = "/var/www/cgi-bin/server_exe/oligowalk/svm/svm_scaler";
my $svm_model = "/var/www/cgi-bin/server_exe/oligowalk/svm/svm_model";

#main flow
####################################
#Read input options
&read_options;
#parse energy talbe
&parse_energy ;
#output_energy ($Form{'energy'});

#run svm to identify efficient siRNA
my $svm_output=$Form{'siRNA'}.".predict"; my $svm_input= $Form{'siRNA'}.".svm"; 
&calc_features;
&run_svm;
&out_siRNA ;
#clean files
unlink $svm_output,$svm_input;

exit;


#subroutines:
##########################################################
sub read_options {
#	$Form{'energy'} = $ARGV[0];
	$Form{'siRNA'} = $ARGV[0];
}

sub parse_energy {
		my @outdata=<STDIN>;
		my $linei=0;#the flag for reading energy table
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
				push(@position,$energies[0]);
				push(@oligoseq,$energies[1]);
				push(@overall, $energies[2]);
				push(@duplex,$energies[3]);
				push(@Tm, $energies[4]);
				push(@breaking,$energies[5]);
				push(@intraoligo, $energies[6]);
				push(@interoligo, $energies[7]);
				push(@end_diff, $energies[8]);
				push(@filter_score, $energies[9]);
			}elsif($outline=~/<table>/) {#if line begins a energy table
				$linei=1;# set $linei flag and prepare to read next
			}	
		}
}

#write the energies read from oligowalk output files to the table
sub output_energy {
	my ($output)=@_;
	my $body;
	$body= <<EOF;
position:	5' end of the target site
overall:	overall free energe change of the equilibrium 
duplex:		free energy change of the sense-antisense hybrid duplex
breaking:	the free energy cost to open the binding site's base pairs
Tm:		melting temperature of sensen-antisense duplex
intraoligo:	free energy change of intra-oligomer structure
interoligo:	free energy change of inter-oligomer (antisense dimer) structure
end_diff:	the difference of the free enegy change of the two ends of siRNA duplex (5' end - 3' endwindow size: 2)
sequence:	the sequence of antisense siRNA (5'-3')
---------------------------------------------------------------------------------

EOF
	$body .= "position\t overall\t duplex\t breaking\t Tm\t intraoligo\t interoligo\t end_diff\t pre_filter_score\t sequence\n";
	my $num=@position;
	for(my $j=0; $j<$num; $j++) {
		$body .= "$position[$j]\t ";
		$body .= "$overall[$j]\t";
		$body .= "$duplex[$j]\t $breaking[$j]\t $Tm[$j]\t";
		$body .= "$intraoligo[$j]\t $interoligo[$j]\t $end_diff[$j]\t $filter_score[$j]\t $oligoseq[$j]\n";
	}
	open (OUT,">$output") or die ("cannot write $output $!");
	print OUT  $body; 
	close OUT;
}

sub run_svm {
	my $svm_input_scale = $svm_input.".scale";
	#svm scaling and prediction
	system ("svm-scale -r $svm_scaler $svm_input > $svm_input_scale");
	system ("svm-predict -b 1 $svm_input_scale $svm_model $svm_output");
	unlink $svm_input_scale;
}

sub out_siRNA {
	#read probablity predicted from svm
	read_prob ($svm_output);
	my @lines ;
	my $num=@position;
	for(my $j=0; $j<$num; $j++) {
		push @lines, "$position[$j]\t $predict_prob[$j]\t $oligoseq[$j]\n";
	}
	#sort the lines by prodited probability from svm
	my @sorted_lines = reverse sort  {
		my @a = split " ",$a;
		my @b = split " ",$b ;
		($a[1]) <=> ($b[1])
	}	@lines;
	#output to file
	my $output= $Form{'siRNA'};
	open (OUT, ">$output") or die ("Cannot write to $output $!");
	print OUT "<table>\n";
	print OUT "<tr><td colspan=\"3\"><span class=\"subheader\">siRNA Candidates</span> </td></tr>\n";
	print OUT "<tr><th>Position on target</th>\t <th>Probability of being efficient siRNA</th>\t <th>siRNA Sequence(5\'->3\')</th></tr>\n";
	my $i=0;
	foreach (@sorted_lines) {
		$i++;
		my @tabs=split " ",$_;
		if ($i % 2 ==1) {
			print OUT "<tr class=\"row-a\"><td>$tabs[0]</td>  <td>$tabs[1]</td>  <td>$tabs[2]</td></tr>\n";
		}	
		else {
			print OUT "<tr class=\"row-b\"><td>$tabs[0]</td>  <td>$tabs[1]</td>  <td>$tabs[2]</td></tr>\n";
		}
	}
	print OUT "</table>\n";
	close OUT;	
}
#read probability predicted from svm
sub read_prob {
	my ($predict)= @_;
	open (Prob, "$predict") or die ("Cannot open $predict");
	my @predictdata=<Prob>;
	close Prob;
	foreach (@predictdata) {
		my @predict_col=split(" ", $_);
		push(@predict_prob, $predict_col[2]);
	}
	my $label=shift @predict_prob;	
	my $j=0;
	foreach (@predict_prob) {
		my $efficacy=0;
		unless ($label)	{	$predict_prob[$j] =1- $predict_prob[$j];}
		$j++;
	}
}

#siRNA features from Ladunga's paper, NAR, 2006.
sub calc_features {
	my $i=0;
	my $tmp;
	#make svm input file
	open (SVM, ">$svm_input") or die ("Cannot write to $svm_input $!");
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
		
		#Criteria 26: DH 1-2 ; 
		push @score,&calc_DH(substr($seq,0,2),0);

		#Criteria 27: DS 1-2 ; not used for svm
#		push @score,(&calc_DH(substr($seq,0,2),0) - &calc_DG(substr($seq,0,2),0) )/0.31015;
	
		#Criteria 28: DH ALL 
		my $DH_ALL=&calc_DH(substr($seq,0,19),1);
		push @score,$DH_ALL;

		#Criteria 29: DS ALL ; not used for svm
#		my $DS_ALL=($DH_ALL-$DG_ALL)/310.15;
#		push @score,$DS_ALL*1000;

		#Criteria 30: DH/DS ALL not used for svm
#		push @score,$DH_ALL/$DS_ALL;

		#Criteria 31: End of two ends this will replace criteria VIII,IX and XVI; already has it in prevoius program
#		push @score,&DG($base[1-1],$base[2-1]) - &DG($base[18-1],$base[19-1]);
	
		print SVM "1 ";
		for (my $j=1;$j<=24;$j++){	print SVM  "$j:$score[$j] ";	}
		print SVM "25:$intraoligo[$i] 26:$interoligo[$i] 27:$end_diff[$i] 28:$breaking[$i]\n";
		$i++;
	}
	close SVM;
}

#basic free energy and enthalp parameters
######################################################################
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
######################################################


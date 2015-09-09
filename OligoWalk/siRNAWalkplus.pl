#!/usr/bin/perl 
use strict;
use warnings;


while(<STDIN>) {
	my $seq;
	$seq=$_;
	chomp $seq;	
	print $seq;
	calc_features($seq);

}



##########################################################
#siRNA features 
sub calc_features {
	my ($seq)=@_;
	my $tmp; my @score;
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

	
		#Criteria 31: End of two ends this will replace criteria VIII,IX and XVI; may have been calculated in prevoius program
		$tmp=&DG($base[1-1],$base[2-1]) - &DG($base[18-1],$base[19-1]);
		push @score,$tmp;	

		foreach (@score) {print "\t",$_};
		print "\n";
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

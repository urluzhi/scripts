#!/usr/bin/perl 
use strict;
use warnings;
use GENOME;

#############################################################################
#print out the usage inforagion
my $usage = <<EOF;
USAGE:
	upSignal_from_bin.pl [loc] [bins] (header) 

EOF
die $usage if (!@ARGV||$ARGV[0] eq '-h'||$ARGV[0] eq '--help');
#############################################################################
#main
	if($ARGV[2]) {	print $ARGV[2],"\n";}	
	
	#Read Bin values
	open (BIN,"grep ^chr1_ $ARGV[1]|") or die ("Error:	Cannot open $ARGV[1]\n");
	my @info;
	push @info,0;
	while (<BIN>) {
		my @t=split;
		push @info, $t[1];
	}
	close BIN;
	my $total=scalar @info;
	#Map to GENEN 
	open (GENE,"grep -P \"chr1\t\" $ARGV[0]|") or die ("Error:	Cannot open $ARGV[0]\n");
	while (<GENE>) {
		my @bins;
		my @t=split;
		my $i=int($t[1]/100);
		my $j=int($t[2]/100);
		my $max=0;
		if($t[3] eq '+' && $i+1-50 >0) {
			for (my $k=1;$k<=50;$k++) {
				my $tmp=$info[$i+1-$k];
				if( $max< $tmp) {	$max=$tmp;	}	
			}
		}
		elsif($t[3] eq '-' && $j+1+50 <$total) {
			for (my $k=1;$k<=50;$k++) {
				my $tmp=$info[$j+1+$k];
				if( $max< $tmp) {	$max=$tmp;	}	
			}
		}
		printf "%.2f\n",$max;
	}
	close GENE;
	exit;

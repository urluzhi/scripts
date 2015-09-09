#!/usr/bin/perl 
use strict;
use warnings;
use GENOME;

#############################################################################
#print out the usage inforagion
my $usage = <<EOF;
USAGE:
	locSignal_from_bin.pl [loc] [bins] (header) 

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
		my $sig= ($info[$i+1]+$info[$j+1])/2;	
		printf "%.2f\n",$sig;
	}
	close GENE;
	exit;

#!/usr/bin/perl 
use strict;
use warnings;
use GENOME;

our $ref='/home9-2/gersteinlab/common/encode/chromatin/hs_bins_ref.test.bedgraph';
our $overlap_prog='/home1/zl222/bin/Overlapper/BinOverlapper';
#############################################################################
#print out the usage inforagion
my $usage = <<EOF;
USAGE:
	genesignal_from_bin.pl [BEDGRAPH]  

EOF
die $usage if ($#ARGV+1!=2||$ARGV[0] eq '-h'||$ARGV[0] eq '--help');
#############################################################################
#main
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
	#Map to GENEN tss and tes
	open (GENE,"grep -P \"chr1\t\" $ARGV[0]|") or die ("Error:	Cannot open $ARGV[0]\n");
	while (<GENE>) {
		my @bins;
		my @t=split;
		my $i=int($t[1]/100);
		if($i-39<=0) {next;}
		for (my $j=-39;$j<=40;$j++) {	push @bins, $info[$i+$j]; 	}
		$i=int($t[2]/100);
		if($i+40>=$total) {next;}
		for (my $j=-39;$j<=40;$j++) {	push @bins, $info[$i+$j]; 	}
		
		(my $line=$_)=~ s/\n//; $line=~s/\t/_/g;	
		print $line,"\t";
		if($t[3] eq '+') {	
			for (my $j=0;$j<=158;$j++) {	print  $bins[$j],"\t"; 	}
			print $bins[159],"\n";
		}
		elsif($t[3] eq '-') {	
			for (my $j=159;$j>=1;$j--) {	print  $bins[$j],"\t"; 	}
			print $bins[0],"\n";
		}
	}
	close GENE;
	exit;

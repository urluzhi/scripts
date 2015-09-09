#!/usr/bin/perl 
use strict;
use warnings;
use GENOME;

our $ref='/home9-2/gersteinlab/common/encode/chromatin/hs_bins_ref.bedgraph';
our $overlap_prog='/home1/zl222/bin/Overlapper/BinOverlapper';
#############################################################################
#print out the usage inforagion
my $usage = <<EOF;
USAGE:
	bedgraph2bin.pl [BEDGRAPH]  

EOF
die $usage if ($#ARGV+1!=1||$ARGV[0] eq '-h'||$ARGV[0] eq '--help');
#############################################################################
#main 
	open (BED,"$ARGV[0]") or die ("Error:	Cannot open $ARGV[0]\n");
	my %info;
	my $line_i=0;
	while (<BED>) {
		$line_i++;
		if($line_i==1) {next};
		my @t=split '\t',$_;
		my $chr=$t[0];
		$chr=~s/chr//;
		my $id=$chr."_".$t[1]."_".$t[2];
		$info{$id}= $t[3];
	}
	close BED;
	#map the information
	open(OVERLAP,"$overlap_prog $ARGV[0] $ref|" );
	while (<OVERLAP>) {
		#fetch the information
		my @t=split;
		my $chr=$t[0];
		$chr=~s/chr//;
		my $sig=0;
		if ($t[3])	{
			my @oligo=split ",",$t[3];
			foreach (@oligo) {
				$_=~/^(\d+)-(\d+)$/;
				my $id=$chr."_".$1."_".$2;
				$sig+=$info{$id}*(calc_coverage($t[1],$t[2],$1,$2)-1);  
			}
		}
		if($sig/100 == 0) {
			printf "chr%s_%d\t%d\n",$chr,$t[2]/100,$sig/100;
		}
		else {
			printf "chr%s_%d\t%.2f\n",$chr,$t[2]/100,$sig/100;
		}
	}
	close OVERLAP;

#############################################################
#caculate overlaped size
##############################################################
sub calc_coverage { # could be num_of_distance bases away
	my ($s1,$e1,$s2,$e2) = @_;	
	my @tmp= ($e1-$s2,$e2-$s1,$e1-$s1,$e2-$s2);
	return minNum(\@tmp)+1; 
}
sub minNum {
	my ($a)=@_;
	my $min=$$a[0];
	foreach (@$a) {
		if($min>$_) {$min=$_;}
	}
	return $min;
}




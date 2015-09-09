#!/usr/bin/perl 
use strict;	use warnings;	# CREATED: Zhi John Lu 08/25/09 12:14:57 EDT
#===============================================================================

my %info;
my ($start,$end,$chr,$strand);
open (IN,"sort -k6r,6 -k1,1 -k4,4n -k5,5n -k2,2n -k3,3n|"); 
my $linei=1;
while (<IN>) {
	if($linei==1){$linei++;next;}
	my @t=split /[\t\n]/, $_;
	my $gene="$t[0]:$t[3]:$t[4]:$t[5]:$t[6]";
	$info{$gene}.= "$t[1]_$t[2],";
}
close(IN);
foreach (keys %info) {
	my $gene=$_;
	my @gene_coor=split /:/,$gene;
	my @info_coor=split /,/,$info{$gene};
	my @exon_start; my @exon_end;
	foreach(@info_coor) {
		my @exon_coor=split /_/,$_;
		push @exon_start,$exon_coor[0];
		push @exon_end, $exon_coor[1];
	}	
	$chr=$gene_coor[0];
	if ($gene_coor[3] eq "1") {$strand="+";}
   	elsif($gene_coor[3] eq "-1") {$strand="-";}	
	
	$start=$gene_coor[1];
	$end=$exon_start[0]-1;
	if($end-$start>10) {print "chr$chr\tENSEMBL_GENE\tintron\t$start\t$end\t.\t$strand\t.\t$gene\n";}
	my $num=scalar @exon_start;
	for (my $i=0;$i<$num-1;$i++){
		$start=$exon_end[$i]+1;
		$end=$exon_start[$i+1]-1;
		if($end-$start>10) {print "chr$chr\tENSEMBL_GENE\tintron\t$start\t$end\t.\t$strand\t.\t$gene\n";}
	}
	$start=$exon_end[$num-1]+1;
	$end=$gene_coor[2];
	if($end-$start>10) {print "chr$chr\tENSEMBL_GENE\tintron\t$start\t$end\t.\t$strand\t.\t$gene\n";}
}

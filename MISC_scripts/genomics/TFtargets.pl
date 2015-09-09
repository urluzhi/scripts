#!/usr/bin/perl 
use strict;	use warnings;	# CREATED: Zhi John Lu 04/07/10 15:13:33 EDT
#===============================================================================
my $known_genes="/data/zl222/Projects_silent/gff/worm/gff4TF/coding-noncoding.gff3";
my $dict="/data/zl222/Projects_silent/gff/worm/gff4TF/elegans_ws170_exon_structures.txt";
my $upstream1=500;
my $upstream2=2000;
my $downstream=300;
my $usage = <<EOF;
USAGE:
	./TFtargets.pl [PEAKS.bed] [TRANSCRIPTS.gff] (add_known_genes) 
EOF
die $usage unless (@ARGV); 
#get gene pulic id for coding transcripts 
my %info; 
open(DICT,"$dict" ) or die("Can't read file: $dict $!\n");
while (<DICT>) {
  my @t=split "\t",$_;
  if ($t[1] && $t[3]) {
    $info{$t[3]}    = $t[1];
  }
}
close DICT;
#search each chomosome
##############################################################
my @Chromosomes=('I','II','III','IV','V','X');
print "Peak_position\tPeak_raw_signal\tPeak_input_raw_signal\tPeak_q-value\tTarget_confidence\tTargeted_Gene_ID\tTargeted_Transcript_ID\tDistance_to_Transcript_Start_Site\tTargeted_Transcript_Position\n";
foreach my $chromosome (@Chromosomes){
	# read peak
	open (PEAKS, "cat $ARGV[0] |sed 's/^chr//g' |cut -f 1-5,8,9|grep -P  \"^$chromosome\\t\" |sort -k2,2n -k3,3n| " );
	while (<PEAKS>) {
		my @p_tab=split;
		my $p_center=int (($p_tab[1]+$p_tab[2])/2);
		my (@targets_id,@targets_dist,@targets_pos);
		my (@targets_id1,@targets_dist1,@targets_pos1);
		my (@targets_id2,@targets_dist2,@targets_pos2);
		my ($found1,$found2,$inside_gene);
		# read transcripts 
		my ($chr,$start,$end,$strand,$id,$dist);
		if($ARGV[2] && $ARGV[2] eq 'add_known_genes') {
			open(TRANS, "cat $known_genes $ARGV[1]  |sed 's/^chr//g'|grep -P \"^$chromosome\\t\" | sort  -k4,4n -k5,5n|" );
		}
		else {
			open(TRANS, "cat $ARGV[1] |sed 's/^chr//g'|grep -P \"^$chromosome\\t\" | sort  -k4,4n -k5,5n|" );
		}
		while(<TRANS>) {
			my @tab=split;
			$chr=$tab[0];
			$start=$tab[3];
			$end=$tab[4];
			$strand=$tab[6];
			$id=$tab[8];
			#calc. the distance between TSS and peak center
			if($strand eq "+") {$dist=$start-$p_center;}
			elsif($strand eq "-") {$dist=$p_center-$end;}
			elsif($strand eq ".") {
				if($p_center<=$start) { $dist=$start-$p_center;}
				if($p_center>=$end) { $dist=$p_center-$end;}
				else { 
					$dist=$p_center-$start;
					if($dist>$end-$p_center) {	$dist=$end-$p_center;}
					$dist=-$dist;
				}
			}
			#pre-check if inside a transcript
			#if($p_center>=$start && $p_center<=$end) {	$inside_gene=1;	}
			#search if the peak is up500nt or down300nt(have to be inside a gene) of TSS
			if($dist <=$upstream1 && $dist >= -$downstream && $dist >= $start-$end ) {     
				push @targets_id1, $id;
				push @targets_dist1, $dist;
				push @targets_pos1, $chr."_".$start."_".$end."_".$strand;
				$found1=1;
			}
			#go ahead to search up501-2000nt upstream if -300-500nt not found yet
			elsif(!$found1 && $dist>$upstream1 && $dist <= $upstream2 ) {     
					push @targets_id2, $id;
					push @targets_dist2, $dist;
					push @targets_pos2, $chr."_".$start."_".$end."_".$strand;
					$found2=1;
			}
		}
		close(TRANS);
		#select the targets
		if($found1) {		# only select if the peak is up500nt or down300nt(have to be inside a gene) of TSS
			@targets_id=@targets_id1;
			@targets_dist=@targets_dist1;
			@targets_pos=@targets_pos1;
		}
		elsif($found2) {	# select 501-2000nt upstream if -300-500 not found
			@targets_id=@targets_id2;
			@targets_dist=@targets_dist2;
			@targets_pos=@targets_pos2;
		}
		#map transcript id to gene id
		my @targets_gene_id; my @out;	
		if (@targets_id) {  
			foreach(@targets_id) {	
				if($info{$_}) {
					push @targets_gene_id, $info{$_};
					$_.="(".$info{$_}.")";	
				}		
				else {  	push @targets_gene_id, $_;   	}		
			}
			my  %saw;	   
			@out = grep(!$saw{$_}++, @targets_gene_id);
		}
		#output peaks and targets
		if($p_tab[6] && $p_tab[6] =~ /_/) { print "$p_tab[6]:";}
		if($p_tab[5]) {	print "$p_tab[0]_$p_tab[1]_$p_tab[2]\t$p_tab[3]\t$p_tab[4]\t$p_tab[5]\t";}
		else {	print "$p_tab[0]_$p_tab[1]_$p_tab[2]\t.\t.\t.\t";}
		if($found1) {print "high\t";}
		elsif($found2) {print "low\t";}
		$"=",";
		print "@out\t@targets_id\t@targets_dist\t@targets_pos\n";
	}
	close(PEAKS);
}
exit;



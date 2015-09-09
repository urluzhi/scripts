#!/usr/bin/perl 
use strict;	use warnings;	# CREATED: Zhi John Lu 03/02/2009 11:50:52 EST
#===============================================================================
my $usage = <<EOF;
OPTIONS:
	-fetch [MAP] [GFF]: run main  
	-t2p <transcript id>  > [public id]
	-t2w <transcript id>  > [wormbase gene id]
	-w2p <wormbase id>  > [public id]
	-s2p <gene sequence name>  > [public id]
	-p2w <public id>  > [wormbase id]
	-p2t <public id>  > [transcript id]
	-p2s <public gene id>  > [sequence name of gene]
	-t2e <transcript id>  > [ENSEMBL id]
	-p2e <public id>  > [ENSEMBL id]
	-ratio <transcript id> [dcpm] > [ratio]
	-net <netid> [DICT] > tid.net
	-orth <ORTHOLOG> [id1] > [id1-id2]
	-ama1: check the correlation between GFP an polII
	-bins size <BED>: split the positive his into bins of \$size nt
	-seqsgr [SGR] <BINS>: average CHIP seq sgr ratio
	-seqwigv [SGR] <BINS>: average CHIP seq sgr ratio from variable wig file
EOF
#options:
my $mydict="/home1/zl222/Projects/gff/worm/gff4TF/elegans_ws170_exon_structures.txt";
my %option_hash = (-fetch=>\&fetch,
   					-t2p=>\&t2p,
   					-t2w=>\&t2w,
   					-w2p=>\&w2p,
   					-s2p=>\&s2p,
   					-p2w=>\&p2w,
   					-p2t=>\&p2t,
   					-p2s=>\&p2s,
   					-net=>\&net,
   					-ratio=>\&ratio,
   					-t2e=>\&t2e,
   					-p2e=>\&p2e,
   					-orth=>\&ortholog,
				   	-ama1=>\&ama1,
					-bins=>\&bins,
					-seqwigv=>\&seqwigv,
					-seqsgr=>\&seqsgr);
##############################################################
die $usage unless (@ARGV); my $option = ($ARGV[0]);
die("Unknown option: $option\n$usage\n") unless (defined($option_hash{$option}));
&{$option_hash{$option}};
exit;

#subroutines
##############################################################
#ama1 
##############################################################
sub ama1{
	open(MAP,"$ARGV[1]" ) or die("Can't read file: $ARGV[1] $!\n");
	open(INCLU,"$ARGV[2]" ) or die("Can't read file: $ARGV[2] $!\n");
	my %info; my %bothinfo;
	#read the information
	while (<INCLU>) {
		my @t=split;
		my $start=$t[1];
		my $end=$t[2];
		my $id1=$t[0]."_".$start."_".$end;
		$info{$id1}= $t[5];
	}
	#map the information
	while (<MAP>) {
		#fetch the information
		my @t=split;	
		my $id=$t[0]."_".$t[1]."_".$t[2];
		if($info{$id}) {
			$bothinfo{$id}=$t[5]."\t".$info{$id};		
		}
	}
	foreach (keys %bothinfo) {
		print $_,"\t",$bothinfo{$_},"\n";
	}
	close MAP;
	close INCLU;
}

##############################################################
#Add fetch 
##############################################################
sub fetch{
	open(MAP,"$ARGV[1]" ) or die("Can't read file: $ARGV[1] $!\n");
	open(INCLU,"$ARGV[2]" ) or die("Can't read file: $ARGV[2] $!\n");
	my %info; my %transcripts;
	#read the information
	while (<INCLU>) {
		my @t=split;
		$_=~/Transcript=\"(.*?)\"/;
		my $start=$t[3];
		my $end=$t[4];
		my $id1=$t[0]."_".$start."_".$end;
		$info{$id1}= $1;
	}
	#map the information
	while (<MAP>) {
		#fetch the information
		my @t=split;	
		my $id=$t[0]."_".$t[1]."_".$t[2];
		$transcripts{$info{$id}}=1;		
	}
	foreach (keys %transcripts) {
		print $_,"\n";
	}
	close MAP;
	close INCLU;
}

##############################################################
#get gene ensembl id 
##############################################################
sub t2e {
	open(DICT,"$mydict" ) or die("Can't read file: $ARGV[1] $!\n");
	my %info; my %transcripts;
	my %info2;
	#read the information
	while (<DICT>) {
		my @t=split "\t",$_;
		if ($t[0] && $t[5]) {
			$info{$t[5]}	= $t[0];
		}
		
	}
	#map the information
	while (<STDIN>) {
		#fetch the information
		my @t=split;	
		if ($info{$t[0]} ) {
			print $info{$t[0]},"\n";
		}
		else {
			die ("No gene Wormbase id was found for $t[0]");
		}
	}
	close DICT;
}
##############################################################
#get gene ensembl id 
##############################################################
sub p2e {
	open(DICT,"$mydict" ) or die("Can't read file: $ARGV[1] $!\n");
	my %info; my %transcripts;
	my %info2;
	#read the information
	while (<DICT>) {
		my @t=split "\t",$_;
		if ($t[1] && $t[0]) {
			$info{$t[1]}	= $t[0];
		}
		
	}
	#map the information
	while (<STDIN>) {
		#fetch the information
		my @t=split;	
		if ($info{$t[0]} ) {
			print $info{$t[0]},"\n";
		}
		elsif ($t[0] eq 'F44E2.6') {
			print "F44E2.6\n";
		}
		else {
			die ("No gene public name was found for $t[0]");
		}
	}
	close DICT;
}
sub ratio{
	open(DICT,"$mydict" ) or die("Can't read file: $ARGV[1] $!\n");
	my %info; my %transcripts;
	my %info2;
	#read the information
	while (<DICT>) {
		my @t=split "\t",$_;
		$info{$t[0]}	= $t[1]."\t".$t[2];
		
	}
	#map the information
	while (<STDIN>) {
		#fetch the information
		my @t=split;	
		if ($info{$t[0]} ) {
			print $t[0],"\t",$info{$t[0]};
		}
		else {
			die ("No gene public name was found for $t[0]");
		}
	}
	close DICT;
}
##############################################################
#get gene pulic id 
##############################################################
sub t2p{
	open(DICT,"$mydict" ) or die("Can't read file: $ARGV[1] $!\n");
	my %info; my %transcripts;
	my %info2;
	#read the information
	while (<DICT>) {
		my @t=split "\t",$_;
		if ($t[1] && $t[5]) {
			$info{$t[5]}	= $t[1];
		}
		
	}
	#map the information
	while (<STDIN>) {
		#fetch the information
		my @t=split;	
		if ($info{$t[0]} ) {
			print $info{$t[0]},"\n";
		}
		elsif ($t[0] eq 'F44E2.6') {
			print "F44E2.6\n";
		}
		else {
			die ("No gene public name was found for $t[0]");
		}
	}
	close DICT;
}

##############################################################
#get gene pulic id 
##############################################################
sub t2w{
	open(DICT,"$mydict" ) or die("Can't read file: $ARGV[1] $!\n");
	my %info; my %transcripts;
	my %info2;
	#read the information
	while (<DICT>) {
		my @t=split "\t",$_;
		if ($t[0] && $t[5]) {
			$info{$t[5]}	= $t[0];
		}
		
	}
	#map the information
	while (<STDIN>) {
		#fetch the information
		my @t=split;	
		if ($info{$t[0]} ) {
			print $info{$t[0]},"\n";
		}
		else {
			die ("No gene public name was found for $t[0]");
		}
	}
	close DICT;
}
##############################################################
#get gene pulic id 
##############################################################
sub w2p{
	open(DICT,"$mydict" ) or die("Can't read file: $ARGV[1] $!\n");
	my %info; my %transcripts;
	my %info2;
	#read the information
	while (<DICT>) {
		my @t=split "\t",$_;
		if ($t[0] && $t[1]) {
			$info{$t[0]}	= $t[1];
		}
		
	}
	#map the information
	while (<STDIN>) {
		#fetch the information
		my @t=split;	
		if ($info{$t[0]} ) {
			print $info{$t[0]},"\n";
		}
		else {
			print "$t[0]\n"
			#die ("No gene public name was found for $t[0]");
		}
	}
	close DICT;
}
##############################################################
#get gene pulic id 
##############################################################
sub s2p{
	open(DICT,"$mydict" ) or die("Can't read file: $ARGV[1] $!\n");
	my %info; my %transcripts;
	my %info2;
	#read the information
	while (<DICT>) {
		my @t=split "\t",$_;
		if ($t[2] && $t[1]) {
			$info{$t[2]}	= $t[1];
		}
		
	}
	#map the information
	while (<STDIN>) {
		#fetch the information
		my @t=split;
		#print $t[0],"\t";	
		if ($info{$t[0]} ) {
			print $info{$t[0]},"\n";
		}
		else {
			print "$t[0]\n"
			#die ("No gene public name was found for $t[0]");
		}
	}
	close DICT;
}
##############################################################
#gene pulic id to worm id
##############################################################
sub p2w{
	open(DICT,"$mydict" ) or die("Can't read file: $ARGV[1] $!\n");
	my %info; my %transcripts;
	my %info2;
	#read the information
	while (<DICT>) {
		my @t=split "\t",$_;
		if ($t[0] && $t[1]) {
			$info{$t[1]}	= $t[0];
		}
		
	}
	#map the information
	while (<STDIN>) {
		#fetch the information
		my @t=split;	
		if ($info{$t[0]} ) {
			print $info{$t[0]},"\n";
		}
		else {
			print "$t[0]\n"
			#die ("No gene public name was found for $t[0]");
		}
	}
	close DICT;
}
sub p2t{
	open(DICT,"$mydict" ) or die("Can't read file: $ARGV[1] $!\n");
	my %info; my %transcripts;
	my %info2;
	#read the information
	while (<DICT>) {
		my @t=split "\t",$_;
		if ($t[1] && $t[5]) {
			$info{$t[1]}	.= $t[5]."\n";
		}
		
	}
	#map the information
	while (<STDIN>) {
		#fetch the information
		my @t=split;	
		if ($info{$t[0]} ) {
			print $info{$t[0]};
		}
		elsif ($t[0] eq 'F44E2.6') {
			print "F44E2.6\n";
		}
		else {
			die ("No gene public name was found for $t[0]");
		}
	}
	close DICT;
}
sub p2s{
	open(DICT,"$mydict" ) or die("Can't read file: $mydict $!\n");
	my %info; 
	#read the information
	while (<DICT>) {
		my @t=split "\t",$_;
		if ($t[1] && $t[2]) {
			$info{$t[1]}	= $t[2]."\n";
		}
		
	}
	#map the information
	while (<STDIN>) {
		#fetch the information
		my @t=split;	
		if ($info{$t[0]} ) {
			print $info{$t[0]};
		}
		elsif ($t[0] eq 'F44E2.6') {
			print "F44E2.6\n";
		}
		else {
			die ("No gene sequence name was found for $t[0]");
		}
	}
	close DICT;
}
sub net{
	open(DICT,"$mydict" ) or die("Can't read file: $ARGV[1] $!\n");
	open(BIOGRID,"$ARGV[2]" ) or die("Can't read file: $ARGV[2] $!\n");
	my $allid;
	while (<STDIN>) {$allid .= $_; }
	my %info;
	#read the information
	while (<DICT>) {
		my @t=split "\t",$_;
		if ($t[2] && $t[5]) {
			$info{$t[2]}	.= $t[5]."_";
		}
		
	}
	print $allid;
	#map the information
	while (<BIOGRID>) {
		#fetch the information
		my @t=split;	
		if ($info{$t[0]} && $info{$t[1]}) {
			my @id1=split "_", $info{$t[0]};
			my @id2=split "_", $info{$t[1]};
			foreach my $tmpid1(@id1) {
				if ($allid =~ /$tmpid1/) {
					foreach my $tmpid2(@id2) {
						if ($allid =~ /$tmpid2/) {
							print "$tmpid1\t $tmpid2\t$t[2]\n";
						}
					}
				}
			}
		}
	}
	close DICT;
	close BIOGRID;
}
#############################################################
#search overlap of two regions
##############################################################
sub overlap { # could be num_of_distance bases away
	my ($polar1,$start1,$end1,$polar2,$start2,$end2,$distance,$reverse) = @_;	
	#check if on the same strand
	unless ($reverse) {
		if ($polar1 ne $polar2) { return 0;	} 
	}
	#check if the on different strands
	else {
		if ($polar1 eq $polar2) { return 0;	} 
	}

	#check if overlaped within $distance away
	if ($start1-$distance < $end2 && $end1+$distance > $start2) {
		return 1;
	}	
	else { return 0;}
}

sub bins{
	my $size=$ARGV[1];
	while (<STDIN>) {
		my @t=split;
		my $chr=$t[0];
		my $start=$t[1];
		my $end=$t[2];
		for (my $i=$start;$i<=$end-$size;$i+=$size) {
			print "chr$chr\t$i\t",$i+$size,"\n";
		}
	}	
}

################################################################################
# seqsgr()
################################################################################
sub seqsgr {
	open(SGR,$ARGV[1]);
	my %wigsgr; my %segsgr;
	my (@chr,@start,@sgr);
	my $linei=0;
	while(<SGR>){
		if($linei==0){$linei++;next;}
		my @t=split;
		my $tmpchr=$t[0];
		push @chr,$tmpchr;
		push @start,$t[1];
		push @sgr,$t[2];
	}
	close SGR;
	
	$linei=0;
	my $totalnum=scalar @chr;
	foreach(@chr) {
		if($linei+1 == $totalnum) {last;}
		if($chr[$linei] eq $chr[$linei+1] && $sgr[$linei] != 0) {
			for (my $i=$start[$linei];$i<$start[$linei+1];$i++) {
				my $name=$chr[$linei]."\t".$i;
				$wigsgr{$name}=$sgr[$linei];
			}
		}	
		$linei++;
	}

	while(<STDIN>) {
		(my $name=$_) =~ s/\n//;
		my @t=split;
		my $chr=$t[0];
		my $start=$t[1];
		my $end=$t[2];
		for(my $i=$start; $i<=$end;$i++) {
			my $wigname=$chr."\t".$i;
			if ($wigsgr{$wigname}) {
				$segsgr{$name}+=$wigsgr{$wigname};
			}
		}
		unless ($segsgr{$name}) {	$segsgr{$name}=0;	}
	}
	foreach (keys %segsgr) {
		print $_,"\t",$segsgr{$_},"\n";
	}
}

sub ortholog {
	#read the information
	my @map=<STDIN>;
	open(ID,"$ARGV[1]" ) or die("Can't read file: $ARGV[1] $!\n");
	while (<ID>) {
		#fetch the information
		my @t=split;
		foreach (@map) {
			if ($_ =~/$t[0]/) {
				print $_;
			}
		}
	}
	close ID;
}

sub ReadWigVariable {
	open(WIG,$ARGV[1]);
	my %wigsgr; 
	my (@chr,@start,@sgr);
	my $chr; my $span;
	while(<WIG>){
		if ($_ =~ /^variablestep/i ) {
			my @t=split;
			$t[1]=~ /chrom\=(.*)/;
			$chr=$1;
			$chr=~ s/chr//;
			$t[2]=~ /span\=(.*)/;
			$span=$1;
		}
		elsif ($chr && !($_=~/^#/)) {
			my @t=split;
			for (my $i=$t[0];$i<$t[0]+$span;$i++) {
				my $name=$chr."\t".$i;
				if ($t[1]>0) {
					if ($wigsgr{$name}) {	$wigsgr{$name}+=$t[1]; }
					else {	$wigsgr{$name}=$t[1]; }
				}
			}
		}
	}
	close WIG;
	return \%wigsgr;	
}
sub seqwigv {
	my $wig=&ReadWigVariable;
	&outsgr($wig);
}
sub outsgr {
	my ($wigsgr) = @_;
	my $linei=0;
	while(<STDIN>) {
		my @t=split;
		my $tmpchr=$t[0];
		$tmpchr =~ s/chr//;
		my $start=$t[1];
		my $end=$t[2];
		my $allreads=0;
		for(my $i=$start; $i<=$end;$i++) {
			my $wigname=$tmpchr."\t".$i;
			if ($$wigsgr{$wigname}) {
				$allreads+=$$wigsgr{$wigname};
			}
		}
		print "chr$tmpchr\t$start\t$end\t",$allreads/($end-$start+1);
		print "\n";	
	}
}

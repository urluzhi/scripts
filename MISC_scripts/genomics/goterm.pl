#!/usr/bin/perl 
use strict;	use warnings;	# CREATED: Zhi John Lu 01/20/2009 13:55:44 EST
#===============================================================================
my $usage = <<EOF;
OPTIONS:
	-o (clean): format output (have >= 3rd level) of gostat with overlap
	-pie: output >=3rd level with lowest p-vlaue without overlap
	-o3 (clean): output only 3rd level of go term with overlap
	-pie3 : output only 3rd level of term for pie without overlap
EOF
#options:
my %option_hash = (-pie=>\&pie, -o=>\&output, -pie3=>\&pie3,-o3=>\&output3);
##############################################################
die $usage unless (@ARGV); my $option = ($ARGV[0]);
die("Unknown option: $option\n$usage\n") unless (defined($option_hash{$option}));
&{$option_hash{$option}};
exit;

#subroutines
##############################################################
# format the >= 3rd levels output for gostat (could having overlap)
##############################################################
sub output {
	print "ID\tName\tLevel\tp-value\tCount\tDatabase-Count";
	unless ($ARGV[1]) {print "\tFullName\tGene";}
	print "\n";
	while (<stdin>) {
		my @t=split "[\t\n]", $_;
		my $pvalue=$t[3];		#if (abs($pvalue) > 0.05) { next;}
		#pick out the name of go at level 3
		$t[0]=~/(GO:.*?),(.*)/;
		my $goid=$1;
		my ($level,$allname,$lastname);
		my @goterms= split ";",$2;
		foreach (@goterms) {
			my $tmplevel = ($_=~tr/%//);
			if (!$level || scalar $level < $tmplevel) {
				$allname=$_;
				$_=~/.*%(.*)$/;
				$lastname=$1;
				$level=$tmplevel;	
			}
		}
		my @num=split "#", $t[2];
		if ($lastname) {
			print "$goid\t$lastname\t$level\t$pvalue\t$num[0]\t$num[1]";
			unless ($ARGV[1]) {print "\t$allname\t$t[1]";}
			print "\n";
		}
	}
}
##############################################################
# output significant 3rd level terms (could having overlap)
##############################################################
sub output3 {
	print "ID\tName\tLevel\tp-value\tCount\tDatabase-Count";
	unless ($ARGV[1]) {print "\tFullName\tGene";}
	print "\n";
	while (<stdin>) {
		my @t=split "[\t\n]", $_;
		my $pvalue=$t[3];		#if (abs($pvalue) > 0.05) { next;}
		#pick out the name of go at level 3
		$t[0]=~/(GO:.*?),(.*)/;
		my $goid=$1;
		my ($allname,$lastname);
		my @goterms= split ";",$2;
		foreach (@goterms) {
			if ( ($_=~tr/%//)==3) {
				$allname=$_;
				$_=~/.*%(.*)$/;
				$lastname=$1;
				last;	
			}
		}
		my @num=split "#", $t[2];
		if ($lastname) {
			print "$goid\t$lastname\t3\t$pvalue\t$num[0]\t$num[1]";
			unless ($ARGV[1]) {print "\t$allname\t$t[1]";}
			print "\n";
		}
	}
}

##############################################################
# pie the >=3rd level
##############################################################
sub pie {
	my (%p4go,%name4go,%go4gene,%allname4go,%level4go);
	while (<stdin>) {
		my @t=split "[\t\n]", $_;
		$t[0]=~/(GO:.*?),(.*)/;
		my $goid=$1;
		#pick out the name of deepest level
		my @goterms= split ";",$2;
		foreach (@goterms) {
			my $tmplevel = ($_=~tr/%//);
			if (!$level4go{$goid} || scalar $level4go{$goid} < $tmplevel) {
				$allname4go{$goid}=$_;
				$_=~/.*%(.*)$/;
				$name4go{$goid}=$1;
				$level4go{$goid}=$tmplevel;
			}
		}
		#skip if level 3 not found
		next unless ($name4go{$goid});
		$p4go{$goid}=$t[3];	#if (abs($pvalue) > 0.05) { next;}
		#assign to gene if pvalue is the smallest
		my @genes=split ";",$t[1];
		foreach (@genes) {
			if(!$go4gene{$_}){
				$go4gene{$_}=$goid;
			}
			else {
				my $tmpid=$go4gene{$_};
				if( abs($p4go{$tmpid}) > abs($p4go{$goid}) ) {
					$go4gene{$_}=$goid;
				}	
			}
		}
	}
	my (%genes4go, %num4go);
	foreach  (keys %go4gene) {
		$genes4go{$go4gene{$_}}.= $_.",";
		$num4go{$go4gene{$_}}++;
	}
	print "ID\tName\tLevel\tp-value\tCount\tFullName\tGenes\n";
	foreach (keys %genes4go) {
		print "$_\t$name4go{$_}\t$level4go{$_}\t$p4go{$_}\t$num4go{$_}\t$allname4go{$_}\t$genes4go{$_}\n";
	}
}
##############################################################
# pie the 3rd level
##############################################################
sub pie3 {
	my (%p4go,%name4go,%go4gene,%allname4go);
	while (<stdin>) {
		my @t=split "[\t\n]", $_;
		$t[0]=~/(GO:.*?),(.*)/;
		my $goid=$1;
		#pick out the name of go at level 3
		my @goterms= split ";",$2;
		foreach (@goterms) {
			if ( ($_=~tr/%//)==3) {
				$allname4go{$goid}=$_;
				$_=~/.*%(.*)$/;
				$name4go{$goid}=$1;
				last;	
			}
		}
		#skip if level 3 not found
		next unless ($name4go{$goid});
		$p4go{$goid}=$t[3];	#if (abs($pvalue) > 0.05) { next;}
		#assign to gene if pvalue is the smallest
		my @genes=split ";",$t[1];
		foreach (@genes) {
			if(!$go4gene{$_}){
				$go4gene{$_}=$goid;
			}
			else {
				my $tmpid=$go4gene{$_};
				if( abs($p4go{$tmpid}) > abs($p4go{$goid}) ) {
					$go4gene{$_}=$goid;
				}	
			}
		}
	}
	my (%genes4go, %num4go);
	foreach  (keys %go4gene) {
		$genes4go{$go4gene{$_}}.= $_.",";
		$num4go{$go4gene{$_}}++;
	}
	print "ID\tName\tLevel\tp-value\tCount\tFullName\tGenes\n";
	foreach (keys %genes4go) {
		print "$_\t$name4go{$_}\t3\t$p4go{$_}\t$num4go{$_}\t$allname4go{$_}\t$genes4go{$_}\n";
	}
}
##############################################################
#deprecated
# assign go term to each gene (at least 3 level, lowest p value)
##############################################################
sub assign {
	my (%pvalue,%go,%level,%name,%id4gene);
	my $start_annot=0;
	while (<STDIN>) {
		#read the deepest name and pvale for each go term
		if ($_ =~ /^GO:/) {
			my @t=split "[\t\n]", $_;
			$t[0]=~/(GO:.*?),/;
			my $id=$1;
			$pvalue{$id}=$t[3];
			$go{$id}=$t[0];
			my @branches = split ";",$t[0];
			foreach (@branches) {
				my @terms=split "%",$_;
				if (!$level{$id} || scalar @terms > $level{$id}) {
				   	$level{$id}=scalar @terms -1;
					$name{$id}=$terms[-1];
				}
			}
		}
		#assign annotation with lowest p-value to 
		elsif($_ =~/# Table with Gene annotations:/) {$start_annot=1;next;}
		elsif($start_annot) {
			my @t=split "[\t\n]", $_;
			if ($t[2]) {
				my $fixed=10;
				my $fixedid;
				my @g=split ";;;",$t[2];
				foreach (@g) {	
					$_=~/(GO:.*?),/;
					my $id=$1;
					if ($pvalue{$id} && $pvalue{$id} < $fixed)	{
					   	$fixedid=$id;
						$fixed=$pvalue{$id};
					}
				}
				if ($fixedid) {$id4gene{$t[1]}=$fixedid;}
			}
		}
	}
	foreach (keys %id4gene) {
		my $goid=$id4gene{$_};
		print "$_\t$goid\t";
		print "$name{$goid}\t$level{$goid}\t",$pvalue{$goid},"\n";
	}
}
##############################################################
#deprecated
# group for pie for the deep level with highest pvalue
##############################################################
sub grp {
	my (%hash,%pvalue,%go);
	while (<stdin>) {
		my @t=split "[\t\n]", $_;
		$pvalue{$t[2]}=$t[4];
		$hash{$t[2]}++;
		$go{$t[2]}=$t[1];
	}
	foreach (keys %hash) {
		print "$_\t$hash{$_}\t$go{$_}\t";
		if ($pvalue{$_}) {
			print abs($pvalue{$_}),"\t";
			if ($pvalue{$_} >0) { print "+\n";}
			else {print "-\n";}
		}
		else {print "1\n";}

	}
}





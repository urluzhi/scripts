#!/usr/bin/perl 
use strict;	use warnings;	# CREATED: Zhi John Lu 04/07/10 15:13:33 EDT
#===============================================================================
my $dict_dir="/home9-2/gersteinlab/common/modencode/RNAseq/polyARNA/Waterston/dcpm/";
my $usage = <<EOF;
USAGE:
	./TFtargets.expression.pl [TARGETS] TF_Stage 
EOF
die $usage unless (@ARGV);

#get the stage
$ARGV[1]=~/_(.*)/;
my $dict;
if($1 eq 'Late_EMB') { $dict="EE";}
elsif ($1 eq 'Early_EMB')	{	$dict="LE";	}
else 	{$dict=$1;}
$dict=$dict_dir.$dict.".wormbase.dcpm.sorted";
#get expresion for each transcript 
my %info;
open(DICT,"$dict" ) or die("Can't read file: $dict $!\n");
while (<DICT>) {
  my @t=split; 
  if ($t[1] && $t[0]) {
    $info{$t[0]}    = $t[1];
  }
}
close DICT;
#get the transcript id and expression for it
##############################################################
open (TARGETS, "$ARGV[0]" );
my $i=0;
while (<TARGETS>) {
	$_=~ s/\n//;
	print $_,"\t";
	if($i==0 && $_=~ /^Peak/) {	$i=1; print "Expression_of_Coding_Transcripts\tExpression_Median\n";next;}	
	my @tab=split;
	if($tab[6]) {
		my @expression;
		my @trans=split /,/, $tab[6];
		foreach(@trans) {
			$_ =~ /(.*)\(.*/;
			if($1 && $info{$1}) { 	push @expression,$info{$1};	}
			else {	push @expression,'na';	}
		}
		$"=",";
		print  "@expression";
		print "\t",median(\@expression);
	}
	print "\n";	
}
close(TARGETS);
exit;

#median
##########
sub median {
	my $rpole = shift;
	my @tmppole = @$rpole;
	my @pole;
	foreach (@tmppole) {
		if($_ ne 'na') {	push @pole,$_;}
	}
	return 'na' unless (@pole);
	my $ret;
	@pole=sort(@pole);
	if( (@pole % 2) == 1 ) {
		$ret = $pole[((@pole+1) / 2)-1];
	} 
	else {
		$ret = ($pole[(@pole / 2)-1] + $pole[@pole / 2]) / 2;
	}
	return $ret;
}



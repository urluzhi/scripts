#!/usr/bin/perl 
use strict;	use warnings;	# CREATED: Zhi John Lu 01/20/2009 13:55:44 EST
#===============================================================================
my $usage = <<EOF;
OPTIONS:
	-m: merge for GSEA

EOF
#options:
my %option_hash = (-m=>\&merge, -o=>\&output);
##############################################################
die $usage unless (@ARGV); my $option = ($ARGV[0]);
die("Unknown option: $option\n$usage\n") unless (defined($option_hash{$option}));
&{$option_hash{$option}};
exit;

#subroutines
##############################################################
# main
##############################################################
sub merge {
	#open(IN,"$ARGV[1]" ) or  die("Can't read file: $ARGV[1] $!\n");
	my %hash;
	my $NR=1;
	print "Sequence_Name\tDescription\tLog2Ratio\tBackground\n";
	while (<STDIN>) {
		if ($NR==1) {$NR++;next;} #skip the first line
		my @t=split "\t", $_;
		if (!$hash{$t[0]} || abs($hash{$t[0]}) < abs($t[2])) {
			$hash{$t[0]}=$t[2];
		}			
		$NR++;
	}
	foreach (keys %hash) {
		print "$_\tna\t$hash{$_}\t0\n";
	}
}


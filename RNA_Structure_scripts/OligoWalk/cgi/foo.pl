#!/usr/bin/perl 
use strict;	use warnings;	# CREATED: Zhi John Lu 09/03/2008 08:08:53 PM EDT
#===============================================================================
my $usage = <<EOF;
OPTIONS:
	-a: run main  
EOF
if (!$ARGV[0]) {	die $usage;	}
elsif ($ARGV[0] eq "-a" && $ARGV[1]) {	&main;	}
else {print "\n\t\tUnrecognized options!!!\n\n\n";	die $usage;}
exit;
#subroutines
#===============================================================================
sub main {
	open(IN,"$ARGV[1]" ) or  die("Can't read file: $ARGV[1] $!\n");
	while (<IN>) {
		my @tab=split "\t", $_;
		push @tab, $_;
	}
	close IN;
}


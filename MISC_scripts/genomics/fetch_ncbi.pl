#!/usr/bin/perl 
use strict;	use warnings;	# CREATED: Zhi John Lu 09/10/09 13:28:43 EDT
#===============================================================================
my $usage = <<EOF;
OPTIONS:
	-a: run main  
EOF
#options:
my %option_hash = (-a=>\&main, -o=>\&output);
##############################################################
die $usage unless (@ARGV); my $option = ($ARGV[0]);
die("Unknown option: $option\n$usage\n") unless (defined($option_hash{$option}));
&{$option_hash{$option}};
exit;

#subroutines
##############################################################
# main
##############################################################
sub main {
	open(IN,"$ARGV[1]" ) or  die("Can't read file: $ARGV[1] $!\n");
	while (<IN>) {
		my @tab=split "\t", $_;
		push @tab, $_;
	}
	close IN;
}


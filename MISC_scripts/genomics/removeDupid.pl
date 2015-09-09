#!/usr/bin/perl 
use strict;	use warnings;	# CREATED: Zhi John Lu 01/20/2009 13:55:44 EST
#===============================================================================
my $usage = <<EOF;
OPTIONS:
	-m: merge same ids

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
	while (<STDIN>) {
		my @t=split;
		if (!$hash{$t[0]})  {$hash{$t[0]}=1;}			
		else { $hash{$t[0]}++;	}
	}
	foreach (keys %hash) {
		if ($hash{$_} >0) {
			print "$_\t$hash{$_}\n";
		}
	}
}


#!/usr/bin/perl 
use strict;
use warnings;
#############################################################################
	open(REF,"$ARGV[0]" ) or die("Can't read file: $ARGV[0] $!\n");
	my %info;
	#read the information
	while (<REF>) {
		my @t=split;
		my $id=$t[0]."_".$t[1]."_".$t[2];
		$info{$id}=$t[4] ;
	}
	close REF;
	#map the information
	while (<STDIN>) {
		#fetch the information
		my @t=split;
		print $t[4],"\t";
		if ($t[8])	{
			my @oligo=split ",",$t[8];
			foreach (@oligo) {
				$_=~/^(\d+)-(\d+)$/;
				my $id=$t[0]."_".$1."_".$2;
				print "$info{$id},"
			}
		}
		print "\n";
	}


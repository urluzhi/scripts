#!/usr/bin/perl 
use strict;
use warnings;
use GENOME;
#############################################################################
#print out the usage inforagion
my $usage = <<EOF;
USAGE:
	overlap.pl [options] (PATTERN) [FILE1 FILE2 ...] 
OPTIONS:
	-frac coverage 0/1/2 <Map>: overlap limited by frac
		coverage: 0.0-1.0
		0: both segments
		1: percentage is from first segment (chr,start,end)
		2: percentage is from second segment (start-end,start-end)
	-size bps <Map>: overlap limited by # of base pairs
EOF
#options:
my %option_hash=(-size=>\&size, -frac=>\&frac	);
#outpust usage information
die $usage unless (@ARGV);	my $option = ($ARGV[0]);
die("Unknown option: $option\n$usage\n") unless (defined($option_hash{$option}));
&{$option_hash{$option}};
exit;


################################################################################
################################################################################
sub size{
	my $annotation;
	#map the information
	while (<STDIN>) {
		#fetch the information
		my @t=split;
		$annotation=0;
		if ($t[3])	{
			my @oligo=split ",",$t[3];
			foreach (@oligo) {
				$_=~/^(\d+)-(\d+)$/;
				#if coverage of query satisfied (first seg.,gff)
				if ( &calc_coverage($t[1],$t[2],$1,$2)  > $ARGV[1] ) {
					$annotation=1;
				}
			}
		}
		print $annotation,"\n";
	}
}
sub frac{
	my $annotation;
	#map the information
	while (<STDIN>) {
		#fetch the information
		my @t=split;
		$annotation=0;
		if ($t[3])	{
			my @oligo=split ",",$t[3];
			foreach (@oligo) {
				$_=~/^(\d+)-(\d+)$/;
				#if coverage of query satisfied (first seg.,gff)
				if ( (&calc_coverage($t[1],$t[2],$1,$2)/($t[2]-$t[1]+1) ) > $ARGV[1] && $ARGV[2]=~/[10]/ ) {
					$annotation=1;
				}
				#if coverage of target satisfied (second seg.,map)
				elsif ( (&calc_coverage($t[1],$t[2],$1,$2)/($2-$1+1) ) > $ARGV[1] && $ARGV[2] =~/[20]/ ) {$annotation=1;}
			}
		}
		print $annotation,"\n";
	}
}


#############################################################
#caculate overlaped size
##############################################################
sub calc_coverage { # could be num_of_distance bases away
	my ($s1,$e1,$s2,$e2) = @_;	
	my @tmp= ($e1-$s2,$e2-$s1,$e1-$s1,$e2-$s2);
	return minNum(\@tmp)+1; 
}

sub minNum {
	my ($a)=@_;
	my $min=$$a[0];
	foreach (@$a) {
		if($min>$_) {$min=$_;}
	}
	return $min;
}

sub log2 {
		my $n = shift;
		return (log($n)/log(2));
}


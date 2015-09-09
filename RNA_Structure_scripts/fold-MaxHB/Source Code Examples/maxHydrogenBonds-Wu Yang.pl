=head1 Usage: 

	Name: maxHydrogenBonds-DP.pl

	Predict RNA structure with the maxmimum hydrogen bonds following such critera.
	Bases pair included:
		GC or CG pair, score = 3.
		AU or UA pair, score = 2.
		GU or UG pair, score = 2.
	Minimum length of internal loops is set to be 3.

	Example:  perl maxHydrogenBonds-DP.pl
		  perl maxHydrogenBonds-DP.pl [sequence] >outfile

	Author:  Yang Wu (wuyang.bnu@139.com)

=cut

#! /usr/bin/perl -w
use strict;
use warnings;

unless (@ARGV) {
    &Usage();  # print help information
    exit;
}

# dealing with the input
my $seq = $ARGV[0];
$seq =~ tr/Tt/Uu/;  # transform to RNA sequence
my @s = split("", $seq);
print "Script: $0\nSequence: $seq\nLength: ", length($seq), "\n";  # general information


# declaration
my $min = 3;   # minimum length of internal loops: $min
my %V;  # max number of H-bonds when i and j are paired
my %W;  # max number of H-bonds
my %T;  # for traceback


# Filling %V, %W and %T using dynamic programming
&fill(@_);

# Traceback to get the optimal structure
print "$seq\n";
&traceback(1, scalar(@s), \%T);
print "\n\n";

# Output tables 
&getTable(1, scalar(@s), \%V, 'V(i,j)');
&getTable(1, scalar(@s), \%W, 'W(i,j)');
#&getTable(1, scalar(@s), \%T, 'T(i,j)');
exit;    


###############################################################################################################
sub fill {
    my $len = scalar(@s);  # length of the whole sequence
    for my $l (1 .. $len) {  # length of the subsequences
	for my $i (1 .. $len-$l+1) {  # start of a subsequence
	    my $j = $i + $l - 1;  # end of a subsequence
	    my $pair = $s[$i-1].$s[$j-1];  # the bases in the end of the subsequcnes
	    #print "length=$l: $i, $j, $pair\n";  # checking
	
	    # fill in V(i,j)
	    if (($j-$i) > $min && $pair =~ /AU|UA|GU|UG|CG|GC/i) {  # canonical pair with a sufficient length
		my $H_bonds = ($pair =~ /CG|GC/i) ? 3 : 2;
		$V{$i}{$j} = $H_bonds + $W{$i+1}{$j-1};
	    } else {  # not paired or too short
		$V{$i}{$j} = 0;
	    }

	    # fill in W(i,j)  # max{ V(i,j), W(i,k)+W(k+1,j) (i<=k<=j-1)}
	    $W{$i}{$j} = $V{$i}{$j};
	    $T{$i}{$j} = 0;
	    for my $k  ($i .. $j-1) {  # search for maxmimum branch
		if ($W{$i}{$k} + $W{$k+1}{$j} > $W{$i}{$j}) {
		    $W{$i}{$j} = $W{$i}{$k} + $W{$k+1}{$j};
		    $T{$i}{$j} = $k;
		}
	    }

	    # fill in T(i,j)
	    $T{$i}{$j} = ($T{$i}{$j} != 0) ? $T{$i}{$j} :  #  k if W(i,j) == W(i,k) + W(k+1,j)
			 ($V{$i}{$j} == 0) ? -1 : -2;      # -1 if V(i,j) == 0
							   # -2 if i,j paired
	    #print "\t=>W=$W{$i}{$j}, T=$T{$i}{$j}\n";  # checking
	}
    }
}

###############################################################################################################
sub traceback {
    my ($i, $j, $ref) = @_;
    if ($i == $j) {  # minimum length
	print '.';
    } elsif ($j > $i) {  # possible structures
	if ($ref->{$i}{$j} == -2) {  # i,j paired
	    print "(";
	    &traceback($i+1, $j-1, $ref);
	    print ")";
	} elsif ($ref->{$i}{$j} == -1) {  # V(i,j) == 0
	    print ".";
	    &traceback($i+1, $j-1, $ref);
	    print ".";
	} else {  # W(i,j) == W(i,k) + W(k+1,j)
	    &traceback($i, $ref->{$i}{$j}, $ref);
	    &traceback($ref->{$i}{$j}+1, $j, $ref);
	}
    }
}

###############################################################################################################
sub getTable {
    my ($s, $e, $ref, $name) = @_;
    my $width = scalar(@s)*2-1;
    
    print "\nTable $name:\nj-> ", join(" ", @s), "\n    ", '-' x $width, "\n";
    for my $i (1 .. scalar(@s)) {
	my $sign = ($i==2) ? 'i' :
	           ($i==3) ? '|' :
		   ($i==4) ? 'v' : " ";
	print "$sign ", $s[$i-1], "|";
	for my $j (1 .. scalar(@s)) {
	    my $score = (exists $ref->{$i}{$j}) ? $ref->{$i}{$j} : " ";
	    print " $score";
	}
	print "\n";
    }
    print "    ", '-' x $width, "\n";
}

###############################################################################################################
sub Usage {
    print "perl maxHydrogenBonds.pl [sequence] >outfile
	
	Example: perl maxHydrogenBonds-DP.pl GCGGGUACCGAUCGUCGC

	Predict RNA structure with the maxmimum hydrogen bonds following such critera.
	Bases pair included:
		GC or CG pair, score = 3.
		AU or UA pair, score = 2.
		GU or UG pair, score = 2.
	Minimum length of internal loops is set to be 3.
	
	Yang Wu \@ Lu lab, Tsinghua University
	wuyang.bnu\@139.com, all rights reserved\n";
}

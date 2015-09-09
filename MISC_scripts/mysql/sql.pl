#!/usr/bin/perl 
#===============================================================================
#
#         FILE:  mysql.pl
#
#        USAGE:  ./mysql.pl 
#
#  DESCRIPTION: Manipulate mysql database 
#
#        NOTES:  ---
#       AUTHOR:  Zhi John Lu <urluzhi@gmail.com>
#      COMPANY:  University of Rochester, Medical Center
#      VERSION:  1.0
#      CREATED:  07/26/2007 14:05:11 EDT
#     REVISION:  ---
#===============================================================================

use strict;
use warnings;
use DBI;



#############################################################################
#print out the usage inforagion
my $usage = <<EOF;

USAGE:
	 mysql.pl [options] [FILE1 FILE2 ...] 
DESCRIPTION:
	Testing the zscore of 5sRNA and tRNA.
OPTIONS:
	-h,--help
		Usage information
	-c DATABASE TABLE
		create a database table structure with mysql
	-w DATABASE TABLE [FILE.in] 
		fetch information in a file to a database table
	-r DTABASE TABLE
		read information in a database table with mysql

Created: Jul. 26, 2007
by John(Zhi Lu)            zhi_lu\@urmc.rochester.edu

EOF

die $usage unless (@ARGV && $ARGV[0] ne "-h" && $ARGV[0] ne "--help");

# Connect to the database.
#############################################################################
our $database= $ARGV[1];
our $table = $ARGV[2];
our $dbh = DBI->connect("DBI:mysql:$database;host=localhost", "zl222", "1234", {'RaiseError' => 1});
#   my $dbh = DBI->connect("DBI:mysql:$database;mysql_read_default_file=/etc/my.cnf;host=localhost", "john", "", {'RaiseError' => 1});

##############################################
#different usages
if ($ARGV[0] eq "-c") {
	&create_mysql; 
}
elsif ($ARGV[0] eq "-w") {
	&write_mysql; 
}
elsif ($ARGV[0] eq "-r") {
	&read_mysql; 
}

# Disconnect from the database.
$dbh->disconnect();
exit;

#subroutines
#####################################################################

sub write_mysql {
	
	my ($fh) = open_file($ARGV[3]);
	my $i=1;
	while (<$fh>) {
#		$_=~ s/\t//g;
#		my @tab=split /\|\|/,$_;
		my @tab=split " ",$_;
		my $seq_file1="sequences/".$tab[3]."__".$tab[0].".seq";

		my %seq1 = read_seq($seq_file1);

#		print "$seq1{'length'},$seq2{'comments'}";
	#	print "$tab[0], $tab[1], $tab[2]\n";
	#       1        2      345 678     9    10    11     12   13 14    15       16        17      
	# -1  energy;min_length;AUC;AUC  start1 end1 polar1 start2 e2 p2 prob_dyn prob_rnaz class_qrna
#		$_ =~ /1:([-.\d]+) 2:([-.\d]+) 3:([-.\d]+) 4:([-.\d]+) 5:([-.\d]+) 6:([-.\d]+) 7:([-.\d]+) 8:([-.\d]+)\t([.\d]+) ([.\d]+) ([-+]) ([.\d]+) ([.\d]+) ([-+]) ([-.\de]+)\t([-.\de]+)\t(\w+)/;
#		$_=~ /overlaps ncRNA \[Ecoli:(.*)\] \[Styphi:(.*)\]/;
#		my $seq1= $1;	my $seq2=$2;
#		while ($seq1 =~ /(\d+),(\d+),([-+]),(\w+)/g) {	print "$3 $1 $2 $4\t";	}
#		while ($seq2 =~ /(\d+),(\d+),([-+]),(\w+)/g) {	print "$3 $1 $2 $4\t";	}
		
#		$dbh->do("update $table set length1=$seq1{'length'}, length2=$seq2{'length'}, seq1='$seq1{'nuc'}', seq2='$seq2{'nuc'}'   where id = $i ");
#		$dbh->do( "insert into $table (ave,sd,zscore) values ($tabs[0],$tabs[1],$tabs[2])" );
#		print "INSERT INTO $table VALUES ($i,'$11',$9,$10,'$14',$12,$13,$3*100,$5*100,(1-$3-$4-$5)*100,$4*100, $6*100,$8*100,(1-$6-$7-$8)*100,$7*100,NULL,NULL,$2,NULL,NULL,$1/10,NULL,$15,$16,'$17')\n";
#		$dbh->do("insert into $table values (NULL, '$tab[0]',$tab[1], $tab[2], '$tab[3]', $tab[4], $tab[5], $tab[6], $tab[7], $tab[8], $tab[9], $tab[10], $tab[11], $tab[12], $tab[13], $tab[14], 100*$tab[15], $tab[16], $tab[17], $tab[18], $tab[19], $tab[20], $tab[21],NULL,NULL,NULL,NULL,NULL,NULL,NULL) ");

		$i++;
	}	
	
	close $fh;
}


sub read_mysql {

  # Now retrieve data from the table.
  my $sth = $dbh->prepare("SELECT * FROM $table");  $sth->execute();

  while (my $ref = $sth->fetchrow_hashref()) {
#    print "id = $ref->{'id'}, position: $ref->{'position'}, zscore = $ref->{'zscore'}\n";
    print "$ref->{'id'}\n";
  }

  $sth->finish();
}


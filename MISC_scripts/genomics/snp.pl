#!/usr/bin/perl 
use strict;	use warnings;	# CREATED: Zhi John Lu 06/16/10 11:30:54 EDT
#===============================================================================
my $usage = <<EOF;
OPTIONS:
	-extract [VCF.gz] [GFF] : count how many snp in GFF file  
	-count [VCF.gz] [GFF] : count how many snpi (averaged by genes) in GFF file (10,000 genes sampled)  
	-count_maf [VCF.gz] [GFF] : count the MAF for all SNPs in GFF file (10,000 genes sampled) 
EOF
#options:
my %option_hash = (-count=>\&count,
	   	-extract=>\&extract,
	   	-count_maf=>\&count_maf,
		-o=>\&output);
##############################################################
die $usage unless (@ARGV); my $option = ($ARGV[0]);
die("Unknown option: $option\n$usage\n") unless (defined($option_hash{$option}));
&{$option_hash{$option}};
exit;

#subroutines
##############################################################
# extract the overlapped SNPs with genes
##############################################################
sub extract {
	my (%out,%chr);
	#define auto chromosomes
	foreach(1 .. 22) { $chr{$_}=$_; };
	#read genes from gff file
	open(GFF,$ARGV[2]) or  die("Can't read file: $ARGV[2] $!\n");
	my @gff=<GFF>;
	close GFF;
	my $i=0;
	foreach(shuffle_array(@gff)) {
		last unless $i< 10000; # only sample 10,000 gene at most from autochromosomes
		my $gene=$_;
		my @gene_s=split;
		$gene_s[0]=~ s/chr//g;
		next unless defined $chr{$gene_s[0]};#only keep auto chromosomes
		my $cmd="/home1/zl222/bin/tabix/tabix ".$ARGV[1]." ".$gene_s[0].":".$gene_s[3]."-".$gene_s[4];
		open (OUT, "$cmd |");
		while (<OUT>) {	$out{$gene}=$_;}
		close OUT;
		$out{$gene}=0 unless ($out{$gene});
		$i++; # count genes from autochromosomes
	}
#	print %out;
	return %out;
}

##############################################################
# cout SNP density and MAF for each gene
##############################################################
sub count {
	print "Freq_of_SNP\tMAF\n";
	my %input=&extract;
	my @num_snps;
	my @maf;
	my @mygenes;
	my $i=0;
	#loop each gene
	foreach (keys %input) {
		$num_snps[$i]=0;
		$maf[$i]=0;
		# if the gene has SNPs
		if($input{$_}) {
			$mygenes[$i]=$_;
			my @gene_s=split;
			my @SNPs=split /\n/,$input{$_};
			$num_snps[$i]=scalar @SNPs;
			#loop each SNP to average MAF
			foreach my $snp(@SNPs) {
				my @snp_s=split /\t/,$snp;
				#minor allele frequency
				$snp_s[7]=~/AN=(\d+);/;
				my $AN=$1;	
				$snp_s[7]=~/AC=([\d,]+);/;
			   	my @AC=split /,/, $1;
				my $num_of_AC=scalar @AC;
				$AC[$num_of_AC]=$AN;
				for(my $j=0;$j<$num_of_AC;$j++) {	$AC[$num_of_AC]-=$AC[$j];}
				$maf[$i]+=min(@AC)/$AN;
			}
			#average SNP density and MAF for each gene
			$maf[$i]=$maf[$i]/$num_snps[$i];
			$num_snps[$i]=$num_snps[$i]/($gene_s[4]-$gene_s[3]+1);
		}
		$i++;
	}
	$i=0;
	#print out
	foreach(@num_snps) {
		print $_,"\t";
		print $maf[$i],"\t";
		print $mygenes[$i];
		$i++;
	}
}
##############################################################
# cout SNP density and MAF for all regions
##############################################################
sub count_maf {
	print "MAF\n";
	my %input=&extract;
	my @num_snps;
	my $i=0;
	#loop each gene
	foreach (keys %input) {
		$num_snps[$i]=0;
		# if the gene has SNPs
		if($input{$_}) {
			my @gene_s=split;
			my @SNPs=split /\n/,$input{$_};
			$num_snps[$i]=scalar @SNPs;
			#loop each SNP to average MAF
			foreach my $snp(@SNPs) {
				my @snp_s=split /\t/,$snp;
				#minor allele frequency
				$snp_s[7]=~/AN=(\d+);/;
				my $AN=$1;	
				$snp_s[7]=~/AC=([\d,]+);/;
			   	my @AC=split /,/, $1;
				my $num_of_AC=scalar @AC;
				$AC[$num_of_AC]=$AN;
				for(my $j=0;$j<$num_of_AC;$j++) {	$AC[$num_of_AC]-=$AC[$j];}
				#push @maf, min(@AC)/$AN;
				print min(@AC)/$AN,"\n";
				#print $snp_s[7],"\n";
			}
		}
		$i++;
	}
}

sub max {
	my ($max, @vars) = @_;
	for (@vars) {
		$max = $_ if $_ > $max;
	}
	return $max;
}
sub min {
	my ($min, @vars) = @_;
	for (@vars) {
		$min = $_ if $_ < $min;
	}
	return $min;
}

#shuffle the array's order
sub shuffle_array {
	my (@old) = @_;  
    my @new = ();
    for( @old ){
        my $r = rand (@new+1); 
		push(@new,$new[$r]);
        $new[$r] = $_;
    }
	return @new;
}


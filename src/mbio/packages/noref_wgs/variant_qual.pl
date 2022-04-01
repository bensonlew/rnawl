#!/usr/bin/env perl
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($snp_vcf,$depth,$qual,@sample,@head);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Cwd qw(abs_path);

my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$snp_vcf,
    "o1:s"=>\$depth,
	) or &USAGE;
&USAGE unless ($snp_vcf and $depth);
my @indi;
open In,$snp_vcf;
my %depth;
while (<In>) {
	chomp;
	next if ($_ eq "" ||/^$/ || /^##/);
	my ($chr,$pos,$id,$ref,$alt,$qual,$Filter,$indo,$format,@geno)=split(/\t/,$_);
	next if ($Filter ne "PASS" && $Filter ne "SNP" && $Filter ne "INDEL" && $Filter ne "FILTER");
	if (scalar @indi ==0) {
		push @indi,@geno;
		next;
	}else{
		my @Geno;
		my @format=split(/\:/,$format);
		for (my $i=0;$i<@geno;$i++) {
			my @info=split(/\:/,$geno[$i]);
			my $dp=-1;
			for (my $j=0;$j<@info;$j++) {
				if ($format[$j] eq "AD") {
					my ($d1,$d2)=split(/\,/,$info[$j]);
					next if ($d1 eq "." || $d2 eq ".");
					$dp+=$d1;
					$dp+=$d2;
				}
			}
			next if ($dp eq "-1");
			$depth{$indi[$i]}{$dp}++;
		}
	}
}
close In;
open Out,">$depth";
print Out "sampleID\tDepth\tnum\n";
my $sum=0;
foreach my $sample (sort keys %depth) {
	foreach my $depth (sort {$a<=>$b} keys %{$depth{$sample}}) {
		$sum+=$depth{$sample}{$depth};
		print Out $sample,"\t",$depth,"\t",$sum,"\n";
	}
}
close Out;
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################


sub USAGE {#
        my $usage=<<"USAGE";
Contact:        xiaohong.wu\@majorbio.com;
Script:			$Script
Description:
	Cumulative depth distribution
	Cumulative GQ Score distribution
eg:
	perl $Script -i vcf -o1  Cumulative_DEP.txt -o2  Cumulative_GQ.txt,

Usage:
  Options:
  -i	<file>	input vcf file
  -o1	<file>	output file of Cumulative_DEP.txt,
  -h         Help

USAGE
        print $usage;
        exit;
}

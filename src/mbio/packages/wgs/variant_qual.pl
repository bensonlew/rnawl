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
	"o2:s"=>\$qual,
	) or &USAGE;
&USAGE unless ($snp_vcf);
my @indi;

open In,$snp_vcf;
if ($snp_vcf=~/.gz$/) {
	close In;
	open In,"gunzip -c $snp_vcf|";
}
my %depth;
my %qual;
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
			my $GQ=-1;
			for (my $j=0;$j<@info;$j++) {
				if ($format[$j] eq "AD") {
					my ($d1,$d2)=split(/\,/,$info[$j]);
					next if ($d1 eq "." || $d2 eq ".");
					$dp+=$d1;
					$dp+=$d2;
				}
				if ($format[$j] eq "GQ") {
					next if ($info[$j] eq ".");
					$GQ=$info[$j];
				}
			}
			next if ($dp eq "-1" || $GQ eq "-1");
			$depth{$indi[$i]}{$dp}++;
			$qual{$indi[$i]}{$GQ}++;
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
open Out,">$qual";
print Out "sampleID\tGQvalue\tnum\n";
my $qualsum=0;
foreach my $sample (sort keys %qual) {
	foreach my $depth (sort {$a<=>$b} keys %{$qual{$sample}}) {
		$qualsum+=$qual{$sample}{$depth};
		print Out $sample,"\t",$depth,"\t",$qualsum,"\n";
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
  -o2   <file>  output file of Cumulative_GQ.txt,
  -h         Help

USAGE
        print $usage;
        exit;
}

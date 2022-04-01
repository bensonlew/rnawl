#!/usr/bin/env perl
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fIn,$fOut,$report,$fmtrix);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$fIn,
	"o:s"=>\$fOut,
	"m:s"=>\$fmtrix,
			) or &USAGE;
&USAGE unless ($fIn and $fOut and $fmtrix);
my %TsTv=(
	"AT"=>"Tv",
	"AG"=>"Ts",
	"AC"=>"Tv",
	"GT"=>"Tv",
	"CT"=>"Ts",
	"CG"=>"Tv",
);
open In,$fIn;
if ($fIn=~/.gz$/) {
	close In;
	open In,"gunzip -c $fIn|";
}
my @indi;
my %vcfstat;
my %diff;
while (<In>) {
	chomp;
	next if ($_ eq "" || /^$/ || /^##/);
	my ($chr,$pos,$id,$ref,$alt,$qual,$Filter,$indo,$format,@geno)=split(/\t/,$_);
	next if ($Filter ne "PASS" && $Filter ne "SNP" && $Filter ne "INDEL" && $Filter ne "FILTER");
	my @Ale=split(/\,/,join(",",$ref,$alt));
	if (scalar @indi ==0) {
		push @indi,@geno;
		next;
	}else{
		my @Geno;
		my %geno;
		my %pchange;
		my $sum=0;
		my @format=split(/\:/,$format);
		for (my $i=0;$i<@geno;$i++) {
			my ($gt,$ad,$dp);
			my @info=split(/\:/,$geno[$i]);
			for (my $j=0;$j<@info;$j++) {
				if ($format[$j] eq "GT") {
					$gt=$info[$j];
				}
				if ($format[$j] eq "AD") {
					$ad=$info[$j];
				}
				if ($format[$j] eq "DP") {
					$dp=$info[$j];
				}
			}
			push @Geno,$gt;
			my ($g1,$g2)=split(/\//,$gt);
			if ($gt eq "./." || $Ale[$g1] eq "*" || $Ale[$g2] eq "*"){
				$vcfstat{$indi[$i]}{miss}++;
				next;
			}
			$geno{$gt}++ if($gt ne "./.");
			if ($gt eq "0/0"){
				$vcfstat{$indi[$i]}{ref}++;
				next;
			}
			# $geno{$gt}++ if($gt ne "./.");
			$sum+=$dp if($gt ne "./.");
			$vcfstat{$indi[$i]}{total}++;
			$vcfstat{$indi[$i]}{dp}+=$dp;
			$vcfstat{$indi[$i]}{homo}++ if($g1 eq $g2);
			my %change;
			$change{$Ale[$g1]}=1;
			$change{$Ale[$g2]}=1;
			$change{$ref}=1;
			$pchange{$Ale[$g1]}++;
			$pchange{$Ale[$g2]}++;
			if (scalar keys %change >= 3) {
				$vcfstat{$indi[$i]}{multi}++;
			}else{
				my $change=join("",sort keys %change);
				$vcfstat{$indi[$i]}{$TsTv{$change}}++;
			}
		}
		$vcfstat{pop}{total}++ if(scalar keys %geno != 1 || scalar @indi ==1);
		$vcfstat{pop}{dp}+=$sum if(scalar keys %geno != 1 || scalar @indi ==1);
		if (scalar keys %pchange >=3) {
			$vcfstat{pop}{multi}++;
		}elsif (scalar keys %pchange !=2) {
			next;
		}else{
			my $change=join("",sort keys %pchange);
			if (!defined $vcfstat{pop}{$TsTv{$change}}) {
				$vcfstat{pop}{$TsTv{$change}}=1;
			}else{
				$vcfstat{pop}{$TsTv{$change}}++;
			}
		}
		for (my $i=0;$i<@Geno;$i++) {
			for (my $j=$i+1;$j<@Geno;$j++) {
				if($Geno[$i] ne $Geno[$j] && $Geno[$i] ne "./." && $Geno[$j] ne "./."){
					$diff{$indi[$j]}{$indi[$i]}++ ;
					$diff{$indi[$i]}{$indi[$j]}++ ;
				};
			}
		}
	}
}
close In;
open Out,">$fmtrix";
my @sample=@indi;
print Out join("\t","",@sample),"\n";
for (my $i=0;$i<@sample;$i++) {
	my @out;
	push @out,$sample[$i];
	for (my $j=0;$j<@sample;$j++) {
		if (!exists $diff{$sample[$i]}{$sample[$j]}) {
			push @out,0;
		}else{
			push @out,$diff{$sample[$i]}{$sample[$j]};
		}
	}
	print Out join("\t",@out),"\n";
}
close Out;
open Out,">$fOut";
print Out "#sampleID\tSNPnumber\tTransition\tTransvertion\tTs/Tv\tHeterozygosity Number\tHomozygosity Number\tAverage Depth\tMiss Number\tRef Number\n";
foreach my $sample (sort keys %vcfstat) {
	next if ($sample eq "pop");
	print Out join("\t",$sample,$vcfstat{$sample}{total},$vcfstat{$sample}{Ts},$vcfstat{$sample}{Tv},sprintf("%.2f",$vcfstat{$sample}{Ts}/$vcfstat{$sample}{Tv}),$vcfstat{$sample}{total}-$vcfstat{$sample}{homo},$vcfstat{$sample}{homo},sprintf("%.2f",$vcfstat{$sample}{dp}/$vcfstat{$sample}{total}),$vcfstat{$sample}{miss},$vcfstat{$sample}{ref}),"\n";
}
print Out join("\t","pop",$vcfstat{pop}{total},$vcfstat{pop}{Ts},$vcfstat{pop}{Tv},sprintf("%.2f",$vcfstat{pop}{Ts}/$vcfstat{pop}{Tv}),"--","--",sprintf("%.2f",$vcfstat{pop}{dp}/$vcfstat{pop}{total}/scalar keys @indi),"--","--"),"\n";
close Out;
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub USAGE {#
        my $usage=<<"USAGE";
Contact:        long.huang\@majorbio.com;
Script:			$Script
Description:
	fq thanslate to fa format
	eg:
	perl $Script -i -o -k -c

Usage:
  Options:
  -i	<file>	input file name
  -o	<file>	output file name
  -m	<file>	output metric file
  -h         Help

USAGE
        print $usage;
        exit;
}

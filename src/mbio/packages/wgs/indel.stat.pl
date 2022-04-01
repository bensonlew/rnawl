#!/usr/bin/env perl
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fIn,$fOut,$report,$fmtrix,$flen);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$fIn,
	"o:s"=>\$fOut,
	"m:s"=>\$fmtrix,
	"l:s"=>\$flen,
			) or &USAGE;
&USAGE unless ($fIn );
open In,$fIn;
if ($fIn=~/.gz$/) {
	close In;
	open In,"gunzip -c $fIn|";
}
open Flen,">$flen";
my @indi;
my %vcfstat;
my %diff;
my %len;
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
		my $sumlen=0;
		for (my $i=0;$i<@geno;$i++) {
			my ($gt,$ad,$dp,undef)=split(/\:/,$geno[$i]);
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
			#$geno{$gt}++ if($gt ne "./.");
			$sum+=$dp if($gt ne "./.");
			$vcfstat{$indi[$i]}{total}++;
			$vcfstat{$indi[$i]}{dp}+=$dp;
			$vcfstat{$indi[$i]}{homo}++ if($g1 eq $g2);
			$vcfstat{$indi[$i]}{lenth}+=length($ref)*2-(length($Ale[$g1])+length($Ale[$g2]));
			$sumlen=length($ref)*2-(length($Ale[$g1])+length($Ale[$g2]));
			print Flen $indi[$i],"\t",(length($Ale[$g1])+length($Ale[$g2]))-length($ref)*2,"\n" if (length($ref)*2-(length($Ale[$g1])+length($Ale[$g2])) !=0);
			if (length($ref)*2 -length($Ale[$g1])-length($Ale[$g2]) < 0) {
				$vcfstat{$indi[$i]}{insert}++;
			}else{
				$vcfstat{$indi[$i]}{delete}++;
			}
		}
		$vcfstat{pop}{dp}+=$sum if(scalar keys %geno != 1 ||scalar @indi == 1);
		$vcfstat{pop}{len}+=$sumlen if(scalar keys %geno != 1 ||scalar @indi == 1);
		 if(scalar keys %geno == 2){
			if ($sumlen < 0) {
				$vcfstat{pop}{insert}++;
			}else{
				$vcfstat{pop}{delete}++;
			}
		 };
		 $vcfstat{pop}{total}++ if(scalar keys %geno !=1|| scalar @indi == 1);
		for (my $i=0;$i<@Geno;$i++) {
			for (my $j=$i+1;$j<@Geno;$j++) {
				if($Geno[$i] ne $Geno[$j] && $Geno[$i] ne "./." && $Geno[$j] ne "./."){
					$diff{$indi[$j]}{$indi[$i]}++ ;
					$diff{$indi[$i]}{$indi[$j]}++ ;
				};
			}
		}
	}
#    print Dumper %vcfstat;
}
close In;
close Flen;
open Out,">$fmtrix";
my @sample=@indi;
print Out join("\t","",@sample),"\n";
for (my $i=0;$i<@sample;$i++) {
	my @out;push @out,$sample[$i];
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
print Out "Sample ID\tInsert Number\tDelete Number\tHeterozygosity Number\tHomozygosity Number\tAverage Length\tAverage Depth\tMiss Number\tRef Number\n";
foreach my $sample (sort keys %vcfstat) {
	$vcfstat{$sample}{lenth}||=0;
	$vcfstat{$sample}{insert}||=0;
	$vcfstat{$sample}{delete}|=0;
	next if($sample eq "pop");
	print Out join("\t",$sample,$vcfstat{$sample}{insert},$vcfstat{$sample}{delete},$vcfstat{$sample}{total}-$vcfstat{$sample}{homo},$vcfstat{$sample}{homo},sprintf("%.2f",$vcfstat{$sample}{lenth}/$vcfstat{$sample}{total}),sprintf("%.2f",$vcfstat{$sample}{dp}/$vcfstat{$sample}{total}),$vcfstat{$sample}{miss},$vcfstat{$sample}{ref}),"\n";
}
print Out join("\t","pop",$vcfstat{pop}{insert},$vcfstat{pop}{delete},"--","--",sprintf("%.2f",$vcfstat{pop}{lenth}/$vcfstat{pop}{total}),sprintf("%.2f",$vcfstat{pop}{dp}/$vcfstat{pop}{total}),"--","--"),"\n";
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
  -o	<file>	output stat name
  -m	<file>	output matrix file
  -l	<file>	output len file
  -h         Help

USAGE
        print $usage;
        exit;
}

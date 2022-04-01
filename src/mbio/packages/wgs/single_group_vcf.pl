#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($group,$vcf,$funtype,$efftype,$fOut,$len,$miss,$maf,$dep,$position);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"gro:s"=>\$group,
	"vcf:s"=>\$vcf,
	"out:s"=>\$fOut,
	"funtype:s"=>\$funtype,
	"efftype:s"=>\$efftype,
	"len:s"=>\$len,
	"dep:s"=>\$dep,
	"miss:s"=>\$miss,
	"maf:s"=>\$maf,
	"pos:s"=>\$position,

		) or &USAGE;
&USAGE unless ( $vcf && $fOut);
$group||="all";
$maf||="0,1.2";
$miss||="0";
$dep||="-1,100000000";
$len||="-1,100000000";
$efftype||="all";
$funtype||="all";
my ($chr,$start,$end);
if (defined $position) {
	($chr,$start,$end)=split(/,/,$position);
}
my ($maf1,$maf2)=split(/,/,$maf);
my ($dep1,$dep2)=split(/,/,$dep);
my ($len1,$len2)=split(/,/,$len);
my (%efftype,%funtype,%sample);
if ($efftype ne "all"){
	my @eff=split(/,/,$efftype);
	foreach my $eff(@eff){
		$efftype{$eff}=1;
	}
}
if ($funtype ne "all"){
	my @funtype=split(/,/,$funtype);
	foreach my $funtype(@funtype){
		$funtype{$funtype}=1;
	}
}
if($group ne "all"){
	my @sam=split(/,/,$group);
	foreach my $sam(@sam){
		$sample{$sam}=1;
	}
}

open In,$vcf;
mkdir $fOut if (!-d $fOut);
open Out,">$fOut/pop.variant";
open Pos,">$fOut/pos.variant";

if ($vcf =~/.gz/) {
	close In;
	open In,"zcat $vcf|";
}
my @Indi;
my %Annstat;
my %EFFstat;
my %Distri;
my %Diff;
my %Length;
my %SampleStat;
while (<In>){
	chomp;
	next if ($_ eq ""||/^$/ || /^##/);
	if (/^#/){
		my ($chrom,$pos,$id,$ref,$alt,$qual,$filter,$info,$format,@indi)=split(/\t/,$_);
		@Indi=@indi;
		if ($group eq "all") {
			foreach my $indi (@indi) {
				$sample{$indi}=1;
			}
		}
		print Out join("\t","#chrom","pos","ref","alt","average","miss","maf","frequence","ANN"),"\n";
		print Pos join("\t","#chrom","pos","ref","alt","average","miss","maf","frequence","ANN"),"\n";
	}else{
		my ($chrom,$pos,$id,$ref,$alt,$qual,$filter,$info,$format,@geno)=split(/\t/,$_);
		next if ($alt =~ /\*/);
		my @ale=split(/,/,join(",",$ref,$alt));
		my %len;
		my $minlen=1000000;
		my $maxlen=0;
		foreach my $ale(@ale){
			$len{length($ale)}=1;
		}
		my $vtype="SNP";
		$vtype="INDEL"if(scalar keys %len>1);
		next if (scalar keys %len > 1 && $len1 eq $len2);
		my %ANNdetail;
		my %FUNdetail;
		my %anninfo;
		my %annout;
		my $efftest=0;
		my $funtest=0;
		my @format=split(/:/,$format);
		my $tdp=0;
		my %gstat;
		my %astat;
		my $asum=0;;
		my %gt;
		my %annstat;
		my %samplestat;
		my %effstat;
		my %lenth;
		if($info=~/ANN=([^;]*)/){
			my $anns=$1;
			my @ann=split(/\,/,$anns);
			foreach my $ann(@ann){
				my @anns=split(/\|/,$ann);
				my @ann1;
				$anns[1]||="--";
				push @ann1,$anns[1];
				if ($anns[1] =~ /\&/) {
					@ann1=split(/\&/,$anns[1]);
				}
				$anns[0]||="--";
				$anns[2]||="--";
				$anns[3]||="--";
				$anns[4]||="--";
				foreach my $ann1 (@ann1) {
					push @{$anninfo{$anns[0]}},join("|",$anns[0],$ann1,$anns[2],$anns[3],$anns[4]);
					$ANNdetail{$anns[0]}{$ann1}++;
				}
				$FUNdetail{$anns[0]}{$anns[2]}++;
			}
		}

		for (my $i=0;$i<@Indi;$i++){
			next if (!exists $sample{$Indi[$i]});
			next if ($geno[$i] eq ".");
			my @info=split(/\:/,$geno[$i]);
			for(my $j=0;$j<@format;$j++){
				if($format[$j] eq "GT"){
					my ($g1,$g2)=split(/\//,$info[$j]);
					my $gt=$info[$j];
					if($g1 ne "."){
						$gt=join("/",$ale[$g1],$ale[$g2]);
						if($g1 !=0){
							if (!defined $anninfo{$ale[$g1]}) {
								print Dumper %anninfo;
								print $ale[$g1];die;
							}
							my $ann=join(";",@{$anninfo{$ale[$g1]}});
							$annout{$ann}=1;
							foreach my $ann(sort keys %{$ANNdetail{$ale[$g1]}}){
								if (exists $efftype{$ann}){
									$efftest=1;
								}
								$annstat{$ann}++;
								$samplestat{$ann}{$Indi[$i]}+=$ANNdetail{$ale[$g1]}{$ann};
							}
							foreach my $ann(sort keys %{$FUNdetail{$ale[$g1]}}){
								if (exists $funtype{$ann}){
								   $funtest=1;
								}
								$effstat{$ann}++;
							    $samplestat{$ann}{$Indi[$i]}+=$FUNdetail{$ale[$g1]}{$ann};
							}
						}
						if($g2 != 0 && $g1 ne $g2){
							if (!defined $anninfo{$ale[$g2]}) {
								print Dumper %anninfo;
								print $ale[$g2];die;
							}
							my $ann=join(";",@{$anninfo{$ale[$g2]}});
							$annout{$ann}=1;
							foreach my $ann(sort keys %{$ANNdetail{$ale[$g2]}}){
								    if (exists $efftype{$ann}){
										$efftest=1;
									}
								$annstat{$ann}++;
								$samplestat{$ann}{$Indi[$i]}+=$ANNdetail{$ale[$g2]}{$ann};
							}
							foreach my $ann(sort keys %{$FUNdetail{$ale[$g2]}}){
								    if (exists $funtype{$ann}){
									    $funtest=1;
									}
								$effstat{$ann}++;
							    $samplestat{$ann}{$Indi[$i]}+=$FUNdetail{$ale[$g2]}{$ann};
							}
						}
						$gt{$Indi[$i]}=$gt;
						$lenth{length($ale[$g2])-length($ref)}++;
						$lenth{length($ale[$g1])-length($ref)}++;
						$minlen=length($ale[$g1]) if(length($ale[$g1]) < $minlen);
						$maxlen=length($ale[$g1]) if(length($ale[$g1]) > $maxlen);
						$minlen=length($ale[$g2]) if(length($ale[$g2]) < $minlen);
						$maxlen=length($ale[$g2]) if(length($ale[$g2]) > $maxlen);
						$asum+=2;
					}
					$gstat{$gt}++;
					$astat{$g1}++ if($g2 ne ".");
					$astat{$g2}++ if($g1 ne ".");
				}
				if($format[$j] eq "AD"){
					my @ad=split(/,/,$info[$j]);
					if ($info[$j] eq ".") {
						print $geno[$i];die;
					}
					foreach my $ad(@ad){
						$tdp+=$ad;
					}
				}
			}
		}
		next if($vtype eq "INDEL" && $minlen > $len2 && $maxlen < $len1);
		next if($asum == 0);
		$gstat{"./."}||=0;
		my $missfilter=sprintf("%.4f", $gstat{"./."}/scalar keys %sample);
		my $average=sprintf("%.4f",$tdp/scalar keys %sample);
		my @mafale=sort{$b<=>$a} values %astat;
		$mafale[1]||=0;
		my $maf=sprintf("%.4f", $mafale[1]/$asum);

		next if ($efftest == 0 && $efftype ne "all");
		next if ($funtest == 0 && $funtype ne "all");
		next if ($average < $dep1 && $average > $dep2);
		next if ($maf < $maf1 && $maf> $maf2);
		next if ($missfilter < $miss);
		foreach my $ann (sort keys %annstat) {
			$Annstat{$ann}++;
		}
		foreach my $ann (sort keys %effstat) {
			$EFFstat{$ann}++;
		}
		foreach my $sample (sort keys %samplestat) {
			foreach my $anns (sort keys %{$samplestat{$sample}}) {
				$SampleStat{$sample}{$anns}+=$samplestat{$sample}{$anns};
			}
		}
		foreach my $len (sort keys %lenth) {
			$Length{$len}+=$lenth{$len};
		}

		my @fre;
		foreach my $gt(sort keys %gstat){
			my $t=sprintf("%.4f",$gstat{$gt}/(scalar keys %sample)*100);
			push @fre,"$gt(".$t."%)";
		}
		my @sample=keys %gt;
		for (my $i=0;$i<@sample;$i++) {
			for (my $j=$i+1;$j<@sample;$j++) {
				$Diff{$sample[$i]}{$sample[$j]}++ if($gt{$sample[$i]} ne $gt{$sample[$j]});
				$Diff{$sample[$j]}{$sample[$i]}++ if($gt{$sample[$i]} ne $gt{$sample[$j]});
			}
		}
		print Out join("\t",$chrom,$pos,$ref,$alt,$average,$missfilter,$maf,join(";",@fre),join(";",keys %annout)),"\n";
		$Distri{$chrom}{$pos}=1;
		next if (defined $position && ($chr ne $chrom || ($pos < $start || $pos > $end)));
		print Pos join("\t",$chrom,$pos,$ref,$alt,$average,$missfilter,$maf,join(";",@fre),join(";",keys %annout)),"\n";

	}
}
close In;
close Out;
open Out,">$fOut/Ann.stat";
foreach my $ann (sort keys %Annstat) {
	next if (!exists $efftype{$ann} && $efftype ne "all");
	print Out join("\t",$ann,$Annstat{$ann}),"\n";
}
close Out;
open Out,">$fOut/Eff.stat";
foreach my $ann (sort keys %EFFstat) {
	next if (!exists $funtype{$ann} && $funtype ne "all");
	print Out join("\t",$ann,$EFFstat{$ann}),"\n";
}
close Out;
open Out,">$fOut/win.stat";
foreach my $chr (sort keys %Distri) {
	foreach my $win (sort {$a<=>$b} keys %{$Distri{$chr}}) {
		print Out join("\t",$chr,$win,$Distri{$chr}{$win}),"\n";
	}
}
close Out;
open Out,">$fOut/sample.eff.stat";
print Out join("\t","SampleID",sort keys %SampleStat),"\n";
foreach my $sample(keys %sample){
	my @out;
	push @out,$sample;
	foreach my $ann(sort keys %SampleStat){
		$SampleStat{$ann}{$sample}||=0;
		push @out,$SampleStat{$ann}{$sample};
	}
	print Out join("\t",@out),"\n";
}
close Out;
open Out,">$fOut/diff.matrix";
for (my $i=0;$i<@Indi;$i++) {
	my @out;
	push @out,$Indi[$i];
	for (my $j=0;$j<@Indi;$j++) {
		$Diff{$Indi[$i]}{$Indi[$j]}||=0;
		push @out,$Diff{$Indi[$i]}{$Indi[$j]};
	}
	print Out join("\t",@out),"\n";
}
close Out;
open Out,">$fOut/Length.txt";
foreach my $l (sort keys %Length) {
	print Out $l,"\t",$Length{$l},"\n";
}
close Out;

#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub USAGE {#
        my $usage=<<"USAGE";
Contact:        minghao.zhang\@majorbio.com;
Script:			$Script
Description:
	fq thanslate to fa format
	eg:
	perl $Script -i -o -k -c

Usage:
  Options:
	"gro:s"=>\$group,
	"vcf:s"=>\$vcf,
	"out:s"=>\$out,
	"funtype:s"=>\$funtype,
	"efftype:s"=>\$efftype,
	"len:s"=>\$len,
	"dep:s"=>\$dep,
	"miss:s"=>\$miss,
	"maf:s"=>\$maf,
  -h         Help

USAGE
        print $usage;
        exit;
}

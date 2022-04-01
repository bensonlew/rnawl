#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
#my ($sample,$sample2,$vcf,$funtype,$efftype,$fOut,$len1,$len2,$miss1,$miss2,$dep1,$dep2,$equal,$position);
my ($sample,$vcf,$funtype,$efftype,$fOut,$len1,$miss1,$miss2,$dep1,$equal,$position);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"sample:s"=>\$sample, 
	"equal:s"=>\$equal,
	"vcf:s"=>\$vcf,
	"out:s"=>\$fOut,
	"funtype:s"=>\$funtype,
	"efftype:s"=>\$efftype,
	"len1:s"=>\$len1,
	"dep1:s"=>\$dep1,
	"pos:s"=>\$position,

		) or &USAGE;
&USAGE unless ( $vcf && $fOut && $sample);
$equal||="0";
$dep1||="-1,100000000";
$len1||="-1,100000000";
$efftype||="all";
$funtype||="all";
my ($chr,$start,$end);
if (defined $position) {
	($chr,$start,$end)=split(/,/,$position);
}

my ($dep1_1,$dep1_2)=split(/,/,$dep1);
my ($len1_1,$len1_2)=split(/,/,$len1);
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

open In,$vcf;
mkdir $fOut if (!-d $fOut);
open Out,">$fOut/pop.variant";
if ($vcf =~/.gz/) {
	close In;
	open In,"zcat $vcf|";
}
open Pos,">$fOut/pos.variant";

my @Indi;
my %Annstat;
my %EFFstat;
my %Distri;
my %Length;
my %SampleStat;
while (<In>){
	chomp;
	next if ($_ eq ""||/^$/ || /^##/);
	if (/^#/){
		my ($chrom,$pos,$id,$ref,$alt,$qual,$filter,$info,$format,@indi)=split(/\t/,$_);
		@Indi=@indi;
		#print Out join("\t","#chrom","pos","ref","alt","$sample\_gt","$sample2\_gt","$sample\_ad","$sample2\_ad","$sample\_dp","$sample2\_dp","ANN"),"\n";
		#print Pos join("\t","#chrom","pos","ref","alt","$sample\_gt","$sample2\_gt","$sample\_ad","$sample2\_ad","$sample\_dp","$sample2\_dp","ANN"),"\n";
        print Out join("\t","#chrom","pos","ref","alt","$sample\_gt","$sample\_ad","$sample\_dp","ANN"),"\n";
        print Pos join("\t","#chrom","pos","ref","alt","$sample\_gt","$sample\_ad","$sample\_dp","ANN"),"\n";
	}else{
		my ($chrom,$pos,$id,$ref,$alt,$qual,$filter,$info,$format,@geno)=split(/\t/,$_);
        my $reftype="$ref/$ref";
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
		next if (scalar keys %len > 1  && $len1_1 eq 1);
		my @format=split(/:/,$format);
		my %tdp;
		my %gstat;
		my %astat;
		my %asum;
		my %gt;
		my %ANNdetail;
		my %FUNdetail;
		my %anninfo;
		my %annout;
		my $efftest=0;
		my $funtest=0;
		my %dep;
		my %minlen;
		my %aledep;
		my %maxlen;
		$minlen{$sample}=1000000000;
		$maxlen{$sample}=-1;
		#$minlen{$sample2}=1000000000;
		#$maxlen{$sample2}=-1;
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
		my %annstat;
		my %samplestat;
		my %effstat;
		my %lenth;
		for (my $i=0;$i<@Indi;$i++){
			next if ($Indi[$i] ne $sample );##delete && $Indi[$i] ne $sample2
            next if ($geno[$i] eq ".");
			my @info=split(/\:/,$geno[$i]);
			for(my $j=0;$j<@format;$j++){
				if($format[$j] eq "GT"){
					my ($g1,$g2)=split(/\//,$info[$j]);
					my $gt=$info[$j];
                    #print $gt,"\n";
					if($g1 ne "."){
						$gt=join("/",$ale[$g1],$ale[$g2]);
						if($g1 !=0){
							if (!defined $anninfo{$ale[$g1]}) {
								print Dumper %anninfo;
								#print $ale[$g1];die;
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
						$gt=join("/",$ale[$g1],$ale[$g2]);
						$lenth{length($ale[$g2])-length($ref)}++;
						$lenth{length($ale[$g1])-length($ref)}++;
						$minlen{$Indi[$i]}=length($ale[$g1]) if(length($ale[$g1]) < $minlen{$Indi[$i]});
						$maxlen{$Indi[$i]}=length($ale[$g1]) if(length($ale[$g1]) > $maxlen{$Indi[$i]});
						$minlen{$Indi[$i]}=length($ale[$g2]) if(length($ale[$g2]) < $minlen{$Indi[$i]});
						$maxlen{$Indi[$i]}=length($ale[$g2]) if(length($ale[$g2]) > $maxlen{$Indi[$i]});
					}
					$gt{$Indi[$i]}=$gt;
				}
				if($format[$j] eq "AD"){
                    if(!defined $info[$j]){print Dumper @info;die;}
					$aledep{$Indi[$i]}=$info[$j];
					my @ad=split(/,/,$info[$j]);
					foreach my $ad(@ad){
						$tdp{$Indi[$i]}+=$ad;
					} 
				}
			}
		}
		$gt{$sample}||="./.";
		#$gt{$sample2}||="./.";
		next if($vtype eq "INDEL" && $minlen{$sample} > $len1_2 && $maxlen{$sample} < $len1_1);
		#next if($vtype eq "INDEL" && $minlen{$sample2} > $len2_2 && $maxlen{$sample2} < $len2_1);
        next if ($gt{$sample} eq "./.");## || $gt{$sample2} eq "./.");
		next if ($efftest == 0 && $efftype ne "all");
		next if ($funtest == 0 && $funtype ne "all");
		next if ($tdp{$sample} < $dep1_1 && $tdp{$sample} > $dep1_2);
		#next if ($tdp{$sample2} < $dep2_1 && $tdp{$sample2} > $dep2_2);
        next if ($gt{$sample} eq "$ref/$ref" and $equal !=1);
		foreach my $ann (sort keys %annstat) {
			$Annstat{$ann}+=$annstat{$ann};
		}
		foreach my $ann (sort keys %effstat) {
			$EFFstat{$ann}+=$effstat{$ann};
		}
		foreach my $sample (sort keys %samplestat) {
			foreach my $anns (sort keys %{$samplestat{$sample}}) {
				$SampleStat{$sample}{$anns}+=$samplestat{$sample}{$anns};
			}
		}
		foreach my $len (sort keys %lenth) {
			$Length{$len}+=$lenth{$len};
		}
		print Out join("\t",$chrom,$pos,$ref,$alt,$gt{$sample},$aledep{$sample},$tdp{$sample},join(";",keys %annout)),"\n";
        $Distri{$chrom}{$pos}=1;
		next if (defined $position && ($chr ne $chrom || ($pos < $start || $pos > $end)));
        print Pos join("\t",$chrom,$pos,$ref,$alt,$gt{$sample},$aledep{$sample},$tdp{$sample},join(";",keys %annout)),"\n";
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
	foreach my $win (sort {$a<=>$b}keys %{$Distri{$chr}}) {
		print Out join("\t",$chr,$win,$Distri{$chr}{$win}),"\n";
	}
}
close Out;
open Out,">$fOut/Length.txt";
foreach my $l (sort keys %Length) {
	print Out $l,"\t",$Length{$l},"\n";
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
	"sample:s"=>\$sample, 
	"equal:s"=>\$equal,"equal or diff"
	"vcf:s"=>\$vcf,
	"out:s"=>\$fOut,
	"funtype:s"=>\$funtype,
	"efftype:s"=>\$efftype,
	"len1:s"=>\$len1,
	"len2:s"=>\$len2,
	"dep1:s"=>\$dep1,
	"dep2:s"=>\$dep2,

  -h         Help

USAGE
        print $usage;
        exit;
}

#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($group1,$group2,$g1name,$g2name,$vcf,$funtype,$efftype,$fOut,$len1,$len2,$miss1,$miss2,$dep1,$dep2,$maf1,$maf2,$position);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"gro1:s"=>\$group1,
	"gro2:s"=>\$group2,
	"g1name:s"=>\$g1name,
	"g2name:s"=>\$g2name,
	"vcf:s"=>\$vcf,
	"out:s"=>\$fOut,
	"funtype:s"=>\$funtype,
	"efftype:s"=>\$efftype,
	"len1:s"=>\$len1,
	"len2:s"=>\$len1,
	"dep1:s"=>\$dep1,
	"dep2:s"=>\$dep2,
	"miss1:s"=>\$miss1,
	"miss2:s"=>\$miss2,
	"maf1:s"=>\$maf1,
	"maf2:s"=>\$maf2,
	"pos:s"=>\$position,
		) or &USAGE;
&USAGE unless ( $vcf && $fOut and $group1 and $group2);
$maf1||="0,1.2";
$maf2||="0,1.2";
$miss1||="0";
$miss2||="0";
$dep1||="-1,100000000";
$dep2||="-1,100000000";
$len1||="-1,100000000";
$len2||="-1,100000000";
my ($chr,$start,$end);
if (defined $position) {
	($chr,$start,$end)=split(/,/,$position);
}
$efftype||="all";
$funtype||="all";
my ($maf1_1,$maf1_2)=split(/,/,$maf1);
my ($dep1_1,$dep1_2)=split(/,/,$dep1);
my ($maf2_1,$maf2_2)=split(/,/,$maf2);
my ($dep2_1,$dep2_2)=split(/,/,$dep2);
my ($len1_1,$len1_2)=split(/,/,$len1);
my ($len2_1,$len2_2)=split(/,/,$len2);
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
my @sam1=split(/,/,$group1);
foreach my $sam(@sam1){
	$sample{$sam}=1;
}
my @sam2=split(/,/,$group2);
foreach my $sam(@sam2){
	$sample{$sam}=2;
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
		print Out join("\t","#chrom","pos","ref","ale","$g1name\_average_dep","$g1name\_miss_ratio","$g1name\_maf","$g1name\_frequence","$g2name\_average_dep","$g2name\_miss_ratio","$g2name\_maf","$g2name\_frequence","ANN"),"\n";
		print Pos join("\t","#chrom","pos","ref","ale","$g1name\_average_dep","$g1name\_miss_ratio","$g1name\_maf","$g1name\_frequence","$g2name\_average_dep","$g2name\_miss_ratio","$g2name\_maf","$g2name\_frequence","ANN"),"\n";
	}else{
		my ($chrom,$pos,$id,$ref,$alt,$qual,$filter,$info,$format,@geno)=split(/\t/,$_);
		next if ($alt =~ /\*/);
		my @ale=split(/,/,join(",",$ref,$alt));
		my %len;
		my %minlen;
		my %maxlen;
		foreach my $ale(@ale){
			$len{length($ale)}=1;
		}
		my $vtype="SNP";
		$vtype="INDEL"if(scalar keys %len>1);
		next if (scalar keys %len > 1 && $len1 eq $len2 && $len1_1 eq 1);
		my %ANNdetail;
		my %FUNdetail;
		my %anninfo;
		my %annout;
		my $efftest=0;
		my $funtest=0;
		my @format=split(/:/,$format);
		my %tdp;
		my %gstat;
		my %astat;
		my %asum;
		my %gt;
		$minlen{1}=1000000000;
		$maxlen{1}=0;
		$minlen{2}=1000000000;
		$maxlen{2}=0;
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

			my $GID=$sample{$Indi[$i]};
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
						$lenth{length($ale[$g2])-length($ref)}++;
						$lenth{length($ale[$g1])-length($ref)}++;
						$minlen{$GID}=length($ale[$g1]) if(length($ale[$g1]) < $minlen{$GID});
						$maxlen{$GID}=length($ale[$g1]) if(length($ale[$g1]) > $maxlen{$GID});
						$minlen{$GID}=length($ale[$g2]) if(length($ale[$g2]) < $minlen{$GID});
						$maxlen{$GID}=length($ale[$g2]) if(length($ale[$g2]) > $maxlen{$GID});
						$asum{$GID}+=2;		
					}
					$gstat{$gt}{$GID}++;	
					$astat{$GID}{$g1}++ if($g2 ne ".");
					$astat{$GID}{$g2}++ if($g1 ne ".");
				}
				if($format[$j] eq "AD"){
					my @ad=split(/,/,$info[$j]);
					foreach my $ad(@ad){
						$tdp{$GID}+=$ad;
					}
				}
			}
		}
		next if($vtype eq "INDEL" && $minlen{1} > $len1_2 && $maxlen{1} < $len1_1);
		next if($vtype eq "INDEL" && $minlen{2} > $len2_2 && $maxlen{2} < $len2_1);
		$gstat{"./."}{1}||=0;
		$gstat{"./."}{2}||=0;
		my $missfilter1=$gstat{"./."}{1}/scalar @sam1;
		my $average1=$tdp{1}/scalar  @sam1;
		my $missfilter2=$gstat{"./."}{2}/scalar  @sam2;
		my $average2=$tdp{2}/scalar  @sam2;
		next if ($missfilter1 < $miss1 || !exists $asum{1});
		next if ($missfilter2 < $miss2 || !exists $asum{2});
		my @mafale1=sort{$b<=>$a} values %{$astat{1}};
		$mafale1[1]||=0;
		my $maf1=$mafale1[1]/$asum{1};
		my @mafale2=sort{$b<=>$a} values %{$astat{2}};
		$mafale2[1]||=0;
		my $maf2=$mafale2[1]/$asum{2};

		next if ($efftest == 0 && $efftype ne "all");
		next if ($funtest == 0 && $funtype ne "all");
		next if ($average1 < $dep1_1 && $average1 > $dep1_2);
		next if ($average2 < $dep2_1 && $average2 > $dep2_2);

		next if ($maf1 < $maf1_1 && $maf1> $maf1_2);
		next if ($maf2 < $maf2_1 && $maf2> $maf2_2);
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
		my (@fre1,@fre2);
		foreach my $gt(sort keys %gstat){
			$gstat{$gt}{1}||=0;
			$gstat{$gt}{2}||=0;
			my $t1=sprintf("%.4f",$gstat{$gt}{1}/(scalar @sam1)*100);
			my $t2=sprintf("%.4f",$gstat{$gt}{2}/(scalar @sam2)*100);
			push @fre1,"$gt(".$t1."%)";
			push @fre2,"$gt(".$t2."%)";
		}
		my @sample=keys %gt;
		for (my $i=0;$i<@sample;$i++) {
			for (my $j=$i+1;$j<@sample;$j++) {
				$Diff{$sample[$i]}{$sample[$j]}++ if($gt{$sample[$i]} ne $gt{$sample[$j]});
				$Diff{$sample[$j]}{$sample[$i]}++ if($gt{$sample[$i]} ne $gt{$sample[$j]});
			}
		}
		print Out join("\t",$chrom,$pos,$ref,$alt,$average1,$missfilter1,$maf1,join(";",@fre1),$average2,$missfilter2,$maf2,join(";",@fre2),join(";",keys %annout)),"\n";
		$Distri{$chrom}{$pos}=1;
		next if (defined $position && ($chr ne $chrom || ($pos < $start || $pos > $end)));
		print Pos join("\t",$chrom,$pos,$ref,$alt,$average1,$missfilter1,$maf1,join(";",@fre1),$average2,$missfilter2,$maf2,join(";",@fre2),join(";",keys %annout)),"\n";

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
	"gro1:s"=>\$group1,
	"gro2:s"=>\$group2,
	"g1name:s"=>\$g1name,
	"g2name:s"=>\$g2name,
	"vcf:s"=>\$vcf,
	"out:s"=>\$fOut,
	"funtype:s"=>\$funtype,
	"efftype:s"=>\$efftype,
	"len1:s"=>\$len1,
	"len2:s"=>\$len1,
	"dep1:s"=>\$dep1,
	"dep2:s"=>\$dep2,
	"miss1:s"=>\$miss1,
	"miss2:s"=>\$miss2,
	"maf1:s"=>\$maf1,
	"maf2:s"=>\$maf2,
	"pos:s"=>\$position,

  -h         Help

USAGE
        print $usage;
        exit;
}

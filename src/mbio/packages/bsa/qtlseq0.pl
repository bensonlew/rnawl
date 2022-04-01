#!/usr/bin/env perl 
use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Data::Dumper;
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fIn,$fOut,$PID,$BID,$Pdep,$Bdep,$popt);
GetOptions(
				"help|?" =>\&USAGE,
				"vcf:s"=>\$fIn,
				"out:s"=>\$fOut,
				"pid:s"=>\$PID,
				"bid:s"=>\$BID,
				"pdep:s"=>\$Pdep,
				"Bdep:s"=>\$Bdep,
				"popt:s"=>\$popt,
				) or &USAGE;
&USAGE unless ($fIn and $fOut and $BID);
$Pdep||=10;
$Bdep||=10;
$PID||="";
$popt||="F2";
my %Indi;
my ($P1,$P2);
if ($PID ne "") {
	($P1,$P2)=split(/\,/,$PID);
	$P2||=$P1;
	$Indi{$P1}="P1";
	$Indi{$P2}="P2";
}
my ($B1,$B2)=split(/\,/,$BID);
$Indi{$B1}="B1";
$Indi{$B2}="B2";
if ($fIn =~ /gz$/) {
	open In,"gunzip -dc $fIn|";
}else{
	open In,$fIn;
}
open Out,">$fOut.index";
#print Out join("\t","#chr","pos","ref","alt","anno","GT","AD","GT","AD","B1GT","B1AD","B2GT","B2AD","INDEX1","INDEX2","DELTA"),"\n";
open Variant,">$fOut.variant";
my @indi;
my %stat;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/ ||/^##/);
	my ($chr,$pos,$ids,$ref,$alt,$qual,$filter,$info,$format,@geno)=split(/\s+/,$_);
	if (/^#/) {
		push @indi,@geno;
		my @out;
		foreach my $indi (@geno) {
			next if (!exists $Indi{$indi});
			push @out,$indi."-GT";
			push @out,$indi."-AD";
		}
		print Out join("\t","#chr","pos","type","ref",@out,"index","ANNOTATION"),"\n";
		print Variant join("\t","#chr","type","pos","ref",@out,"ANNOTATION","HIGH","MODERATE","LOW","MODIFIER"),"\n";

	}else{
		my %info;
		my @alle=split(",",join(",",$ref,$alt));
		my %len;
		for (my $i=0;$i<@alle;$i++) {
			$len{$alle[$i]}=1;
		}
		my $type="SNP";
		if (scalar keys %len > 1) {
			$type="INDEL";
		}
		my @format=split(/:/,$format);
		my %ginfo;
		my @outvariant;
		for (my $i=0;$i<@indi;$i++) {
			next if (!exists $Indi{$indi[$i]});
			my @gt=split(/\//,$geno[$i]);
			my $id=$Indi{$indi[$i]};
			my @info=split(/:/,$geno[$i]);
			for (my $j=0;$j<@info;$j++) {
				$info{$id}{gt}=$info[$j] if ($format[$j] eq "GT");
				$info{$id}{ad}=$info[$j] if ($format[$j] eq "AD");
				$info{$id}{dp}=$info[$j] if ($format[$j] eq "DP");
			}
			if ($info{$id}{gt} eq "./.") {
				$ginfo{$indi[$i]}{gt}="--";
				$ginfo{$indi[$i]}{ad}="0";
			}else{
				my ($g1,$g2)=split(/\//,$info{$id}{gt});
				my @ad=split(/\,/,$info{$id}{ad});
				$ginfo{$indi[$i]}{gt}=join("/",$alle[$g1],$alle[$g2]);
				if ($g1 eq $g2) {
					$ginfo{$indi[$i]}{ad}=$ad[$g1];
				}else{
					$ginfo{$indi[$i]}{ad}=join(",",$ad[$g1],$ad[$g2]);
				}
			}
			push @outvariant,$ginfo{$indi[$i]}{gt};
			push @outvariant,$ginfo{$indi[$i]}{ad};
		}
		my @out;
					my %ann;

		if($info=~/ANN=([^\;]*)/g){
			my @ann=split(/\,/,$1);
			for (my $i=0;$i<@ann;$i++) {
				my @str=split(/\|/,$ann[$i]);
				my $ann=join("|",$str[1],$str[2],$str[3],$str[4]);
				$ann{$str[2]}++;
				push @out,$ann;
			}
		}
		$ann{HIGH}||=0;
		$ann{MODERATE}||=0;
		$ann{LOW}||=0;
		$ann{MODIFIER}||=0;

		if (scalar @out ==0) {
			print Variant join("\t",$chr,$pos,$type,"$ref"."$ref",@outvariant,"--",$ann{HIGH},$ann{MODERATE},$ann{LOW},$ann{MODIFIER}),"\n";
		}else{
			print Variant join("\t",$chr,$pos,$type,"$ref"."$ref",@outvariant,join(";",@out),$ann{HIGH},$ann{MODERATE},$ann{LOW},$ann{MODIFIER}),"\n";
		}

		next if ($info{B1}{gt}  eq "./." || $info{B2}{gt} eq "./." ||$info{B1}{dp} < $Bdep || $info{B2}{dp} < $Bdep);
		my @b1=split(/\/|\|/,$info{B1}{gt});
		my @b2=split(/\/|\|/,$info{B2}{gt});
		my @ad1=split(/\,/,$info{B1}{ad});
		my @ad2=split(/\,/,$info{B2}{ad});
		my $sum1=$ad1[$b1[0]];
		my $sum2=$ad2[$b2[0]];
		$sum1+=$ad1[$b1[1]] if($b1[1] ne $b1[0]);
		$sum2+=$ad2[$b2[1]] if($b2[1] ne $b2[0]);
		next if ($sum1 == 0 || $sum2 == 0);
		my %stat;
		$stat{$b1[0]}++;
		$stat{$b1[1]}++;
		$stat{$b2[0]}++;
		$stat{$b2[1]}++;
		my @geno=sort {$stat{$a}<=>$stat{$b}} keys %stat;
		my $index1;
		my $index2;
		my $delta;
		if ($popt eq "F2") {
			next if (scalar @geno!=2);
			if ($PID eq "") {
				$index1=$ad1[$geno[0]]/$sum1;
				$index2=$ad2[$geno[0]]/$sum2;
				$delta=abs($index1-$index2);
			}else{
				my ($p1,$p2)=split(/\/|\|/,$info{P1}{gt});
				if (!exists $info{P2}) {
					if ($p1 eq $geno[0]) {
						$info{P2}{gt}="$geno[0]q/$geno[0]";
						$info{P2}{dp}=10;
					}else{
						$info{P2}{gt}="$geno[1]/$geno[1]";
						$info{P2}{dp}=10;
					}
				}
				my ($p3,$p4)=split(/\/|\|/,$info{P2}{gt});
				next if ($p1 ne $p2 || $p3 ne $p4 || $p1 eq "." || $p2 eq ".");
				next if ($info{P1}{gt} eq "./." || $info{P2}{gt} eq "./.");
				next if	($info{P2}{dp}< $Pdep || $info{P1}{dp} < $Bdep);
				next if ($p1 eq $p3);
				next if (!exists $stat{$p1} || !exists $stat{$p3});
				$index1=$ad1[$p1]/$sum1;
				$index2=$ad2[$p1]/$sum2;
				$delta=$index1-$index2;
			}
		}else{
			next if (scalar @geno!=2);
			my ($p1,$p2)=split(/\/|\|/,$info{P1}{gt});
			my ($p3,$p4)=split(/\/|\|/,$info{P2}{gt});
			next if ($p1 eq $p2 || $p3 ne $p4);
			next if ($p1 ne $p2 || $p3 eq $p4);
			next if	($info{P2}{dp}< $Pdep || $info{P1}{dp} < $Bdep);
			if ($p3 eq $p1 && $p1 eq $p2) {#nnxnp
				$index1=$ad1[$p1]/$sum1;
				$index2=$ad2[$p1]/$sum2;
				$delta=$index1-$index2;
			}elsif ($p3 eq $4 && $p3 eq $p1) {#lmxll
				$index1=$ad1[$p2]/$sum1;
				$index2=$ad2[$p2]/$sum2;
				$delta=$index1-$index2;
			}else{
				next;
			}
		}
		$info{B1}{gt}||="--";
		$info{B1}{ad}||="--";
		$info{B2}{gt}||="--";
		$info{B2}{ad}||="--";
		$info{P1}{gt}||="--";
		$info{P1}{ad}||="--";
		$info{P2}{gt}||="--";
		$info{P2}{ad}||="--";
		my $lchr=$chr;
		if ($lchr !~ /chr/) {
			$lchr = "scaffords";
		}
		$stat{$lchr}{$type}++;
		if (scalar @out ==0) {
			print Out join("\t",$chr,$pos,$type,"$ref/"."$ref",@outvariant,$index1,$index2,$delta,"--",),"\n";
		}else{
			print Out join("\t",$chr,$pos,$type,"$ref/"."$ref",@outvariant,$index1,$index2,$delta,join(";",@out)),"\n";
		}
	}
}
close In;
close Out;
open Out,">$fOut.stat";
print Out "#chr\tsnp\tindel\n";
foreach my $chr (sort keys %stat) {
	print Out join("\t",$chr,$stat{$chr}{SNP},$stat{$chr}{InDel}),"\n";
}
close Out;
#######################################################################################
print "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

sub USAGE {#
	my $usage=<<"USAGE";
Description: 
Version:  $Script
Contact: long.huang

Usage:
  Options:
	-vcf	<file>	input file 
	-out	<file>	output file
	-pid	<str>	wild parent id may not
	-bid	<str>	mut bulk id
	-pdep	<str>	parent depth default 10
	-bdep	<str>	mut build depth default 10
	-popt	<str>	population type default F2
USAGE
	print $usage;
	exit;
}

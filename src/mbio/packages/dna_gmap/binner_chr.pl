#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($bin,$vcf,$mark,$fOut);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"vcf:s"=>\$vcf,
	"bin:s"=>\$bin,
	"mark:s"=>\$mark,
	"out:s"=>\$fOut
			) or &USAGE;
&USAGE unless ($vcf and $bin and $mark );
my %mark;
open In,$mark;
while (<In>){
	next if ($_=~/^#/|| /^$/);
	chomp;
	my ($id,@undi)=split(/\t/,$_);
	my ($chr,$pos)=split(/\_/,$id);
	$chr=(split(/\D+/,$chr))[-1];
	$mark{$chr}{$pos}=$pos;
}
close In;

my %stat;
open In,$vcf;
if($vcf=~/gz$/){
	close In;
	open In,"gunzip -c $vcf|";
}
while (<In>){
	chomp;
	next if ($_ eq ""|| /^$/|| /^#/|| /^CHROM/);
	my ($chrid,$posid,$id,$ref,$alt,@undi)=split(/\t/,$_,7);#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
	$chrid=(split(/\D+/,$chrid))[-1];
	$alt=(split(/\,/,$alt))[0];
	foreach my$chr(keys %mark){
		if($chrid eq $chr){
			foreach my $pos(keys %{$mark{$chr}}){
				if($posid eq $pos){
					if(length($ref)=="1" && length($alt)=="1"){
						$stat{$chr}{snp}++ ;
					}else{
						$stat{$chr}{indel}++ ;
					}
				}
			}
		}
	}
	#$stat{$chr}{$pos}{ref}=length($ref);
	#$stat{$chr}{$pos}{alt}=length($alt);
}
close In;
open In,$bin;
while(<In>){
	chomp;
	next if($_ eq ""|| /^$/|| /^#/);
	my($binnr,@undi)=split(/\t/,$_);
	my ($chr,$wind)=split(/\_/,$binnr);
	$chr=(split(/\D+/,$chr))[-1];
	$stat{$chr}{winnr}++ ;
}
close In;
open Out,">$fOut";
print Out "#Chromosome ID\tBin Number\tSNP Number Per Bin\tIndel Number Per Bin\n";
foreach my $chr (sort{$a<=>$b}keys %stat){
	my $persnp=sprintf("%.0f",$stat{$chr}{snp}/$stat{$chr}{winnr});
	my $perindel=sprintf("%.0f",$stat{$chr}{indel}/$stat{$chr}{winnr});
	print Out "chr$chr\t$stat{$chr}{winnr}\t$persnp\t$perindel\n";
	#	if($mark{$lg}{$id}{pos} eq $mark{$lg}{$id}{vcpos}){
	#		if ($mark{$lg}{$id}{ref} eq $mark{$lg}{$id}{alt}){
	#			$mark{$lg}{snp}++ ;
	#		}else{
	#			$mark{$lg}{indel}++ ;
	#		}
	#	}
	#}
	#print Vcf "$lg\t$mark{$lg}{snp}\t$mark{$lg}{indel}\n";
}
close Out;
#######################################################################################
#print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub USAGE {#
        my $usage=<<"USAGE";
Contact:        chongqing.shi\@majorbio.com;
Script:			$Script
Description:
	fq thanslate to fa format
	eg:
	perl $Script -i -o -k -c

Usage:
  Options:
  -vcf	<file>	input pop.vcf file
  -bin	<stri>	gmap's total.bin.marker	
  -mark	<file>  input filtered.marker
  -out	<file>	output result file
  -h         Help

USAGE
        print $usage;
        exit;
}

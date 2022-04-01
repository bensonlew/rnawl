#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($vcf,$mark,$fOut,$type);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"vcf:s"=>\$vcf,
	"mark:s"=>\$mark,
	"out:s"=>\$fOut,
	"type:s"=>\$type,
			) or &USAGE;
&USAGE unless ($vcf and $mark );
$type||="ALL";
my %info;
my %stat;
my %reg;
open In,$mark;
while (<In>){
	chomp;
	next if ($_ =~ /^#/|| /^$/);
	my ($ids,$genotype,@undi)=split(/\t/,$_,3);
	$reg{$ids}{genotype}=$genotype;
	if($genotype eq "abxcd"){
		$stat{abxcd}++;
	}elsif($genotype eq "aaxbb"){
		$stat{aaxbb}++;
	}elsif($genotype eq "ccxab"){
		$stat{ccxab}++;
	}elsif($genotype eq "abxcc"){
		$stat{abxcc}++;
	}elsif($genotype eq "efxeg"){
		$stat{efxeg}++;
	}elsif($genotype eq "hkxhk"){
		$stat{hkxhk}++;
	}elsif($genotype eq "nnxnp"){
		$stat{nnxnp}++;
	}elsif($genotype eq "lmxll"){
		$stat{lmxll}++;
	}
}
close In;
open Out,">$fOut";
if($type eq "ALL"){
	my %for;
	print Out "Type\tSNP Number\tInDel Number\n";
	open In,$vcf;
	if ($vcf=~/gz$/){
		close In;
		open In,"gunzip -c $vcf|";
	}
	while (<In>){
		chomp;
		next if ($_ eq ""|| /^$/|| /^#/);
		next if ($_ =~/CHROM/);
		my ($chr,undef,$id,$ref,$alt,@undi)=split(/\t/,$_,7);#CHROM     POS     ID      REF     ALT     QUAL    FILTER  INFO
		$for{$id}{format}=join("\,",$ref,$alt);
	}
	close In;
	foreach my$id(keys %for){
		foreach my$ids(keys %reg){
			if($ids eq $id){
				my($ref,$alt)=split(/\,/,$for{$id}{format},2);
				my $marktype;
				if(length$ref eq length$alt){
					$marktype="SNP";
				}else{
					$marktype="InDel";
				}
				$info{aaxbb}{$marktype}++ if($reg{$ids}{genotype} eq "aaxbb");
				$info{abxcd}{$marktype}++ if($reg{$ids}{genotype} eq "abxcd");
				$info{abxcc}{$marktype}++ if($reg{$ids}{genotype} eq "abxcc");
				$info{ccxab}{$marktype}++ if($reg{$ids}{genotype} eq "ccxab");
				$info{efxeg}{$marktype}++ if($reg{$ids}{genotype} eq "efxeg");
				$info{hkxhk}{$marktype}++ if($reg{$ids}{genotype} eq "hkxhk");
				$info{nnxnp}{$marktype}++ if($reg{$ids}{genotype} eq "nnxnp");
				$info{lmxll}{$marktype}++ if($reg{$ids}{genotype} eq "lmxll");
			}else{
				next;
			}
		}
	}
	foreach my $info (sort keys %info) {
		print Out $info,"\t";
		my @out=();
		foreach my $marktype(sort keys %{$info{$info}}){
			#print Out $marktype,"\t",$info{$info}{$marktype},"\t";
			my $num=join(",",$marktype,$info{$info}{$marktype});
			push @out,$num;
		}
		if(scalar@out eq "2"){
			my($snp,$indel);
			for(@out){
				my($markty,$number)=split(/\,/,$_);
				if($markty eq "SNP"){
					$snp=$number;
				}else{
					$indel=$number;
				}
			}
			print Out "$snp\t$indel\n";
		}else{
			my $out=shift@out;
			my($markty,$number)=split(/\,/,$out);
			if($markty eq "SNP"){
				print Out $number,"\t0\n";
			}else{
				print Out "0\t$number\n";
			}
		}
	}
}else{
	print Out "TYPE\t$type Number\n";
	foreach my $stat (sort keys %stat) {
		print Out "$stat\t$stat{$stat}\n";
	}
}
close Out;

#######################################################################################
#print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub USAGE {#
        my $usage=<<"USAGE";
Contact:        caixia.tian\@majorbio.com;
Script:			$Script
Description:
	fq thanslate to fa format
	eg:
	perl $Script -i -o -k -c

Usage:
  Options:
  -vcf	<file>	input pop.vcf file
  -mark	<file>  input marker file 
  -out	<file>	output result file
  -type	<str>	ALL or SNP or InDel
  -h         Help

USAGE
        print $usage;
        exit;
}

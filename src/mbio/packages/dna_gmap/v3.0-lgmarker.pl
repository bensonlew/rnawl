#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($lg,$vcf,$mark,$out);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"vcf:s"=>\$vcf,
	"lg:s"=>\$lg,
	"out:s"=>\$out,
			) or &USAGE;
&USAGE unless ($vcf and $lg and $out);

my %stat;
open VCF,$vcf;
if($vcf=~/gz$/){
	close VCF;
	open VCF,"gunzip -c $vcf|";
}
while (<VCF>){
	chomp;
	next if ($_ eq ""|| /^$/|| /^#/|| /^CHROM/);
	my (undef,undef,$id,$ref,$alt,@undi)=split(/\t/,$_,7);#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
	$stat{$id}{info}=join("\,",$ref,$alt);
}
close VCF;			

my %for;

open OUT,">$out";
print OUT "LG ID\tMarker Number\tSNP Number\tInDel Number\n";
open In,$lg;
$/=">";
while(<In>){
	chomp;
	next if ($_ eq ""|| /^#/ || /^$/);
	my($lgid,$lgnum,@marker)=split(/\s+/,$_);
	my $formation=join(",",$lgnum,@marker);
	$for{$lgid}=$formation;
}
close In;

foreach my$lgid(sort{$a<=>$b} keys%for){
	my($num,@marker)=split(/\,/,$for{$lgid});
	#print $num,"\n",join(",",@marker),"\n";
	my $snp=0;
	my $indel=0;
	print OUT "$lgid\t$num\t";
	foreach my$id(keys %stat){
			foreach(@marker){
				next if($_ ne $id);
				my($ref,$alt)=split(/\,/,$stat{$id}{info},2);
				if(length$ref eq length$alt){
					$snp++ ;
				}else{
					$indel++ ;
				}
			}
		}
		print OUT "$snp\t$indel\n";
	}
close OUT;
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
  -lg	<file>	Total.lg file
  -mark	<file>  if bin,input filtered.marker file
  -out	<file>	output result file
  -h         Help

USAGE
        print $usage;
        exit;
}

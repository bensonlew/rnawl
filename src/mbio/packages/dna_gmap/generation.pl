#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fIn,$fOut,$mark);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"vcf:s"=>\$fIn,
	"mark:s"=>\$mark,
	"out:s"=>\$fOut,
			) or &USAGE;
&USAGE unless ($fIn and $fOut and $mark);
my %stat;
my (@sample,@marker);
open In,$mark;
while(<In>){
	chomp;
	next if($_ eq ""|| /^$/);
	if(/^#Mark/){
		(undef,undef,@sample)=split(/\t/,$_);
	}else{
		my ($marker,$ptype,@samtype)=split(/\t/,$_);
		for(my$n=0;$n<scalar@sample;$n++){
			my $sample=$sample[$n];
			$stat{$marker}{$sample}=$samtype[$n];
		}
		push @marker,$marker;
	}
}
close In;
my %region;
my @Indi;
open In,$fIn;
if($fIn=~/gz$/){
        close In;
        open In,"gunzip -c $fIn|";
}
while (<In>){
	chomp;
	next if($_ eq ""|| /^$/|| /^##/);
	my ($chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$format,@indi)=split(/\t/,$_);
	if(/^#CHROM/){
		@Indi=@indi;
	}
	foreach my $marker(sort keys%stat){
		if($marker eq $id){
			for(my$i=0;$i<scalar@Indi;$i++){
				my $sampleid=$Indi[$i];
				foreach my $sample(sort keys%{$stat{$marker}}){
					if($sample eq $sampleid){
						if($stat{$marker}{$sample} eq "--"){
							$region{$sample}{miss}++ ;
						}
						my @format=split(/\:/,$format);
						for(my$j=0;$j<scalar@format;$j++){
							my @info=split(/\:/,$indi[$i]);
							if($format[$j] eq "DP"){
								$region{$sample}{depth} += $info[$j];
								#print $sample,"\t",$region{$sample}{miss},$info[$j],"\n";
							}
						}
					}
				}
			}
		}
	}
}
close In;
open Out,">$fOut";
print Out scalar@marker,"\t",scalar@sample,"\n";
foreach my$sample(sort keys%region){
	$region{$sample}{depth}||=0;
	$region{$sample}{miss}||=0;
	print Out "$sample\t",sprintf("%.2f",$region{$sample}{depth}/(scalar@marker)),"\t",sprintf("%.2f",100*$region{$sample}{miss}/(scalar@marker)),"\n";
	#print Out "$sample\t$depth\t$miss\n";
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
  -vcf	<file>	input pop.final.vcf
  -mark	<file>	pop.filtered.markr	
  -out	<file>	output result file
  -h         Help

USAGE
        print $usage;
        exit;
}

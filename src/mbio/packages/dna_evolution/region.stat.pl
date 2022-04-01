#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($vcf,$out,$anno,$region,$pvalue);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"vcf:s"=>\$vcf,
	"out:s"=>\$out,
	"ann:s"=>\$anno,
	"region:s"=>\$region,
	"p:s"=>\$pvalue,
			) or &USAGE;
&USAGE unless ($out);
open IN,$region;
my %region;
while (<IN>) {
	chomp;
	next if ($_ eq "" || /^$/ || /^#/);
	my ($chr,$start,$end)=split/\s+/;
	$region{$chr}{join("\t",sort{$a<=>$b}($start,$end))}=1;
}
close IN;
my %stat;
open IN,$pvalue;
my $PVAL =1;
while (<IN>) {
	chomp;
	s/\"//g;
	s/\s+//g;
	next if ($_ eq ""|| /^$/ || /^SNP/);
	my($chr,undef,$pos,undef,$pval)=split/\,/;
	next if ($pval eq "NA");
	#print $chr,"\t",$pos,"\t",$pval,"\n";
	foreach my $region (sort keys %{$region{$chr}}) {
		my ($start,$end)=split/\s+/,$region;
		if ($pos >= $start && $pos <= $end) {
			#print $chr,"\t",$pos,"\t",$pval,"\n";
			$PVAL = $pval if ($pval < $PVAL);
		}
		$stat{$chr}{$region}{pval}=$PVAL;
	}
}
close IN;
#print Dumper \%stat;die;
open IN,$vcf;
while (<IN>) {
	chomp;
	next if ($_ eq "" || /^$/ || /^#/);
	my ($chr,$pos,$id,$ref,$alt,$info)=split/\s+/,$_,6;
	my @ale = split(/\,/,join(",",$ref,$alt));
	my %len;
	foreach my $ale (@ale) {
		$len{length($ale)}=1;
	}
	my $type = "SNP";
	$type = "INDEL" if (scalar keys %len > 1);
	foreach my $region (sort keys %{$region{$chr}}) {
		my ($start,$end)=split/\s+/,$region;
		if ($pos >= $start && $pos <= $end) {
			$stat{$chr}{$region}{SNP}++ if ($type eq "SNP");
			$stat{$chr}{$region}{INDEL}++ if ($type eq "INDEL");
		}
	}
}
close IN;
if ($anno=~/pop.summary/) {
	open IN,$anno;
	while (<IN>) {
		chomp;
		next if ($_ eq "" || /^$/ || /^#/);
		my ($Gene_name,$Gene_id,$Transcript_id,$Bio_Type,$chr,$Pos1,$Pos2,$High,$Moderate,$Low,$Modifier,$nrid,$nranno,$uniid,$unianno,$koid,$koanno,$goid,$goanno,$egid,$express)=split(/\t/,$_);
		foreach my $region (sort keys %{$region{$chr}}) {
			my ($pos3,$pos4)=split(/\s+/,$region);
			if (($Pos1 > $pos3 && $Pos1 <$pos4) ||($Pos2 > $pos3 && $Pos2 < $pos4) || ($pos3 > $Pos1 && $pos3 < $Pos2) || ($pos4 > $Pos1 && $pos4 < $Pos2)) {
				$stat{$chr}{$region}{gene}++;
			}
		}
	}
}elsif ($anno=~/anno.summary/) {
	open IN,$anno;
	while (<IN>) {
		chomp;
		next if ($_ eq "" || /^$/ || /^#/);
		my ($geneid,$nrid,$nranno,$uniid,$unianno,$koid,$koanno,$goid,$goanno,$egid,$express)=split(/\t/,$_);
		my $REGION=(split(/\|/,$geneid))[-1];
		$REGION=~s/^://g;
		my ($chr,$Pos1,$Pos2)=split/\:/,$REGION;
		foreach my $region (sort keys %{$region{$chr}}) {
			my ($pos3,$pos4)=split(/\s+/,$region);
			if (($Pos1 > $pos3 && $Pos1 <$pos4) ||($Pos2 > $pos3 && $Pos2 < $pos4) || ($pos3 > $Pos1 && $pos3 < $Pos2) || ($pos4 > $Pos1 && $pos4 < $Pos2)) {
				$stat{$chr}{$region}{gene}++;
			}
		}
	}
	close IN;
}
open OUT,">$out";
foreach my $chr (sort keys %stat) {
	foreach my $region (sort keys %{$stat{$chr}}) {
		$stat{$chr}{$region}{gene}||=0;
		$stat{$chr}{$region}{INDEL}||=0;
		$stat{$chr}{$region}{SNP}||=0;
		print OUT join("\t",$chr,$region,$stat{$chr}{$region}{pval},$stat{$chr}{$region}{gene},$stat{$chr}{$region}{INDEL},$stat{$chr}{$region}{SNP}),"\n";
	}
}
close OUT;
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
	"vcf:s"=>\$vcf,   vcf file
	"out:s"=>\$out,   output file 
	"ann:s"=>\$anno,	
	"region:s"=>\$region,
	"p:s"=>\$pvalue,
  -h         Help

USAGE
        print $usage;
        exit;
}

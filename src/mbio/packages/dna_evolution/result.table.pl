#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($vcf,$out,$ann,$region,$fst,$pi1,$pi2,$tajimad1,$tajimad2,$win);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"vcf:s"=>\$vcf,
	"out:s"=>\$out,
	"ann:s"=>\$ann,
	"region:s"=>\$region,
	"fst:s"=>\$fst,
	"pi1:s"=>\$pi1,
	"pi2:s"=>\$pi2,
	"tajimad1:s"=>\$tajimad1,
	"tajimad2:s"=>\$tajimad2,
	"win:s"=>\$win,
			) or &USAGE;
&USAGE unless ($out);
$win||=10000;
my $pop1=(split/\./,basename($pi1))[0];
my $pop2=(split/\./,basename($pi2))[0];
$pop1 = "pop".$pop1;
$pop2 = "pop".$pop2;

open IN,$region;
my %region;
while (<IN>) {
	chomp;
	next if ($_ eq "" || /^$/ ||/^#/);
	my ($chr,$pos1,$pos2)=split/\s+/;
	$region{$chr}{join("\t",$pos1,$pos2)}=1;
}
close IN;
open IN,$vcf;
my %stat;
while (<IN>) {
	chomp;
	next if ($_ eq "" || /^$/ || /^#/);
	my ($chr,$pos,$id,$ref,$alt,undef)=split/\s+/;
	my @ale = split(/\,/,join(",",$ref,$alt));
	my %len;
	foreach my $ale (@ale) {
		$len{length($ale)}=1;
	}
	foreach my $region (sort keys %{$region{$chr}}) {
		my ($pos1,$pos2)=split/\s+/,$region;
		if ($pos >= $pos1 && $pos <= $pos2) {
			if (scalar keys %len eq 1) {
				$stat{$chr}{$region}{SNP}++;
			}else{
				$stat{$chr}{$region}{INDEL}++;
			}
		}
	}
}
close IN;
#print Dumper \%stat;die;
if (defined $ann) {
	open IN,$ann;
	while (<IN>) {
		chomp;
		next if ($_ eq ""|| /^$/ || /^#/);
		my ($gene,$id,$tid,$type,$chr,$pos1,$pos2,undef)=split/\s+/;
		foreach my $region (sort keys %{$region{$chr}}) {
			my ($pos3,$pos4)=split/\s+/,$region;
			if (($pos1 > $pos3 && $pos1 <$pos4)||($pos2 > $pos3 && $pos2 < $pos4)||($pos3 > $pos1 && $pos3 < $pos2)||($pos4 > $pos1 && $pos4 < $pos2)) {
				$stat{$chr}{$region}{gene}++;
			}
		}
	}
	close IN;
}

#print Dumper \%stat;die;
open IN,$fst;
while (<IN>) {
	chomp;
	next if ($_ eq ""|| /^$/ || /^CHROM/);
	my ($chr,$pos1,$pos2,$vnum,$Fst,$mean_fst)=split/\s+/;
	foreach my $region (sort keys %{$region{$chr}}) {
		my ($pos3,$pos4)=split/\s+/,$region;
		if (($pos1 > $pos3 && $pos1 <$pos4)||($pos2 > $pos3 && $pos2 < $pos4)||($pos3 > $pos1 && $pos3 < $pos2)||($pos4 > $pos1 && $pos4 < $pos2)) {
			$stat{$chr}{$region}{fstnum}++;
			$stat{$chr}{$region}{totalfst}+=$mean_fst;
		}
	}
}
close IN;
#print Dumper \%stat;die;
open IN,$pi1;
while (<IN>) {
	chomp;
	next if ($_ eq ""|| /^$/ || /^CHROM/);
	my ($chr,$pos1,$pos2,$vnum,$pi)=split/\s+/;
	foreach my $region (sort keys %{$region{$chr}}) {
		my ($pos3,$pos4)=split/\s+/,$region;
		if (($pos1 > $pos3 && $pos1 <$pos4)||($pos2 > $pos3 && $pos2 < $pos4)||($pos3 > $pos1 && $pos3 < $pos2)||($pos4 > $pos1 && $pos4 < $pos2)) {
			$stat{$chr}{$region}{pi1num}++;
			$stat{$chr}{$region}{totalpi1}+=$pi;
		}
	}
}
close IN;
#print Dumper \%stat;die;
open IN,$pi2;
while (<IN>) {
	chomp;
	next if ($_ eq ""|| /^$/ || /^CHROM/);
	my ($chr,$pos1,$pos2,$vnum,$pi)=split/\s+/;
	foreach my $region (sort keys %{$region{$chr}}) {
		my ($pos3,$pos4)=split/\s+/,$region;
		if (($pos1 > $pos3 && $pos1 <$pos4)||($pos2 > $pos3 && $pos2 < $pos4)||($pos3 > $pos1 && $pos3 < $pos2)||($pos4 > $pos1 && $pos4 < $pos2)) {
			$stat{$chr}{$region}{pi2num}++;
			$stat{$chr}{$region}{totalpi2}+=$pi;
		}
	}
}
close IN;
#print Dumper \%stat;die;
open IN,$tajimad1;
while (<IN>) {
	chomp;
	next if ($_ eq ""|| /^$/ || /^CHROM/);
	my ($chr,$pos1,$vnum,$tajimad)=split/\s+/;
	my $pos2 = $pos1 + 10000;
	$pos1 = $pos1 -1;
	$tajimad = 0 if ($tajimad =~/nan/);
	foreach my $region (sort keys %{$region{$chr}}) {
		my ($pos3,$pos4)=split/\s+/,$region;
		if ( ($pos1 > $pos3 && $pos1 <$pos4)||($pos2 > $pos3 && $pos2 < $pos4)||($pos3 > $pos1 && $pos3 < $pos2)||($pos4 > $pos1 && $pos4 < $pos2)) {
				$stat{$chr}{$region}{tajimad1num}++;
				$stat{$chr}{$region}{totaltajimad1}+=$tajimad;
		}
	}
}
close IN;
#print Dumper \%stat;die;
open IN,$tajimad2;
while (<IN>) {
	chomp;
	next if ($_ eq ""|| /^$/ || /^CHROM/);
	my ($chr,$pos1,$vnum,$tajimad)=split/\s+/;
	my $pos2 = $pos1 + 10000;
	$pos1 = $pos1 -1;
	$tajimad = 0 if ($tajimad =~/nan/);
	foreach my $region (sort keys %{$region{$chr}}) {
		my ($pos3,$pos4)=split/\s+/,$region;
		if (($pos1 > $pos3 && $pos1 <$pos4)||($pos2 > $pos3 && $pos2 < $pos4)||($pos3 > $pos1 && $pos3 < $pos2)||($pos4 > $pos1 && $pos4 < $pos2)) {
			$stat{$chr}{$region}{tajimad2num}++;
			$stat{$chr}{$region}{totaltajimad2}+=$tajimad;
		}
	}
}
close IN;
#print Dumper \%stat;die;
open OUT,">$out";
print OUT "Chr ID\tStart\tEnd\tGene Number\tSNP Number\tIndel Number\tAverage FST\t$pop1 PI\t$pop2 PI\t$pop1 Tajima\'D\t$pop2 Tajima\'D\n";
foreach my $chr (sort keys %stat) {
	foreach my $region (sort keys %{$stat{$chr}}) {
		$stat{$chr}{$region}{gene}||=0;
		$stat{$chr}{$region}{SNP}||=0;
		$stat{$chr}{$region}{INDEL}||=0;
		if (exists $stat{$chr}{$region}{fstnum}){
			my $mean_fst =$stat{$chr}{$region}{totalfst}/$stat{$chr}{$region}{fstnum};
			my $mean_pi1 =$stat{$chr}{$region}{totalpi1}/$stat{$chr}{$region}{pi1num};
			my $mean_pi2 =$stat{$chr}{$region}{totalpi2}/$stat{$chr}{$region}{pi2num};
			my $mean_tajimad1 =$stat{$chr}{$region}{totaltajimad1}/$stat{$chr}{$region}{tajimad1num};
			my $mean_tajimad2 =$stat{$chr}{$region}{totaltajimad2}/$stat{$chr}{$region}{tajimad2num};
			my ($pos1,$pos2)=split/\s+/,$region;
			print  OUT join("\t",$chr,$pos1,$pos2,$stat{$chr}{$region}{gene},$stat{$chr}{$region}{SNP},$stat{$chr}{$region}{INDEL},$mean_fst,$mean_pi1,$mean_pi2,$mean_tajimad1,$mean_tajimad2),"\n";
		}
		# my $mean_fst =$stat{$chr}{$region}{totalfst}/$stat{$chr}{$region}{fstnum};
		# my $mean_pi1 =$stat{$chr}{$region}{totalpi1}/$stat{$chr}{$region}{pi1num};
		# my $mean_pi2 =$stat{$chr}{$region}{totalpi2}/$stat{$chr}{$region}{pi2num};
		# my $mean_tajimad1 =$stat{$chr}{$region}{totaltajimad1}/$stat{$chr}{$region}{tajimad1num};
		# my $mean_tajimad2 =$stat{$chr}{$region}{totaltajimad2}/$stat{$chr}{$region}{tajimad2num};
		# my ($pos1,$pos2)=split/\s+/,$region;
		# print  OUT join("\t",$chr,$pos1,$pos2,$stat{$chr}{$region}{gene},$stat{$chr}{$region}{SNP},$stat{$chr}{$region}{INDEL},$mean_fst,$mean_pi1,$mean_pi2,$mean_tajimad1,$mean_tajimad2),"\n";
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
	"vcf:s"=>\$vcf,
	"out:s"=>\$out,
	"ann:s"=>\$ann,
	"region:s"=>\$region,
	"fst:s"=>\$fst,
	"pi1:s"=>\$pi1,
	"pi2:s"=>\$pi2,
	"tajimad1:s"=>\$tajimad1,
	"tajimad2:s"=>\$tajimad2,
  -h         Help

USAGE
        print $usage;
        exit;
}

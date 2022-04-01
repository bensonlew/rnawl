#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($vcf,$out,$tag,$key,$config);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"vcf:s"=>\$vcf,
	"out:s"=>\$out,
	"tag:s"=>\$tag,
	"key:s"=>\$key,
	"config:s"=>\$config,
			) or &USAGE;
&USAGE unless ($vcf and $out and $key);
mkdir $out if(!-d $out);

my(@Alle_Number,$mode,$gtp,$mtp,$maxdepth,$mindepth);
my %stat;
my $round;
my $test;
open IN,$config;
while(<IN>){
	chomp;
	next if($_ eq ""|| /^$/|| /^#/|| /^\{/|| /^\}/);
	if($_=~/Alle_Number/){
		@Alle_Number = split(/,/,(split(/\:/,$_))[-1]);
	}
	next if($_=~/Model/);
	next if($_=~/Selecttion/);
	if($_=~/mode:(\w+)/){
		#print $1,"\n";die;
		$mode=$1;next;
	}
	if($_ =~ /marktype:(\w+)/){
		#print $1,"\n";die;
		$mtp=$1;next;
	}
	if($_ =~ /genetype:(\w+)/){
		#print $1,"\n";die;
		$gtp=$1;next;
	}
	if($_ =~ /depth:(.*),(.*)/){
		$mindepth=$1;
		$maxdepth=$2;
		#print "$1\t$2\n";die;
		next;
	}
	if($_ =~ /Sample:(\S+),(\S+)/){
		$round++;
		$stat{$round}{name1}=$1;
		$stat{$round}{name2}=$2;
		#print "$1\t$2\n";die;
		next;
	}
}
close IN;

my %seq;
$/=">";
open SEQ,$tag;
while(<SEQ>){
	chomp;
	next if($_ eq ""|| /^$/|| /^#/);
	my($seqid,$squence,$undi)=split(/\n/,$_);
	$seq{$seqid}=$squence;
}
close SEQ;

$/="\n";
my @indi;
my %snp;
open In,$vcf;
open OUT,">$out/$key.table.xls";
print OUT "Tag ID\tPos\tSNP ID\tREF\tALT\t";
open OUT2,">$out/$key.vcf";
foreach my$round(sort{$a<=>$b}keys%stat){
	print OUT "$stat{$round}{name1} Genotype\t$stat{$round}{name1} Allele Depth\t$stat{$round}{name2} Genotype\t$stat{$round}{name2} Allele Depth";
}
print OUT "\tTag seq\n";
#open STAT,">$out/$key.stat";
#print STAT "Analysis ID\tSNP Number\tAverage Depth\tMiss Ratio\n";
my($outsnp,$outnum,$outad);
while(<In>){
	chomp;
	next if($_ eq ""|| /^$/|| /^##/);
	my ($chr,$pos,$id,$ref,$alt,$qual,$Filter,$indo,$format,@info)=split(/\t/,$_);#CHROM  POS     ID      REF     ALT     QUAL	FILTER  INFO    FORMAT  A01     A02     A03
	if (scalar @indi ==0) {
		push @indi,@info;
		for(my$i=0;$i<scalar@indi;$i++){
			foreach my$round(keys%stat){
				$stat{$round}{name1}=join("\t",$stat{$round}{name1},$i) if($indi[$i] eq $stat{$round}{name1});
				$stat{$round}{name2}=join("\t",$stat{$round}{name2},$i) if($indi[$i] eq $stat{$round}{name2});
			}
		}
		print OUT2 $_,"\n";
		#print join("\n",scalar@indi,@info),"\n";die;
		next;
	}else{
		my @alle=split(/\,/,join(",",$ref,$alt));
		my $allenumber=scalar@alle;
		my $d = grep {$_ eq $allenumber }@Alle_Number;#检测是否符合需要的等位基因数
		#print "alle number:$allenumber\n",join("\,",@Alle_Number),"\n",$d,"\n";die;
		next if($d eq "0");
		my($gt,$genotype);
		my @out;
		foreach my$round(sort{$a<=>$b}keys%stat){
			my $sroundmiss=0;
			my($s1gt,$s2gt,$s1ad,$s2ad);
			my $s1nr=(split(/\t/,$stat{$round}{name1}))[1];
			my $s2nr=(split(/\t/,$stat{$round}{name2}))[1];
			#print "$stat{$round}{name1}\t$s1nr\n$stat{$round}{name2}\t$s2nr\n";die;
			my @infos1=split(/\:/,$info[$s1nr]);
			my @infos2=split(/\:/,$info[$s2nr]);
			next if($infos1[0] eq "./.");#样本间计较，缺失直接过滤掉
			next if($infos2[0] eq "./.");
			next if(($mtp eq "diff") && ($infos1[0] eq $infos2[0]));
			next if(($mtp eq "same") && ($infos1[0] ne $infos2[0]));
			#next if(($mtp eq "same") && ($infos1[0] eq "0/0"));
			
			next if($infos1[1] < $mindepth || $infos1[1] > $maxdepth);
			next if($infos2[1] < $mindepth || $infos1[1] > $maxdepth);

			my @geno1=split(/\//,$infos1[0]);
			my @geno2=split(/\//,$infos2[0]);
			#push @out,join("\t",join("",$alle[$geno1[0]],$alle[$geno1[1]]),$infos1[2],join("",$alle[$geno2[0]],$alle[$geno2[1]]),$infos2[2]);

			#print $gtp,"\n";
			#die;
			next if(($gtp eq "hete") && ($geno1[0] eq $geno1[1]));
			next if(($gtp eq "hete") && ($geno2[0] eq $geno2[1]));
			next if(($gtp eq "homo") && ($geno1[0] ne $geno1[1]));
			next if(($gtp eq "homo") && ($geno2[0] ne $geno2[1]));
			push @out,join("\t",join("",$alle[$geno1[0]],$alle[$geno1[1]]),$infos1[2],join("",$alle[$geno2[0]],$alle[$geno2[1]]),$infos2[2]);
		}
		next if(scalar@out ne $round);
		my ($tagid,$snppos)=split(/\_/,$id);
		print OUT join("\t",$tagid,$snppos,$id,$ref,$alt,@out,$seq{$tagid}),"\n";
		print OUT2 $_,"\n";
	}
}
close In;
close OUT;
close OUT2;
#print STAT join("\t",$key,$outsnp,sprintf("%.2f",$outad/($outsnp*$sround*2)),sprintf("%.2f",$outnum/($outsnp*$gnum))),"\n";
#close STAT;
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub uniq {
  my %seen;
  return grep { !$seen{$_}++ } @_;
  #grep { ++$count{ $_ } < 2; };
}
#my @uniq_times = grep { ++$count{ $_ } < 2; } @array;#数组元素去重
#my %hash_a = map{$_=>1} @a;
#my %hash_b = map{$_=>1} @b;
#my %merge_all = map {$_ => 1} @a,@b;
#my @common = grep {$hash_a{$_}} @a;#交集
#my @a_only = grep {!$hash_b{$_}} @a;#a独有的
#my @b_only = grep {!$hash_a{$_}} @b;#b独有的
#my @merge = keys (%merge_all);     #并集
sub USAGE {#
        my $usage=<<"USAGE";
Contact:        long.huang\@majorbio.com;
Script:			$Script
Description:

Usage:
  Options:
  -vcf	<file>	input vcf list file
  -out	<dir>	output dir
  -tag	<file>	tag file
  -key	<str>	analysis id
  -config	<file>	config file
  -h         Help

USAGE
        print $usage;
        exit;
}

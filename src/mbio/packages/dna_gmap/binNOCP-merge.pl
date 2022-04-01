#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($female,$male,$dOut,$Key);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"f:s"=>\$female,
	"m:s"=>\$male,
	"o:s"=>\$dOut,
	"k:s"=>\$Key,
			) or &USAGE;
&USAGE unless ($female and $male and $dOut and $Key);
open In,$female;
my @Head;
my %female;
my %pos;#��¼���еĶϵ�
my %fpos;
while (<In>) {
	chomp;
	next if ($_ eq "" || /^$/);
	if (/^#/) {
		(undef,undef,undef,undef,@Head)=split;
	}else {
		my ($chr,$binID,$bin,$phase,@info)=split;
		my ($p1,$p2)=split(/\-/,$bin);
		$p2++;
		$bin=$p1."-".$p2;
		$female{$chr}{$bin}=join("\t",$phase,@info);
		$pos{$chr}{$p1}=$bin;
		$pos{$chr}{$p2}=$bin;
		$fpos{$chr}{$p1}=$bin;
		$fpos{$chr}{$p2}=$bin;
	}
}
close In;
open In,$male;
my %male;
my %mpos;
while (<In>) {
	next if ($_ eq "" || /^$/ || /^#/);
	my ($chr,$binID,$bin,$phase,@info)=split;
	my ($p1,$p2)=split(/\-/,$bin);
	$p2++;
	$bin=$p1."-".$p2;
	$male{$chr}{$bin}=join("\t",$phase,@info);
	$pos{$chr}{$p1}=$bin;
	$pos{$chr}{$p2}=$bin;
	$mpos{$chr}{$p1}=$bin;
	$mpos{$chr}{$p2}=$bin;
}
close In;
my %merge;
open Bin,">$dOut/$Key.merge.detail";
print Bin "#chr\tpos\ttype\t",join("\t",@Head),"\n";
open Geno,">$dOut/$Key.bin.marker";
print Geno "#MarkerID\ttype\t",join("\t",@Head),"\n";
foreach my $chr (sort keys %pos) {
	my @pos=sort {$a<=>$b} keys %{$pos{$chr}};
	#shift @pos;pop @pos;
	my @mpos=sort {$a<=>$b} keys %{$mpos{$chr}};
	my @fpos=sort{$a<=>$b}keys %{$fpos{$chr}};
	my @range;
	my $binID=0;
		my %mtest;
		for (my $i=0;$i<@mpos;$i++) {#��male��bin�������նϵ�洢
			for (my $j=0;$j<@pos-1;$j++) {
				if ($pos[$j] <= $mpos[$i] && !exists $mtest{$pos[$j]}) {
					$mtest{$pos[$j]}=1;
					my $range="$pos[$j]-$pos[$j+1]";
					if ($pos[$j+1]==$mpos[$i]) {#male��ÿ���ϵ㱣������Զϵ�Ϊ��ʼ��bin�������Ҫ��1
						$merge{$chr}{$range}{male}=$male{$chr}{$mpos{$chr}{$mpos[$i-1]}};
					}else{
						$merge{$chr}{$range}{male}=$male{$chr}{$mpos{$chr}{$mpos[$i]}};
					}
					push @range,$range;
				}else{
					next;
				}
			}
		}
		my %ftest;
		my @range1;
		for (my $i=0;$i<@fpos;$i++) {#��female��bin�������նϵ�洢
			for (my $j=0;$j<@pos-1;$j++) {
				if ($pos[$j] <= $fpos[$i] && !exists $ftest{$pos[$j]}) {
					$ftest{$pos[$j]}=1;
					my $range="$pos[$j]-$pos[$j+1]";
					if ($pos[$j+1]==$fpos[$i]) {#female��ÿ���ϵ㱣������Զϵ�Ϊ��ʼ��bin�������Ҫ��1
						$merge{$chr}{$range}{female}=$female{$chr}{$fpos{$chr}{$fpos[$i-1]}};
					}else{
						$merge{$chr}{$range}{female}=$female{$chr}{$fpos{$chr}{$fpos[$i]}};
					}
					push @range1,$range;
				}else{
					next;
				}
			}
		}
		foreach my $bin (@range) {#���ս�nnxnp��lmxll����Ϊabxcd
			$binID++;
			my @out;
			my @outG;
			
			my ($phase1,@g1)=split(/\t/,$merge{$chr}{$bin}{male});
			my ($phase2,@g2)=split(/\t/,$merge{$chr}{$bin}{female});
			for (my $i=0;$i<@g1;$i++) {
				if ($g1[$i] == 1 && $g2[$i] == 1) {#np,lm
					push @out,"aa $g1[$i] $g2[$i]";
					push @outG,"aa";
				}elsif ($g1[$i] == 1 && $g2[$i] == -1) {#lm,nn
					push @out,"ab $g1[$i] $g2[$i]";
					push @outG,"ab";
				}elsif ($g1[$i] == 0 || $g2[$i] == 0) {#--
					push @out,"-- $g1[$i] $g2[$i]";
					push @outG,"--";
				}elsif ($g1[$i] == -1 && $g2[$i] == -1) {#ll,np
					push @out,"bb $g1[$i] $g2[$i]";
					push @outG,"bb";
				}elsif ($g1[$i] == -1 && $g2[$i] == 1) {#lm,nn
					push @out,"ab $g1[$i] $g2[$i]";
					push @outG,"ab";
				}
			}
			#print Bin scalar @out,"\t";
			print Bin $chr,"\t",$bin,"\t","aaxbb","\t",join("\t",@out),"\n";
			print Geno "$chr\_$bin\taaxbb\t",join("\t",@outG),"\n";
#			print Pos "Bin$binID\t$chr\t$bin\n";
		}
}
close Bin;
close Geno;
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub USAGE {#
        my $usage=<<"USAGE";
Contact:        long.huang\@majorbio.com;
Script:			$Script
Description:
	eg:
	perl $Script -i -o 

Usage:
  Options:
  -f	<file>	input female binphase
  -m	<file>	input male binphase
  -o	<dir>	output dir
  -k	<str>	output key of filename
  -h         Help

USAGE
        print $usage;
        exit;
}

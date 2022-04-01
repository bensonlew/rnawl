#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fIn,$dOut,$Key,$popt);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$fIn,
	"o:s"=>\$dOut,
	"k:s"=>\$Key,
	"p:s"=>\$popt,
			) or &USAGE;
&USAGE unless ($fIn and $dOut and $Key);
$popt||="F2";
open In,$fIn;
my %info;
my %pos;
my @Head;
while (<In>) {
	chomp;
	next if ($_ eq "" ||/^$/);
	if (/^#/) {
		(undef,undef,undef,@Head)=split(/\s+/,$_);
		next;
	}
	my ($chr,$win,undef,@info)=split;
	my $p1=(split(/\-/,$win))[0];
	$info{$chr}{$p1}=join("\t",@info);
	$pos{$chr}{$p1}=$win;
}
close In;
open Out,">$dOut/$Key.bin.phase";
open Info,">$dOut/$Key.bin.info";
print Out "#Chr\tBinID\tpos\tphase\t",join("\t",@Head),"\n";
print Info "#Chr\twinID\tbinID\tpos\tphase\t",join("\t",@Head),"\n";
foreach my $chr (sort keys %info) {
	if (scalar keys %{$info{$chr}}==1) {#一个scaffold只是一个windows，不作处理，直接为bin
		my $binID=0;
		my $p=(keys %{$info{$chr}})[0];
		print Out "$chr\tbin$binID\t$pos{$chr}{$p}\t1\t",$info{$chr}{$p},"\n";
	}else{#一个scaffold多个windows，首先判断两两间重组事件最少的连锁相（may be not correct）
		my @win=sort {$a<=>$b} keys %{$info{$chr}};
		my %phase;
		my %break;
		for (my $i=0;$i<@win;$i++) {##确定每个marker的连锁相，用当前染色体所有marker与其连锁相的加和，少数服从多数原则
			my @g1=split(/\s+/,$info{$chr}{$win[$i]});
			for (my $j=$i+1;$j<@win;$j++) {
				my @g2=split(/\s+/,$info{$chr}{$win[$j]});
				my $diff1=0;
				my $diff2=0;
				for (my $k=0;$k<@g1;$k++) {
					$diff1++ if ($g1[$k] != $g2[$k]  && $g1[$k]!=0 && $g2[$k]!=0);
					$diff2++ if ($g1[$k] != $g2[$k]*-1 && $g1[$k]!=0 && $g2[$k]!=0);
				}
				if ($diff1 < $diff2 || $popt ne "CP") {#不需要进行连锁相转换
					$phase{$win[$i]}++;
					$phase{$win[$j]}++;
					if ($j-$i==1 && $diff1 > 0) {
						$break{$i}=1;
					}
				}else {#需要进行连锁相转换
					$phase{$win[$i]}--;
					$phase{$win[$j]}--;
					if ($j-$i == 1 && $diff2 > 0) {
						$break{$i}=1;
					}
				}
			}
		}
		my @break=sort{$a<=>$b} keys %break;
		push @break,scalar @win-1;
		my $lastpos=0;
		my $start=0;
		my $binid=0;
		for (my $i=0;$i<@break;$i++) {
			if (!defined $break[$i]) {
				print scalar @break,"\t",$i;
				die;
			}
			my $end=(split("-",$pos{$chr}{$win[$break[$i]]}))[-1]-1;
			my $pos=$start."-".$end;
			$binid++;
			my %stat;
			my @bin;
			my $binphase=0;
			for (my $j=$lastpos;$j<=$break[$i];$j++) {
				$phase{$win[$j]}=-1 if ($phase{$win[$j]}<0);
				$phase{$win[$j]}=1 if ($phase{$win[$j]}>0);
				$binphase+=$phase{$win[$j]};
				print Info $chr,"\t",$win[$j],"\t",$binid,"\t",$pos{$chr}{$win[$j]},"\t",$phase{$win[$j]},"\t",$info{$chr}{$win[$j]},"\n";
				my @g1=split(/\s+/,$info{$chr}{$win[$j]});
				push @bin,$win[$j];
				for (my $k=0;$k<@g1;$k++) {
					if ($phase{$win[$j]} < 0) {
						$stat{$Head[$k]}-= $g1[$k]
					}else{
						$stat{$Head[$k]}+= $g1[$k]	
					}
				}
				$lastpos=$j+1;
			}
		#	print Info $chr,"\tBin$binid\t",join(":",@bin),"\t",join(":",%stat),"\n";
			#die Dumper @bin if ($chr eq "135");
			my @out;
			for (my $k=0;$k<@Head;$k++) {
				if (!exists $stat{$Head[$k]}) {
					die $chr;
				}
				if ($stat{$Head[$k]} > 0) {
					push @out,1;
				}elsif ($stat{$Head[$k]} < 0) {
					push @out,-1;
				}else{
					push @out,0;
				}
			}
			if ($binphase > 0) {
				$binphase = 1;
			}elsif ($binphase<0) {
				$binphase =-1;
			}else{
				$binphase = 0
			}
			if ($i == scalar @break -1) {
				print Out "$chr\tBin$binid\t",$pos,"\t",$binphase,"\t",join("\t",@out),"\n";
			}else{
				print Out "$chr\tBin$binid\t",$pos,"\t",$binphase,"\t",join("\t",@out),"\n";
			}
			$start=$end+1;
		}
	}
}
close Out;
close Info;

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
  -i	<file>	input file name
  -o	<dir>	output dir
  -k	<str>	output keys of filename
  -h         Help

USAGE
        print $usage;
        exit;
}

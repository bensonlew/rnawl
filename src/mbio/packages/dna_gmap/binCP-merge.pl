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
my %pos;#记录所有的断点
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
	if (scalar @mpos == 0 && scalar @fpos != 0) {#单个scaffold上面只有nnxnp型
		pop @mpos;
		foreach my $pos (@fpos) {
			my $bin=$fpos{$chr}{$pos};
			$binID++;
			my ($phase,@info)=split("\t",$female{$chr}{$bin});
			my @out;
			my @outG;
			foreach my $geno (@info) {
				if ($geno eq "1") {
					push @out,"np 1";
					push @outG,"np";
				}elsif ($geno eq "0") {
					push @out,"-- 0";
					push @outG,"--";
				}else{
					push @out,"nn -1";				
					push @outG,"nn";
				}
			}
			print Bin $chr,"\t",$bin,"\t","nnxnp","\t",join("\t",@out),"\n";
			print Geno "$chr\_$bin\tnnxnp\t",join("\t",@outG),"\n";
		}
	}elsif (scalar @fpos ==0 && scalar @mpos !=0) {#单个scaffold上面只有lmxll型
		pop @fpos;
		foreach my $pos (@fpos) {
			$binID++;
			my $bin=$mpos{$chr}{$pos};
			my ($phase,@info)=split("\t",$male{$chr}{$bin});
			my @out;
			my @outG;
			foreach my $geno (@info) {
				if ($geno eq "1") {
					push @out,"lm 1";
					push @outG,"lm";
				}elsif ($geno eq "0") {
					push @out,"-- 0";
					push @outG,"--";
				}else{
					push @out,"ll -1";				
					push @outG,"ll";
				}
			}
			print Bin $chr,"\t",$bin,"\t","lmxll","\t",join("\t",@out),"\n";
			print Geno "$chr\_$bin\tlmxll\t",join("\t",@outG),"\n";
		}
	}else{
		my %mtest;
		for (my $i=0;$i<@mpos;$i++) {#将male的bin按照最终断点存储
			for (my $j=0;$j<@pos-1;$j++) {
				if ($pos[$j] <= $mpos[$i] && !exists $mtest{$pos[$j]}) {
					$mtest{$pos[$j]}=1;
					my $range="$pos[$j]-$pos[$j+1]";
					if ($pos[$j+1]==$mpos[$i]) {#male中每个断点保存的是以断点为起始的bin，因此需要减1
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
		for (my $i=0;$i<@fpos;$i++) {#将female的bin按照最终断点存储
			for (my $j=0;$j<@pos-1;$j++) {
				if ($pos[$j] <= $fpos[$i] && !exists $ftest{$pos[$j]}) {
					$ftest{$pos[$j]}=1;
					my $range="$pos[$j]-$pos[$j+1]";
					if ($pos[$j+1]==$fpos[$i]) {#female中每个断点保存的是以断点为起始的bin，因此需要减1
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
		foreach my $bin (@range) {#最终将nnxnp和lmxll整合为abxcd
			$binID++;
			my @out;
			my @outG;
			if (!exists $merge{$chr}{$bin}{female} && exists $merge{$chr}{$bin}{male}){#某个window中只包含父本数据，未包含母本数据
				my ($phase1,@g1)=split(/\t/,$merge{$chr}{$bin}{male});
				for (my $i=0;$i<@g1;$i++) {
					if ($g1[$i] == 1) {#lm
						push @out,"lm $g1[$i]";
						push @outG,"lm";
					}elsif ($g1[$i] == -1 ) {#ll
						push @out,"ll $g1[$i]";
						push @outG,"ll";
					}else {#--
						push @out,"-- $g1[$i]";
						push @outG,"--";
					}
				}
				print Bin $chr,"\t",$bin,"\t","lmxll","\t",join("\t",@out),"\n";
				print Geno "$chr\_$bin\tlmxll\t",join("\t",@outG),"\n";
			}elsif (!exists $merge{$chr}{$bin}{male} && exists $merge{$chr}{$bin}{female}) {#某个window中只包含父本数据，未包含母本数据
				my ($phase1,@g1)=split(/\t/,$merge{$chr}{$bin}{female});
				for (my $i=0;$i<@g1;$i++) {
					if ($g1[$i] == 1) {#lm
						push @out,"np $g1[$i]";
						push @outG,"np";
					}elsif ($g1[$i] == -1 ) {#ll
						push @out,"nn $g1[$i]";
						push @outG,"nn";
					}else {#--
						push @out,"-- $g1[$i]";
						push @outG,"--";
					}
				}
				print Bin $chr,"\t",$bin,"\t","nnxnp","\t",join("\t",@out),"\n";
				print Geno "$chr\_$bin\tnnxnp\t",join("\t",@outG),"\n";
			}elsif (!exists $merge{$chr}{$bin}{male} && !exists $merge{$chr}{$bin}{female}) {#某个bin中父本母本均无数据
				next;
			}else{
				my ($phase1,@g1)=split(/\t/,$merge{$chr}{$bin}{male});
				my ($phase2,@g2)=split(/\t/,$merge{$chr}{$bin}{female});
				for (my $i=0;$i<@g1;$i++) {
					if ($g1[$i] == 1 && $g2[$i] == 1) {#np,lm
						push @out,"bd $g1[$i] $g2[$i]";
						push @outG,"bd";
					}elsif ($g1[$i] == 1 && $g2[$i] == -1) {#lm,nn
						push @out,"bc $g1[$i] $g2[$i]";
						push @outG,"bc";
					}elsif ($g1[$i] == 0 || $g2[$i] == 0) {#--
						push @out,"-- $g1[$i] $g2[$i]";
						push @outG,"--";
					}elsif ($g1[$i] == -1 && $g2[$i] == 1) {#ll,np
						push @out,"ad $g1[$i] $g2[$i]";
						push @outG,"ad";
					}else{
						push @out,"ac $g1[$i] $g2[$i]";
						push @outG,"ac";
					}
				}
			#	print Bin scalar @out,"\t";
				print Bin $chr,"\t",$bin,"\t","abxcd","\t",join("\t",@out),"\n";
				print Geno "$chr\_$bin\tabxcd\t",join("\t",@outG),"\n";
			}
		}
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

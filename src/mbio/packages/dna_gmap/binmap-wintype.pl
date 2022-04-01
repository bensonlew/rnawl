#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fgeno,$dir,$Key,$fpos,$winsize,$stepsize);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$fgeno,
	"o:s"=>\$dir,
	"k:s"=>\$Key,
	"win:s"=>\$winsize,
	"step:s"=>\$stepsize,
			) or &USAGE;
&USAGE unless ($fgeno and $dir and $Key);
mkdir $dir if (!-d $dir);
$winsize||=500;
$stepsize||=$winsize/5;
#print $winsize,"\t",$stepsize,"\n";
$winsize = $winsize *1000;
$stepsize = $stepsize * 1000;
open In,$fgeno;
my @Head;
my %Phase;
my %pos;
while (<In>) {#��ȡmatrix�����ÿ�������Դ
	chomp;
	next if ($_ eq ""||/^$/);
	if (/^#/) {
		(undef,undef,@Head)=split(/\t/,$_);
	}else{
		my ($id,$type,@info)=split(/\t/,$_);
		my @posinfo=split(/\_/,$id);
		my $chr=join("_",@posinfo[0..$#posinfo-1]);
		my $pos=$posinfo[-1];
		if ($pos - $stepsize < 0) {
			$pos{$chr}{0}{$id}=$pos;
		}elsif ($pos - $winsize < 0) {
			my $end=int($pos/$stepsize)+1;
			for (my $i=0;$i<$end;$i++) {
				$pos{$chr}{$i}{$id}=$pos;
			}
		}else{
			my $a=$pos-$winsize;
			$a =0 if ($a < 0);
			my $start=int(($a+$stepsize)/$stepsize);
			my $end=int($pos/$stepsize)+1;
			for (my $i=$start;$i<$end;$i++) {
				$pos{$chr}{$i}{$id}=$pos;
			}
		}
		for (my $i=0;$i<@info;$i++) {
			my ($geno,$phase)=split(/\s+/,$info[$i]);
			$Phase{$Head[$i]}{$id}=$phase;
		}
	}
}
close In;
open Out,">$dir/$Key.sort.markertype";
print Out "#MarkerID\tchr\tpos\twinID\t",join("\t",@Head),"\n";
foreach my $chr (sort keys %pos) {
	foreach my $winID (sort {$a<=>$b} keys %{$pos{$chr}}) {
		my @id=sort{$pos{$chr}{$winID}{$a}<=>$pos{$chr}{$winID}{$b}}keys %{$pos{$chr}{$winID}};
		foreach my $id (@id) {
			print Out $id,"\t",$chr,"\t",$pos{$chr}{$winID}{$id},"\t",$winID,"\t";
			my @out;
			for (my $i=0;$i<@Head;$i++) {
				if (!exists $Phase{$Head[$i]}{$id}) {
					next;
				}
				push @out,$Phase{$Head[$i]}{$id};
			}
			print Out join("\t",@out),"\n";
		}
	}
}
close Out;
open Pri,">$dir/$Key.wintype";
print Pri "#chr\trange\tNo.SNP\t",join("\t",@Head),"\n";
my %wintype;
foreach my $chr (sort keys %pos) {
	if (scalar keys %{$pos{$chr}} == 1) { #һ��scaffold��ֻ��һ��windows�����,�޷������ݣ�ȱʧ����
		my $winID=(sort {$a<=>$b} keys %{$pos{$chr}})[0];
		my @id=sort{$pos{$chr}{$winID}{$a}<=>$pos{$chr}{$winID}{$b}}keys %{$pos{$chr}{$winID}};
		my @out;
		for (my $i=0;$i<@Head;$i++) {
			my $sumPhase=0;
			for (my $j=0;$j<@id;$j++) {
				$sumPhase+=$Phase{$Head[$i]}{$id[$j]};
			}
			$sumPhase =1 if ($sumPhase > 0);
			$sumPhase =-1 if ($sumPhase <0);
			push @out,$sumPhase; 
		}
		print Pri $chr,"\t",0,"-",$pos{$chr}{$winID}{$id[-1]},"\t",scalar @id,"\t",join("\t",@out),"\n";
	}else{
		my @winID=sort {$a<=>$b}keys %{$pos{$chr}};
		for (my $i=0;$i<@winID;$i++) {
			my @id=sort{$pos{$chr}{$winID[$i]}{$a}<=>$pos{$chr}{$winID[$i]}{$b}}keys %{$pos{$chr}{$winID[$i]}};
			my $sumPhase=0;
			my @out;
			for (my $j=0;$j<@Head;$j++) {
				my $sumPhase=0;
				for (my $k=0;$k<@id;$k++) {
					$sumPhase+=$Phase{$Head[$j]}{$id[$k]};
				}
				$sumPhase =1 if ($sumPhase > 0);
				$sumPhase =-1 if ($sumPhase <0);
				push @out,$sumPhase; 
				push @{$wintype{$chr}{$winID[$i]}},$sumPhase;
			}
			print Pri $chr,"\t",$winID[$i]*$stepsize,"-",$winID[$i]*$stepsize+$winsize,"\t",scalar @id,"\t",join("\t",@out),"\n";
		}
	}
}
close Pri;
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub USAGE {#
        my $usage=<<"USAGE";
Contact:        long.huang\@majorbio.com;
Script:			$Script
Description: wintype
	eg:
	perl $Script -i -o 

Usage:
  Options:
  Options:
  -i	<file>	input genotype file name
			#MakrerID\tTYPE\tSample
  -o	<dir>	output dir
  -k	<str>	output file name 
  -win	<num>	input window size(kb)	default 500
  -step	<num>	input step size(kb)	default w/5
  -h         Help
  

USAGE
        print $usage;
        exit;
}

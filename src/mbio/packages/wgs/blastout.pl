#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fIn1,$fIn2,$key,$fOut);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$fIn1,
	"k:s"=>\$fIn2,
	"key:s"=>\$key,
	"o:s"=>\$fOut
			) or &USAGE;
&USAGE unless ($fIn1 and $fIn2);
open In1,$fIn1;
open In2,$fIn2;
open Out1,">$fOut.result-1.out";
open Out2,">$fOut.result-2.out";
my %hash;
while(<In1>){
	my @line=split(/\t/);
	if (exists $hash{$line[0]}){
		my @long=split(/\t/,$hash{$line[0]});
		if ($line[3]> $long[3]){
			$hash{$line[0]}=$_;
		}
	}
	else{
		$hash{$line[0]}=$_;
	}
}

while(<In2>){
	next if (/^$/ or "");
	my @line=split(/\s+/);
	my @str=split(/\s+/,$_,2);
	next if (scalar(@line) < 6);
	if (exists $hash{$line[0]}){
		my @long=split(/\t/,$hash{$line[0]});
		if ($line[3]>= $long[6] and $line[4]<=$long[7]){
			my $dis=int($line[4]-$long[6]);
			if ($long[8] > $long[9]){
				my $pos2=int($long[9]+$dis);
				print Out1"$long[1]\t$str[-1]";
			}
			else{
				my $pos2=int($long[8]+$dis);
				print Out1"$long[1]\t$str[-1]";
			}
		}
		else{
			print Out2"$_";
		}
	}
	else{
		print Out2"$_";
	}
}
close In1;
close In2;
close Out1;
close Out2;
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
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
  -i	<file>	blast out file
  -k	<file>	newmisa file
# -key	<key>	key of blastdb(ref or seq)
  -o	<file>	out file
  -h         Help

USAGE
        print $usage;
        exit;
}

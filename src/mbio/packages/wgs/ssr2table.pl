#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fIn,$fOut);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$fIn,
	"o:s"=>\$fOut
			) or &USAGE;
&USAGE unless ($fIn );
open In,$fIn;
open Out,">",$fOut;
my (%hash,%hash1);
my @ssrtype=("c","c*","p1","p2","p3","p4","p5","p6");
while (<In>){
	next if (/^#/ or /^$/ or "");
	my @line=split(/\t/);
	$hash{$line[0]}{$line[2]}+=1;
	$hash1{$line[0]}+=1;
}
close In;
print Out"#Chr\tSSR Number\tc\tc*\tp1\tp2\tp3\tp4\tp5\tp6\n";
foreach my $key1 (sort {$a cmp $b} keys %hash){
	print Out"$key1\t$hash1{$key1}";
	foreach my $tp (@ssrtype){
		if (exists $hash{$key1}{$tp}){
			print Out"\t$hash{$key1}{$tp}";
		}
		else{
			print Out"\t0";
		}
	}
	print Out"\n";
}
close Out;
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
  -i	<file>	input file name
  -o	<file>	split windows sh
  -h         Help

USAGE
        print $usage;
        exit;
}

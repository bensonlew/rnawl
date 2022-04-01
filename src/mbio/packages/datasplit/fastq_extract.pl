#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fq1,$fq2,$num,$fOut);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"1:s"=>\$fq1,
	"2:s"=>\$fq2,
	"n:s"=>\$num,
	"o:s"=>\$fOut,
			) or &USAGE;
&USAGE unless ($fq1 and $fOut);
open Out,">$fOut";
if ($fq1 =~ /gz/){
	open Fq1,"gzip -dc $fq1|" || die "error open $fq1";
}else{
	open Fq1,$fq1 || die "error open $fq1";
}
my $Sname=basename($fq1);
$Sname=~s/.1.fq//g;
$Sname=(split("-",$Sname))[-1];
my $seqnum=1;
my $read;
my $seq;
my $plus;
my $qual;
while ($read=<Fq1>) {
	$seq=<Fq1>;
	$plus = <Fq1>;
	$qual = <Fq1>;
	last if ($seqnum > $num);
	if (rand() < 0.01) {
		print Out ">$Sname"."1\_"."$seqnum\n$seq";
		$seqnum++;
	}
}
close Fq1;
if (-f $fq2){
	if ($fq2=~/gz/){
		open Fq2,"gzip -dc $fq2|" || die "error open $fq2";
	}else{
		open Fq2,$fq2 || die "error open $fq2";
	}
	$Sname=basename($fq1);
	$Sname=~s/.1.fq//g;
	$Sname=(split("-",$Sname))[-1];
	$seqnum=1;
	while ($read=<Fq2>) {
		$seq=<Fq2>;
		$plus = <Fq2>;
		$qual = <Fq2>;
		last if ($seqnum > $num);
		if (rand() < 0.01) {
			print Out ">$Sname"."2\_"."$seqnum\n$seq";
			$seqnum++;
		}
		print $seqnum,"\n";
	}
	close Fq2;

	close Out;
}else{
}



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
  -1	<file>	input fq1 name
  -2	<file>	input fq2 name
  -o	<file>	output file name
  -n	<num>	input seq num to extract
  -h         Help

USAGE
        print $usage;
        exit;
}

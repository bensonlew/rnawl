#!/usr/bin/env perl
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($bamstat,$metric,$fOut,$insert,$depth,$Key,$dictfile);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"b:s"=>\$bamstat,
	"d:s"=>\$dictfile,
	"m:s"=>\$metric,
	"i:s"=>\$insert,
	"c:s"=>\$depth,
	"o:s"=>\$fOut,
	"k:s"=>\$Key,
			) or &USAGE;
&USAGE unless ($fOut and $bamstat and $insert and $depth and $dictfile);
open In,$bamstat;
open IS,">$insert";
open COV,">$depth";
my %mapstat;
my %cov;
while (<In>) {
	chomp;
	next if ($_ eq "" || /^$/ || /^#/ );
	if (/^SN/) {
		if (/raw total sequences:/) {
			$mapstat{total}=(split(/\t/,$_))[2];
		}
		if (/reads mapped\:/) {
			$mapstat{mapped}=(split(/\t/,$_))[2];
		}
		if (/reads mapped and paired:/) {
			$mapstat{paired}=(split(/\t/,$_))[2];
		}
		if (/reads properly paired:/) {
			$mapstat{properly}=(split(/\t/,$_))[2];
		}
		if (/insert size average:/) {
			$mapstat{insert}=(split(/\t/,$_))[2];
		}
		if (/total length:/) {
			$mapstat{totallen}=(split(/\t/,$_))[2];
		}
		if (/bases mapped:/) {
			$mapstat{coverage}=(split(/\t/,$_))[2];
		}
	}
	if (/^IS/) {
		my (undef,$insert,$depth,undef,undef)=split("\t",$_);
		print IS $insert,"\t",$depth,"\n";
	}
	if (/^COV/) {
		my (undef,undef,$deps,$cova)=split("\t",$_);
		$cov{1}+=$cova if($deps >= 1);
		$cov{5}+=$cova if($deps >= 5);
		print COV $deps,"\t",$cova,"\n";
	}
}
close IS;
close COV;
close In;
if ($metric ne "NULL") {
	open In, $metric;
	while (<In>) {
		chomp;
		next if ($_ eq "" || /^$/ || /^#/);
		my @info=split(/\t/,$_);
		$mapstat{dup}=$info[-2] if (scalar @info > 3);
	}
	close In;
}else{
	$mapstat{dup}="-";
}
my $reflen=0;
open In,$dictfile;
while (<In>) {
	chomp;
	next if ($_ eq ""|| /^$/);
	if (/LN:(\d+)/) {
		$reflen+=$1;
	}
}
close In;
open Out,">$fOut";
print Out "#type\t$Key\n";
print Out "mapped ratio(%)\t",sprintf("%.2f",$mapstat{mapped}/$mapstat{total}*100),"\n";
print Out "proper ratio(%)\t",sprintf("%.2f",$mapstat{properly}/$mapstat{total}*100),"\n";
print Out "duplicate ratio(%)\t",sprintf("%.2f",$mapstat{dup}*100),"\n";
print Out "average insert size\t",$mapstat{insert},"\n";
#print Out "genome coverage\t",sprintf("%.2f",$mapstat{mapbase}/$reflen),"\n";
print Out "average depth\t",sprintf("%.2f",$mapstat{coverage}/$reflen),"\n";
print Out "real depth\t",sprintf("%.2f",$mapstat{coverage}/$cov{1}),"\n";
print Out "cover base\t",$cov{1},"\n";
print Out "genome coverage(1X)\t",sprintf("%.2f",$cov{1}/$reflen*100),"\n";
print Out "genome coverage(5X)\t",sprintf("%.2f",$cov{5}/$reflen*100),"\n";
close Out;
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub ABSOLUTE_DIR #$pavfile=&ABSOLUTE_DIR($pavfile);
{
	my $cur_dir=`pwd`;chomp($cur_dir);
	my ($in)=@_;
	my $return="";
	if(-f $in){
		my $dir=dirname($in);
		my $file=basename($in);
		chdir $dir;$dir=`pwd`;chomp $dir;
		$return="$dir/$file";
	}elsif(-d $in){
		chdir $in;$return=`pwd`;chomp $return;
	}else{
		warn "Warning just for file and dir \n$in";
		exit;
	}
	chdir $cur_dir;
	return $return;
}

sub USAGE {#
        my $usage=<<"USAGE";
Contact:        long.huang\@majorbio.com;
Script:			$Script
Description:
	fq thanslate to fa format
	eg:
	perl $Script -i -o -k -c

Usage:
  Options:
  -b	<file>	samtools flagstat file
  -m	<file>	picard metric file
  -r	<file>	picard dict file
  -i	<file>	output insert file
  -d	<file>	output cover depth file
  -o	<file>	output stat file
  -k	<srt>	sample name
  -h         Help

USAGE
        print $usage;
        exit;
}

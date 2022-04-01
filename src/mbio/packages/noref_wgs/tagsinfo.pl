#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
my ($in,$out,$dep);        #1
GetOptions(
	"help|?" =>\&USAGE,
	"in:s"=>\$in,     #2
	"out:s"=>\$out,   #3
	"dep:s"=>\$dep,   #3
			) or &USAGE;
&USAGE unless ($in and $out); #4
#######################################################################################
open OUT,">$out";
print OUT "#Sample ID\tTags Number\tAverage Length\tAverage Depth\tTag Coverage(5X)(%)\tTag Coverage(10X)(%)\n";
open DEP,">$dep";
print DEP "#Sample ID\tDepth\tTag Num\n";
my @files = glob "$in/*.clustS.gz";
foreach my $file (@files){
	open IN,"gzip -dc $file|" or die "can't open $in\n";
	$/ = "//\n//";
	my $tagnum = 0;
	my $taglength = 0;
	my $tagdepth = 0;
	my $tagnum5x = 0;
	my $tagnum10x = 0;
	my $temp = basename($file);
	$temp =~/^(.*)\.clustS\.gz$/;
	my $name = $1;
	while (<IN>){
		chomp;
		next if (/^\s*$/);
		$tagnum++;
		my $count = 0;
		my $length = 0;
		my $depth = 0;
		my @lines = split /\n/;
		foreach my $line (@lines){
			$count++;
			if ($line=~/size=([0-9]+)/){
				$depth+=$1;
			}else{
				$length+=length $line;
			}
		}
		$taglength+=($length/$count);
		print DEP join("\t",$name,$depth),"\t1\n";
		$tagdepth+=$depth;
		if ($depth >=5){
			$tagnum5x++;
		}
		if ($depth >=10){
			$tagnum10x++;
		}
	}
	close IN;
	# my $temp = basename($file);
	# $temp =~/^(.*)\.clustS\.gz$/;
	# my $name = $1;
	print OUT join("\t",$name,$tagnum,sprintf("%.2f",$taglength/$tagnum),sprintf("%.2f",$tagdepth/$tagnum),sprintf("%.2f",$tagnum5x/$tagnum*100),sprintf("%.2f",$tagnum10x/$tagnum*100) ),"\n";
}
close OUT;
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
		warn "Warning! just for existed file and dir \n$in";
		exit;
	}
	chdir $cur_dir;
	return $return;
}

sub USAGE {           #5
        my $usage=<<"USAGE";
Contact:	tong.wang\@majorbio.com
Version:	$version
Script:		$Script
Description:
		for ipyrad result, calculate tag clust info
Usage:
  Options:
  -in	<dir>	input dir, e.g. ./data_clust_0.85/
  -out	<file>	output file
	-dep	<file>	depth file
  -h		Help

USAGE
        print $usage;
        exit;
}

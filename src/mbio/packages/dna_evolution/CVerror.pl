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
	"o:s"=>\$fOut,
			) or &USAGE;
&USAGE unless ($fIn and $fOut);
open In,$fIn;
open Out,">$fOut/cv.error";
my %CV;
my $bestfile;
my $min=1;
my $K=1;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/);
	my ($i,$log,$xls)=split(/\s+/,$_);
	next if ($i == 1);
	my $CV=`grep \"CV error\" $log`;
	chomp $CV;
	$CV=(split(/\s+/,$CV))[-1];
	print Out "$i\t$CV\n";
	if ($CV < $min) {
		$min=$CV;
		$bestfile=$xls;
		$K=$i;
	}
}
close In;
close Out;
`cp $bestfile $fOut/best.$K.xls`;
# my $job="Rscript $Bin/cverror.R --infile $fOut/cv.error --outfile $fOut/cv.error";
# print $job;
# `$job`;
open In,$bestfile;
open Out,">$fOut/group.list";
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/);
	my ($id,@info)=split(/\s+/,$_);
	my $max=0;
	my $group=0;
	for (my $i=0;$i<@info;$i++) {
		if ($info[$i] > $max) {
			$max=$info[$i];
			$group=$i+1;
		}
	}
	print Out $id,"\t",$group,"\n";
}
close In;
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
  -i	<file>	input file name
  -o	<file>	split windows sh
  -h         Help

USAGE
        print $usage;
        exit;
}

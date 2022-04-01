#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($tags,$fOut,$sampleID);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$tags,
	"o:s"=>\$fOut,
	"k:s"=>\$sampleID,
			) or &USAGE;
&USAGE unless ($tags and $fOut and $sampleID );
open In,"zcat $tags.tags.tsv.gz|";
my $pid="";
my $dep=0;
my $dep5=0;
my $dep10=0;
my $total=0;
my $all_dep=0;
my $len;
while (<In>) {
	chomp;
	next if ($_ eq "" || /^$/ || /^#/);
	my (undef,$tid,$type,undef,undef,$tag,undef)=split(/\t/,$_);
	if ($tid ne $pid) {
			$dep5++ if ($dep >=5);
			$dep10++ if($dep >=10);
			$all_dep+=$dep;
			$dep=0;
			$len+=length($tag);
			$total++;
			$pid=$tid;
	}
	if (/primary/ || /secondary/) {
		$dep++;
	}
}
close In;
open Out,">$fOut";
print Out "#sampleID\ttotal tags\ttotal depth\taverage depth\taverage length\tdep5\tdep10\n";
print Out "$sampleID\t$total\t$all_dep\t",sprintf("%.2f",$all_dep/$total),"\t",sprintf("%.2f",$len/$total),"\t",$dep5,"\t",$dep10,"\n";
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
  -o	<file>	output file
  -h         Help

USAGE
        print $usage;
        exit;
}

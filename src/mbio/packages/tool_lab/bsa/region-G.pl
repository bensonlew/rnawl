#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($input,$output,$bootstrap,$pvalue,$nsites);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"input:s"=>\$input,
	"pvalue:s"=>\$pvalue,
	"output:s"=>\$output,
	"nsites:s"=>\$nsites
			) or &USAGE;
&USAGE unless ($input and $output);
$nsites||=10;
open In,$input;
open G,">$output.G";
my @head;
my $chr="";
my $startI="";
my $endI="";
my $numI=0;
my $startE="";
my $endE="";
my $numE=0;
my $startG="";
my $endG="";
my $numG=0;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/);
	if (/X\.chr/) {
		@head=split;
	}else{
		my @info=split;
		my %info;
		for (my $i=0;$i<@head;$i++) {
			$info{$head[$i]}=$info[$i];
		}
		if($info{"X.chr"} ne $chr){
			print G  join("\t",$chr,$startG,$endG,$numG),"\n"   if ($numG !=0);
			$chr=$info{"X.chr"};
			$startI="";
			$startG="";
			$startE="";
			$endI="";
			$endG="";
			$endE="";
			$numI=0;
			$numG=0;
			$numE=0;
		}
		if ($info{qvalue} > $pvalue  ) {
			if ($numG !=0) {
				print G join("\t",$chr,$startG,$endG,$numG),"\n";
			}
			$chr=$info{"X.chr"};
			$startG="";
			$endG="";
			$numG=0;
		}else{
				$chr=$info{"X.chr"};
				$startG=$info{pos} if ($numG == 0);
				$endG=$info{pos};
				$numG++;

		}
	}
}
print G  join("\t",$chr,$startG,$endG,$numG),"\n"   if ($numG !=0);
close G;
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
Usage:
  Options:
  -input	<file>	input reference file name
  -output	<file>	input gff file name
	"bootstrap:s"=>\$bootstrap,
	"pvalue:s"=>\$pvalue,
	"output:s"=>\$output,
 
-h         Help

USAGE
        print $usage;
        exit;
}

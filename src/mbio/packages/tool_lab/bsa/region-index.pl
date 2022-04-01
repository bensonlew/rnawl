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
	"bootstrap:s"=>\$bootstrap,
	"pvalue:s"=>\$pvalue,
	"output:s"=>\$output,
	"nsites:s"=>\$nsites
			) or &USAGE;
&USAGE unless ($input and $output and $bootstrap);
$nsites||=10;
(open In,$bootstrap)||die;;
my %index;
my @head1;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/);
	if (/depth/) {
		(undef,@head1)=split;
	}else{
		my ($depth,@info)=split(/\s+/,$_);
		for (my $i=0;$i<@info;$i++) {
			$index{$depth}{$head1[$i]}=$info[$i];
		}
	}
}
close In;
open In,$input;
open Index,">$output.index";
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
print(1-$pvalue);
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
		$info{CI}=$index{$info{mdepth}}{$pvalue};
		print $info{slidingD},"\t",$info{CI},"\n";
		$info{pos}=($info{pos1}+$info{pos2})/2;
		if($info{"X.chr"} ne $chr){
			print Index join("\t",$chr,$startI,$endI,$numI),"\n"if ($numI !=0);
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
		if ($info{slidingD} < $info{CI} ) {
			if ($numI !=0) {
				if($chr eq "chr14"){
					print $_,"\n\n";
					print join("\t",$chr,$startI,$endI,$numI),"\n\n";
					print join("\t",$chr,$startI,$info{mdepth},$info{slidingD}),"\t",$info{CI},"\n" ;
					die;
				}
				print Index join("\t",$chr,$startI,$endI,$numI),"\n";
			}
			$chr=$info{"X.chr"};
			$startI="";
			$endI="";
			$numI=0;
		}else {
				$chr=$info{"X.chr"};
				$startI=$info{pos} if ($numI == 0);
				$endI=$info{pos};
				$numI+=$info{ncount};
		}
	}
}
print Index join("\t",$chr,$startI,$endI,$numI),"\n"if ($numI !=0);
close Index;
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

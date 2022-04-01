#!/usr/bin/env perl
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
$/="#";
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/);
	next if (!/Count by effects/ && !/Count by genomic region/);
	if (/Count by effects/) {
		open Out,">$fOut.effects";
		my @line=split(/\n/,$_);
		foreach my $l (@line) {
			my @info=split(/ \, /,$l);
			next if (scalar @info != 3);
			print Out join("\t",@info),"\n";
		}
		close Out;
	}
	if (/Count by genomic region/) {
		open Out,">$fOut.region";
		my @line=split(/\n/,$_);
		foreach my $l (@line) {
			my @info=split(/ \, /,$l);
			next if (scalar @info != 3);
			print Out join("\t",@info),"\n";
		}
		close Out;
	}

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

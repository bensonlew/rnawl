#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fIn,$fOut,$select);
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
open Out,">$fOut";
print Out "#marker\tchr\tpos\tlod\tvar\tpm1\tpm2\tstart\tend\tmark1\tmark2\n";
my %stat;
while (<In>) {
	chomp;
	s/\"//g;
	next if ($_=~ ""|| /^$/ || /#/|| /marker/);
	my ($marker,$chr,$pos,$lod,$var,$pm1,$pm2,$start,$end,$mark1,$mark2)=split(/\s+/,$_);
	my ($chr1,$start1)=split(/\_/,$mark1);
	my ($chr2,$start2)=split(/\_/,$mark2);
	if ($chr1 ne $chr2) {
		die "error region!";
	}else{
		#my $dis=$start2 - $start1;
		#if($dis < 5000){
		#	$start1=$start1 - 4000 ;
		#	$start2=$start2 + 4000 ;
		#}
		#my $region=join("\t",$start1,$start2);
		$chr=$chr1;
		print Out "$marker\t$chr\t$pos\t$lod\t$var\t$pm1\t$pm2\t$start\t$end\t$start1\t$start2\n";
	}
}
close Out;
close In;

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
  -o	<file>	output file name
  -h         Help

USAGE
        print $usage;
        exit;
}

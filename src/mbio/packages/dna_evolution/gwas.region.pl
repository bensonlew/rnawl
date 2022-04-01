#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($input,$output,$threshold,$adjust,$distance,$Chr);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"input:s"=>\$input,
	"output:s"=>\$output,
	"distance:s"=>\$distance,
	"threshold:s"=>\$threshold,
	"adjust:s"=>\$adjust,
			) or &USAGE;
&USAGE unless ($input and $output and $threshold);
$distance||=200000;
$threshold||=0.05;
$adjust||=1;
$Chr||="chr";
my $nsnp=`wc -l $input`;
chomp $nsnp;
$nsnp=(split(/\s+/,$nsnp))[0];
$nsnp=$nsnp-1;
if ($adjust == 1) {
	$threshold=$threshold/$nsnp;
}
open In,$input;
my %region;
while (<In>) {
	s/\"//g;
	s/\s+//g;
	chomp;
	next if ($_ eq ""||/^$/||/Chrom/);
	my ($chr,$CHR,$pos,$eff,$pvalue)=split(/\,/,$_);
	if ($pvalue eq "NA") {
		next;
	}
	if ($pvalue < $threshold) {
		my $start=$pos-$distance;
		$start=1 if ($start <= 0);
		my $end=$pos+$distance;
		my $regioned=0;
		
		foreach my $region (sort keys %{$region{$chr}}) {
			my ($chr1,$start1,$end1)=split(/\s+/,$region);
			if (($start1 > $end && $end1 > $end) || ($end1 < $start && $start1 < $start)) {
		
			}else{
				my @region=sort {$a<=>$b} ($start,$end,$start1,$end1);
				delete $region{$chr}{$region};
				$region{$chr}{join("\t",$chr,$region[0],$region[-1])}=1;
				$regioned=1;
			}
		}
		if ($regioned == 0) {
			$region{$chr}{join("\t",$chr,$start,$end)}=1;
		}
	}
}
close In;
#print Dumper \%region;die;
open Out,">$output";
print Out "#chr\tstart\tend\n";
foreach my $chr(sort keys %region) {
	foreach my $region (sort keys %{$region{$chr}}) {
		print Out $region,"\n";
	}
}
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
Usage:
  Options:
  -input	<file>	input reference file name
  -output	<file>	input gff file name
  -distance     region size (default 200000)
  -threshold	threshold (default 0.05)
  -adjust       default 1
  }
  -h         Help

USAGE
        print $usage;
        exit;
}

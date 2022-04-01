#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fIn,$fOut);
use Data::Dumper;
use FindBin qw($Bin $Script);
use SVG;
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$fIn,
	"o:s"=>\$fOut,
			) or &USAGE;
&USAGE unless ($fIn and $fOut);
open In,$fIn;
my $id;
my %group;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/);
	if (/group/) {
		$id=(split(/\s+/,$_))[-1];
	}else{
		my ($markerID,$pos)=split(/\s+/,$_);
		$group{$id}{$pos}++;
	}
}
close In;
open Out,">$fOut";
print  Out "#LGID\tNumber Marker\tNumber Uniq\tTotal Distance\tAvarage Distance\tGap < 5cM(%)\tMax Gap\n";
foreach my $lgid (sort{$a<=>$b} keys %group) {
	my @pos=sort {$a<=>$b} keys %{$group{$lgid}};
	my $distance=$pos[-1];
	my $nuniq=scalar @pos;
	my $nloc=0;
	my $gap=scalar @pos -1;
	my $gap5=0;
	my $maxGap=0;
	for (my $i=0;$i<@pos;$i++) {
		$nloc+=$group{$lgid}{$pos[$i]};
		if ($i > 0 && abs($pos[$i]-$pos[$i-1]) >5) {
			$gap5++;
		}
		if ($i > 0 && abs($pos[$i]-$pos[$i-1]) > $maxGap) {
			$maxGap = abs($pos[$i]-$pos[$i-1]);
		}
	}
	my $gapratio=sprintf("%.2f",100*($nloc - $gap5)/$nloc);
	#print Out join("\t",$lgid,$nloc,$nuniq,$distance,sprintf("%.2f",$distance/($nloc-1)),$gap5,$maxGap),"\n";
	print Out join("\t",$lgid,$nloc,$nuniq,$distance,sprintf("%.2f",$distance/($nloc-1)),$gapratio,$maxGap),"\n";
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

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
			) or &USAGE;
&USAGE unless ($tags and $fOut  );
open In,"$tags";
open Out,">$fOut";
print Out "#sampleID\ttotal tags\ttotal depth\taverage depth\taverage length\tdep5\tdep10\n";
while (<In>) {
	chomp;
	next if ($_ eq "" || /^$/ || /^#/);
	my ($sample,$ustacks)=split(/\s+/,$_);
	open STAT,"$ustacks.tags.stat";
	while (<STAT>) {
		next if ($_ eq "" || /^$/ || /^#/);
		print Out $_;
	}
	close STAT;
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
  -o	<file>	output file
  -h         Help

USAGE
        print $usage;
        exit;
}

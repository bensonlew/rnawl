#!/usr/bin/perl -w
use 5.10.0;
use strict;
use warnings;
use Getopt::Long;
my ($BEGIN_TIME,$fIn,$fOut,$fMet,$Key);
$BEGIN_TIME=time();
use Cwd;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"r|result:s"=>\$fIn,
	"s|sort:s"=>\$fOut,
			); 		#;/or &USAGE;
&USAGE unless ($fIn);
$fIn=DIR($fIn);
open In,$fIn;
open Out,">$fOut";
my %ha;
while (<In>) {
	chomp;
	next if ($_ eq "" ||/^$/);
	if(/^#/){print Out $_,"\n";
	}else{
		my($chr,$pos,@ssr)=split/\t/,$_;
		$chr=~/([a-z]*)(\d+)/;
		my $ch=$1;my $r=$2;
		$ha{$ch}{$r}{$pos}=$_;
	}

}
close In;
foreach my $ch (sort keys %ha){
	foreach my $r (sort {$a<=>$b} keys %{$ha{$ch}}){
		foreach my $pos (sort {$a<=>$b} keys %{$ha{$ch}{$r}}){
			print Out $ha{$ch}{$r}{$pos},"\n";
		}
	}
	
}
close Out;
#######################################################################################
# print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub DIR{
	my $cur_dir=`pwd`;
	chomp($cur_dir);
	my ($in)=@_;
	my $return="";
	if(-f $in){
		my $dir=dirname($in);
		my $file=basename($in);
		chdir $dir;
		$dir=`pwd`;
		chomp $dir;
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
sub USAGE {
        my $usage=<<"USAGE";
Contact:        qingmei.cui\@majorbio.com;
Script:			$Script
Description:
	fq thanslate to fa format
	eg:
	perl $Script -i -o
Usage:
  Options:
  -i	<file>	input file name
  -o	<file>	output file name
  -h         Help

USAGE
        print $usage;
        exit;
}

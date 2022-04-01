#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fin,$fout,$adjust,$threshold,$Chr,$dsh,$list);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"input:s"=>\$fin,
	"output:s"=>\$fout,
	"adjust:s"=>\$adjust,
	"threshold:s"=>\$threshold,
			) or &USAGE;
&USAGE unless ($fout);
$threshold||=0.05;
$adjust||=1;
my $nsnp=`wc -l $fin`;
chomp $nsnp;
$nsnp=(split(/\s+/,$nsnp))[0];
$nsnp=$nsnp-1;
if ($adjust == 1) {
	$threshold = $threshold/$nsnp;
}
open IN,$fin;
open OUT,">$fout";
while (<IN>) {
	chomp;
	s/\"//g;
	s/ //g;
	next if ($_ eq ""||/^$/||/Chrom/);
	my ($chr,$CHR,$pos,$eff,$pvalue)=split(/\,/,$_);
	next if ($pvalue eq "NA");
	next if ($pvalue > $threshold);
	print OUT "$chr\_$pos\n";
}
close IN;
close OUT;
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
Contact:        minghao.zhang\@majorbio.com;
Script:			$Script
Description:
	fq thanslate to fa format
	eg:
	perl $Script -i -o -k -c

Usage:
  Options:
	"input:s"=>\$fin,  step02  gwas.dir.list
	"output:s"=>\$fout,
	"adjust:s"=>\$adjust,
	"threshold:s"=>\$threshold,
	"chr:s"=>\$Chr,
	"dsh:s"=>\$dsh 
  -h         Help

USAGE
        print $usage;
        exit;
}

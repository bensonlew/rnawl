#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($trt,$out);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"trt:s"=>\$trt,
	"out:s"=>\$out,
			) or &USAGE;
&USAGE unless ($trt);
sub USAGE {#
        my $usage=<<"USAGE";
Contact:        nanshan.yang\@majorbio.com;
Script:			$Script
Description:
	
	eg:
	perl $Script -i -o -k -c

Usage:
  Options:
  -trt	<file>	input file name
  -out	<dir>	
  -h         Help

USAGE
        print $usage;
        exit;
}
open IN,$trt;
my @trits;
my @trits_tas;
my $line=0;
while(<IN>){
    $_=~s/[\n\r]//g;
    if($line==0){
         my ($sampleid,@trit)=split/\t/,$_;
         push(@trits,@trit);
        }
        $line++;
}
close IN;

for(my $i=0;$i<=$#trits;$i++){
    my $num=$i+2;
    `cat $trt|cut -f 1,$num > $out/$trits[$i].trt`;
    `cat $trt|cut -f 1,$num|sed "s/SampleID/<Trait>/g" > $out/$trits[$i].tassel.trt`;
    }

#`cat $trt|sed "s/SampleID/<Trait>/g" > $out/trit.TASSEL.trt`;

open OUT1,">$out/trit.list";
open OUT2,">$out/trit.tassel.list";
for(@trits){
    print OUT1 "$_\t$out/$_.trt\n";
    print OUT2 "$_\t$out/$_.tassel.trt\n";
    }
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub ABSOLUTE_DIR
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


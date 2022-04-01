#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($in,$out);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$in,
	"o:s"=>\$out,
			) or &USAGE;
&USAGE unless ($in and $out);

#misa
#ID	SSR nr.	SSR type	SSR	size	start	end
#chr1	1	p2	(TA)7	14	27213	27226
open IN,$in;
open Out,">$out";
print Out "ID\tSSR nr.\tSSR type\tSSR\tsize\tstart\tend\n";
my ($type,$ssr,$size,$start,$end,$n);
while (<IN>) {
	chomp;
	next if ($_ eq ""|| /^$/|| /^#/);
	my($chr,$pos,$ssrbase,$unit,$dep,$total)=split(/\t/,$_);#chr	pos	SSRbase	Unit	Dep	Totaldep chr1	27213	TA	4:5:8	3:1:3	7
	next if($total < "10");
	if($unit=~/\:/){
		my@unit=split(/\:/,$unit);
		my@dep=split(/\:/,$dep);
		for(my$i=0;$i<scalar@unit;$i++){
			my $nowdep=$dep[$i];
			next if($nowdep < "2");
			$type=length$ssrbase;
			$ssr="($ssrbase)$unit[$i]";
			$size=(length$ssrbase)*$unit[$i];
			next if($size < "8");
			$start=$pos;
			$end=$start + $size -1 ;
			$n++ ;
			print Out "$chr\t$n\tp$type\t$ssr\t$size\t$start\t$end\n";
		}
	}else{
		$type=length$ssrbase;
		$ssr="($ssrbase)$unit";
		$size=$type*$unit;
		next if($size < "8");
		$start=$pos;
		$end=$start + $size -1 ;
		$n++ ;
		print Out "$chr\t$n\tp$type\t$ssr\t$size\t$start\t$end\n";
	}
}
close Out;
close IN;
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
Contact:        caixia.tian\@majorbio.com;
Script:			$Script
Description:
	fq thanslate to fa format
	eg:
	perl $Script -i -o -k -c

Usage:
  Options:
  -i	<file>	input ssr.result file
  -o	<dir>	output dir
  -h         Help

USAGE
        print $usage;
        exit;
}

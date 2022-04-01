#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fqlist,$out);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"f|fqlist:s"=>\$fqlist,
	"o|out:s"=>\$out,
			) or &USAGE;
&USAGE unless ($fqlist and $out);
$fqlist=ABSOLUTE_DIR($fqlist);
mkdir $out if (!-d "$out");
open In,$fqlist;
my %sample;
while (<In>) {
	chomp;
	next if ($_ eq ""|| /^$/);
	my ($id,$fq1,$fq2)=split(/\s+/,$_);
	if(exists $sample{$id}){
		open Out,">>$out/$id.config";
		$sample{$id}++;
	}else{
		open Out,">$out/$id.config";
		$sample{$id}=1;
	}
	print Out "max_rd_len=150\n";
	print Out "[LIB]\n";
	print Out "avg_ins=400\n";
	print Out "reverse_seq=0\n";
	print Out "asm_flags=3\n";
	print Out "rank=$sample{$id}\n";
	print Out "q1=$fq1\n";
	print Out "q2=$fq2\n";
	close Out;
}
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
  -fqlist	<file>	input fqlist file
  -out	<dir>	output dir
  -h         Help

USAGE
        print $usage;
        exit;
}

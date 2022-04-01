#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($vcf,$fOut,$catalog);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"vcf:s"=>\$vcf,
	"catalog:s"=>\$catalog,
	"out:s"=>\$fOut,
			) or &USAGE;
&USAGE unless ($vcf and $fOut and $catalog  );
open In,"$vcf";
my %variant;
while (<In>) {
	chomp;
	next if ($_ eq "" || /^$/ || /^#/);
	my ($chr,$pos,undef)=split(/\s+/,$_);
	$variant{$chr}{$pos}=1;
}
close In;
open In,"zcat $catalog|";
open Out,">$fOut";
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/ || /^#/);
	my (undef,$cataid,undef,undef,undef,$seq,undef)=split(/\s+/,$_);
	print Out ">$cataid\n";
	print Out $seq,"\n";
	my @out;
	for (my $i=0;$i<length($seq);$i++) {
		if(!exists $variant{$cataid}{$i+1}){
			push @out," ";
		}else{
			push @out,"*";
		}
	}
	print Out join("",@out),"\n";
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
  -vcf	<file>	input file name
  -catalog	<file>	catalog file
  -out	<file>	output file
  -h         Help

USAGE
        print $usage;
        exit;
}

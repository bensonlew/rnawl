#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($ref,$gff,$gene,$out);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"ref:s"=>\$ref,
	"gff:s"=>\$gff,
	"gene:s"=>\$gene,
	"out:s"=>\$out
			) or &USAGE;
&USAGE unless ($ref and $gff and $gene and $out);
open In,$ref;
$/=">";
my %seq;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/);
	my ($id,@seq)=split(/\n/,$_);
	$seq{$id}=join("",@seq);
}
close In;
$/="\n";
my $check=0;
open In,$gff;
my %out;
my %pos;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/ ||/^#/);
	my ($chr,$source,$type,$start,$end,undef,$flag,undef,$info)=split(/\t/,$_);
	next if ($type eq "region");
	if ($type eq "gene" && $info=~/ID=([^;]*)/) {
		last if ($check == 1);
		if ($1 eq $gene) {
			$check=1;
		}
	}
	next if ($check == 0);
	next if ($type ne "exon");
	my $id;
	if ($info =~ /Parent=([^;]*)/) {
		$id=$1;
		$pos{join("\t",$id,$start,$end)}=join("\t",$chr,$type);
	}
	my $seq=substr($seq{$chr},$start,$end-$start);
	if ($flag eq "-") {
		$seq =~ tr/a-z/A-Z/;
		$seq =~ tr/ATGC/TACG/;
		$seq=reverse($seq);
	}
	$out{$id}.=$seq;
}
close In;
open Out,">$out/$gene".".exon.fa";
foreach my $id (sort keys %out) {
	print Out ">$id\n$out{$id}\n";
}
close Out;
open Out,">$out/$gene".".mirna.xls";	# 画图用
foreach my $id (sort keys %pos) {
	print Out ">$pos{$id}"."\t"."$id\n";
}

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
  -ref	<file>	input reference file name
  -gff	<file>	input gff file name
  -gene	<str>	input geneid
  -out	<file>	output filename
  -h         Help

USAGE
        print $usage;
        exit;
}

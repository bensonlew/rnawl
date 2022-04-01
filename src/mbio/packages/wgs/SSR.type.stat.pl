#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fmisa,$fOut);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"m|misa:s"=>\$fmisa,
	"o:s"=>\$fOut,
			); 		#;/or &USAGE;
&USAGE unless ($fmisa and $fOut);
open In,$fmisa;
my $dir=DIR($fmisa);
my %stat;
my %total;
while (<In>) {
	chomp;
	next if ($_ eq "" || /^ID/ || /^#/);
	my @info=split /\t/;
	$info[0] =~ /([a-z]*)([0-9]*)/i;
	my ($chrtype,$num)=($1,$2);
	$stat{$chrtype}{$num}{$info[2]}++;
	$total{$info[0]}++;
}
close In;
my $chrtype;
my $num;
open Out,">$fOut";
my @header=qw/Chr SSR_Number c c* p1 p2 p3 p4 p5 p6/;
print Out join("\t",@header),"\n";
foreach  $chrtype (sort keys %stat){
	foreach $num (sort {$a<=>$b} keys %{$stat{$chrtype}}){
		my $chr=$chrtype.$num;
		$stat{$chrtype}{$num}{c}||=0;$stat{$chrtype}{$num}{"c*"}||=0;$stat{$chrtype}{$num}{"p1"}||=0;$stat{$chrtype}{$num}{p2}||=0;$stat{$chrtype}{$num}{p3}||=0;$stat{$chrtype}{$num}{p4}||=0;$stat{$chrtype}{$num}{p5}||=0;$stat{$chrtype}{$num}{p6}||=0;
		print Out join("\t",$chr,$total{$chr},$stat{$chrtype}{$num}{c},$stat{$chrtype}{$num}{"c*"},$stat{$chrtype}{$num}{"p1"},$stat{$chrtype}{$num}{p2},$stat{$chrtype}{$num}{p3},$stat{$chrtype}{$num}{p4},$stat{$chrtype}{$num}{p5},$stat{$chrtype}{$num}{p6}),"\n";
	}
}
close Out;
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
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
	eg:
	perl $Script -m *misa  -o ***
Usage:
  Options:
  -i	<file>	input misa name
  -o	<file!>	output file name
  -h         Help

USAGE
        print $usage;
        exit;
}

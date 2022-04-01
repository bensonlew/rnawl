#!/usr/bin/env perl
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fIn,$fOut);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$fIn,
	"o:s"=>\$fOut,
			) or &USAGE;
&USAGE unless ($fIn );
open In,$fIn;
my %stat;
my %length;
my %Alen;
while (<In>) {
	chomp;
	next if ($_ eq "" ||/^$/ ||/^#/);
	my ($chr1,$pos1,$chr2,$pos2,$length,$type,$pvalue,$depth,$genenum,$geneinfo)=split(/\t/,$_);
	$stat{$type}{total}++;
	$stat{$type}{gene}++ if ($genenum !=0);
	$length{$length}{$type}++;
	$Alen{$type}+=$length;
}
close In;
open Out,">$fOut.sv.stat";
print Out "#type\ttotal\tgene\talen\n";
foreach my $type (sort keys %stat) {
	my @out;
	push @out,$type;
	$stat{$type}{total}||=0;
	$stat{$type}{gene}||=0;
	push @out,$stat{$type}{total};
	push @out,$stat{$type}{gene};
	push @out,$Alen{$type}/$stat{$type}{total};
	print Out join("\t",@out),"\n";
}
close Out;
open Out,">$fOut.sv.length";
print Out "#length\t",join("\t",sort keys %stat),"\n";
foreach my $len (sort {$a<=>$b} keys %length) {
	my @out;
	push @out,$len;
	foreach my $type (sort keys %stat) {
		$length{$len}{$type}||=0;
		push @out,$length{$len}{$type};
	}
	print Out join("\t",@out),"\n";
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

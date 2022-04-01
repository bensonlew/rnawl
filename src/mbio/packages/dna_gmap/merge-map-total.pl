#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($list,$map,$fOut);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"list:s"=>\$list,
	"mark:s"=>\$map,
	"out:s"=>\$fOut
			) or &USAGE;
&USAGE unless ($list and $map );
#$dOut||="./";
my %stat;
open In,$list;
while (<In>){
	chomp;
	next if ($_ eq ""|| /^$/);
	my($chr,$len)=split/\t/;
	my $lg=(split(/\D+/,$chr))[-1];
	$stat{$lg}{length}=$len;
}
close In;
my $min=50000000;
my $max=0;
my $maxpos=0;
open In,$map;
$/="\n";
while (<In>){
	chomp;
	next if ($_ eq ""|| /^$/);
	next if ($_=~/group/);
	my $lg=(split(/\D+/,(split(/\_/,(split(/\t/,$_))[0]))[0]))[-1];
	$stat{$lg}{id}=(split(/\t/,$_))[0];
	$stat{$lg}{pos}=(split(/\t/,$_))[1];
	$stat{$lg}{now}=(split(/\_/,$stat{$lg}{id}))[-1];
	$stat{$lg}{max}||=0;
	$stat{$lg}{min}||=0;
	$stat{$lg}{maxpos}||=0;
	if($stat{$lg}{max} < $stat{$lg}{now}){
		$stat{$lg}{max}=$stat{$lg}{now};
	}
	if($stat{$lg}{min} > $stat{$lg}{now}){
		$stat{$lg}{min}=$stat{$lg}{now};
	}
	if($stat{$lg}{maxpos} < $stat{$lg}{pos}){
		$stat{$lg}{maxpos}=$stat{$lg}{pos};
	}
}
close In;
open Out,">$fOut";
print Out"#LG ID\tGenetic Distance(cM)\tPhysical Coverage(%)\tcM/Mb\n";
foreach my $lg (sort {$a<=>$b}keys %stat){
	my $genlen=$stat{$lg}{max} - $stat{$lg}{min};
	my $cover=100*$genlen/$stat{$lg}{length} ;
	print Out join("\t",$lg,$stat{$lg}{maxpos},$cover,sprintf("%.2f",1000000*$stat{$lg}{maxpos}/$stat{$lg}{length})),"\n";
}
close Out;
#######################################################################################
#print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub USAGE {#
        my $usage=<<"USAGE";
Contact:        chongqing.shi\@majorbio.com;
Script:			$Script
Description:
	fq thanslate to fa format
	eg:
	perl $Script -i -o -k -c

Usage:
  Options:
  -list	<file>	input ref.chrlist file
  -mark	<file>  input total.map file 
  -out	<file>	output result file
  -h         Help

USAGE
        print $usage;
        exit;
}

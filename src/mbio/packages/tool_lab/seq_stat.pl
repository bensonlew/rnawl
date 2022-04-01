#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($input,$output,$split);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"input:s"=>\$input,
	"output:s"=>\$output,
	"split:s"=>\$split
			) or &USAGE;
&USAGE unless ($input and $output);

open In,$input;
open Out,">$output.length";
$/=">";
my %length;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/);
	my ($id,@line)=split(/\n/,$_);
	my $base=join("",@line);
	print Out $id,"\t",length($base),"\n";
	my $l=length($base);
	$length{$l}++;
}
close In;
close Out;
my @l=sort {$a<=>$b} keys %length;
my $step=int(($l[-1]-$l[0])/$split*100)/100;
my %stat;
for (my $i=0;$i<@l;$i++) {
	my $win=int($l[$i]/$step);
	my $winname=$win*$step."-".($win+1)*$step;
	$stat{$win}{name}=$winname;
	$stat{$win}{value}+=$length{$l[$i]};
}
open Out,">$output.distribution";
for (my $i=0;$i<$split;$i++) {
	my $winname=$i*$step."-".($i+1)*$step;
	$stat{$i}{value}||=0;
	print Out join("\t",$winname,$stat{$i}{value}),"\n";
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
Usage:
  Options:
  -input	<file>	input reference file name
  -output	<str>	output keys of filename
  -split	<num>	split
  -h         Help

USAGE
        print $usage;
        exit;
}

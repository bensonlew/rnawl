#!/usr/bin/perl -w
use 5.10.0;
use strict;
use warnings;
use Getopt::Long;
my ($BEGIN_TIME,$fIn,$fOut,$fstat);
$BEGIN_TIME=time();
use Cwd;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"seq:s"=>\$fIn,
	"stat:s"=>\$fstat,
	"o:s"=>\$fOut,
			); 		#;/or &USAGE;
&USAGE unless ($fIn );
mkdir $fOut if(!-d $fOut);
$fIn=DIR($fIn);
$fOut=DIR($fOut);
open In,$fIn;
open Out,">$fOut/denovo.stat.xls";
print Out join("\t","#Largest Scaffold","Largest Length","Large Scaffolds(>1000bp)","Bases in Large Scaffolds","N50 Scaffold","N50 Length","N90 Scaffold","N90 Length","GC Content","N Rate"),"\n";
$/=">";
my $max=0;
my ($sum1000,$sum);
my %ha;
while (<In>) {
	chomp;
	next if ($_ eq "" ||/^$/);
	my ($id,@seq)=split /\n/,$_;
	$id=(split /\s+/,$id)[0];
	my $seq=join("",@seq);
	if(length($seq)>$max){$max=length($seq);$ha{1}="$id\t".length($seq);}
	if(length($seq)>=1000){$sum1000+=length($seq);$sum++;}
}
close In;
$/="<--";
open Stat,$fstat;
my ($N,$GC,$N50,$L50,$N90,$L90);
while(<Stat>){
	chomp;
	next if ($_ eq "" ||/^$/ || /Contig/);
	my @info=split/\n/;
	foreach my $a (@info){
		if($a ne ""){
			if($a=~ /GapContent_N\s+([0-9]*)\s+(.*)/){$N=$2;}
			if($a=~ /GC_Content/){$GC=(split/\s+/,$a)[1];}
			if($a=~ /N50\s+([0-9]*)\s+([0-9]*)/){$N50=$1;$L50=$2;}
			if($a=~ /N90\s+([0-9]*)\s+([0-9]*)/){$N90=$1;$L90=$2;}
		}	
	}
	
}
close Stat;
if($N eq "0.00%"){$N="0%";}
print Out join("\t",$ha{1},$sum,$sum1000,$L50,$N50,$L90,$N90,$GC,$N),"\n";
#######################################################################################
# print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
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
	fq thanslate to fa format
	eg:
	perl $Script -i -o
Usage:
  Options:
	"seq:s"=>\$fIn,
	"stat:s"=>\$fstat,
	"o:s"=>\$fOut,

USAGE
        print $usage;
        exit;
}

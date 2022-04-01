#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fIn,$dOut,$Split,$Key);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$fIn,
	"o:s"=>\$dOut,
	"n:s"=>\$Split,
			) or &USAGE;
&USAGE unless ($fIn and $dOut and $Split);
$Split||=50;
open In,$fIn;
mkdir $dOut if (!-d $dOut);
$dOut=ABSOLUTE_DIR($dOut);
my %handsh;
$/=">";
my $nseq=0;
open Falist,">$dOut/fasta.list";
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/);
	$nseq++;
	my $filehand=$nseq % $Split;
	if (!exists $handsh{$filehand}) {
		open $handsh{$filehand},">$dOut/sub.$filehand.fa";
		print Falist "$dOut/sub.$filehand.fa\n";
	}
	my ($id,@line)=split(/\n/,$_);
	$id=(split(/\s+/,$id))[0];
	print {$handsh{$filehand}} ">$id\n",join("\n",@line),"\n";
}
close In;
close Falist;
foreach my $key (sort keys %handsh) {
	close $handsh{$key};
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
	this Script will split fasta file into n files

	eg:
	perl $Script -i demo.fa -d Nr -o ./ -k demo

Usage:
  Options:
  -i	<file>	input fa file name
  -o	<dir>	output dir 
  -n	<num>	split file number
  -h         Help

USAGE
        print $usage;
        exit;
}

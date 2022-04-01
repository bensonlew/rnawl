#!/usr/bin/perl -w
use 5.10.0;
use strict;
use warnings;
use Getopt::Long;
my ($BEGIN_TIME,$fIn,$fOut,$fMet,$Key);
$BEGIN_TIME=time();
use Cwd;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$fIn,
	"o:s"=>\$fOut,
			); 		#;/or &USAGE;
&USAGE unless ($fIn );
open In,$fIn;
open Out,">$fOut";
$fIn=DIR($fIn);
my @header;
my %hash;
while (<In>) {
	chomp;
	next if ($_ eq "" ||/^$/);
	my ($marker,$type,@info)=split /\t/,$_;
	if (/^#/ || /MarkerID/){
		@header = @info;
		print Out "$marker\t$type\t",join("\t",sort @header),"\n";
	}
	else{
		for my $idx (0..$#info){
			$hash{$marker}{$type}{$header[$idx]}=$info[$idx];
		}
	}
}
close In;
# print Dumper \%hash;
foreach my $i (sort keys %hash) {
	foreach my $j (keys %{$hash{$i}}){
		print Out join("\t",$i,$j),"\t";
		my @a;
		foreach my $m (sort keys %{$hash{$i}{$j}}){
			push @a,$hash{$i}{$j}{$m};
		}
		print Out join("\t",@a),"\n";
	}
}
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
	perl $Script -i marker_upload -o 
Usage:
  Options:
  -i	<file>	input file name
  -o	<file>	output file name
  -h         Help

USAGE
        print $usage;
        exit;
}

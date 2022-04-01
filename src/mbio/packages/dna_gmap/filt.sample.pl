#!/usr/bin/perl -w
use 5.10.0;
use strict;
use warnings;
use Getopt::Long;
my ($BEGIN_TIME,$fIn,$fOut,$id);
$BEGIN_TIME=time();
use Cwd;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$fIn,
	"sample:s"=>\$id,	# ,分开
	"o:s"=>\$fOut,
			); 		#;/or &USAGE;
&USAGE unless ($fIn );
open In,$fIn;
$fIn=DIR($fIn);
#####################################存样本
my %samples;
foreach my $sample (split /,/,$id) {
	$samples{$sample} = 1;
}
my @sample_list;
my @sample;
my $num=0;
my %hash;
print Dumper \%samples;
while (<In>) {
	chomp;
	next if ($_ eq "" || /^$/);
	my ($chrid, $type, @info)=split /\t/;
	if (/^#/ || /ID$/) {
		@sample = @info;
		 foreach my $id (@info) {
			next if (!exists $samples{$id});
			push @sample_list,$id;
		}
	}else{
		print join("\t",@sample),"\n";
		$num ++;
		for(my $idx=0;$idx<=@info;$idx++) {
			my $i = $sample[$idx];
			if (exists $samples{$i}){
				$hash{$num}{join("\t",$chrid, $type)}{$sample[$idx]} = $info[$idx+2];
			}
		}
	}
}
open Out,">$fOut";
print Out "#MarkID","\t","Type","\t",join("\t",@sample_list),"\n";
foreach my $nums (sort {$a<=>$b} keys %hash) {
	foreach my $type (keys %{$hash{$nums}}) {
		print Out "$type\t";
		foreach my $sample (keys %{$hash{$nums}{$type}}) {
			print Out "$hash{$nums}{$type}{$sample}\t";			
		}
		print Out "\n";
	}
}
close In;
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
	perl $Script 
Usage:
  Options:
	"i:s"=>\$fIn,	input marker 
	"sample:s"=>\$id,	# ,分开
	"o:s"=>\$fOut,	output file name

USAGE
        print $usage;
        exit;
}

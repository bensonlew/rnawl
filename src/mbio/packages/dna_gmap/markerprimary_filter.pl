#!/usr/bin/perl -w
use 5.10.0;
use strict;
use warnings;
use Getopt::Long;
my ($BEGIN_TIME,$fIn,$fOut,$ParentDepth,$OffDepth,$markertype,$indi);
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
	"Pdep=i"=>\$ParentDepth,
	"Odep=i"=>\$OffDepth,
	"mtype=s"=>\$markertype,
	"indi:s"=>\$indi,
			); 		#;/or &USAGE;
&USAGE unless ( $fIn and $fOut and $ParentDepth and $OffDepth and $markertype);
#check mtype;
my @ParentDepth = split /_/,$ParentDepth;
$ParentDepth[0]||=0;
my @OffDepth = split /_/,$OffDepth;
$OffDepth[0]||=0;
$markertype = uc($markertype);
my @tt=("SNP","INDEL","ALL");
my $val = 0;
foreach my $tru (@tt) {
	$val = 1 if ($tru eq $markertype);
}
die "\-markertype 必须输入SNP或者INDEL或者ALL\n" if ($val != 1);
####################################################存所有样本
my %indi;
if (defined $indi) {
	my @indi = split /,/,$indi;
	foreach my $i (@indi) {
		$indi{$i} = 1;
	}
}

open In,$fIn;
open Geno,">$fOut.primary.marker.xls";
# #MarkerID       marker_type     ParentDepth     type    16S456_1 
# #MarkerID       type    16S456_1        16S456_10
my @filter_num;
my @new_sample;
while (<In>) {
	chomp;
	next if ($_ eq "" || /^$/ );
	my ($id,$marker_type,$PDepth,$type,@info)=split(/\s+/,$_);
	if (/^#/) {
		while (my($idx, $val) = each @info){
			if (exists $indi{$val}) {
				push @new_sample,$info[$idx];
				push @filter_num, $idx;
			}
		}
		print Geno join("\t",$id,$type,@new_sample),"\n";
		next;
	}
		my @new_info;
		next if ($marker_type eq "SNP" && $markertype eq "INDEL");
		next if ($marker_type eq "INDEL" && $markertype eq "SNP");
		my @PDepth =split /_/,$PDepth;
		next if $PDepth[0] < $ParentDepth[0];		# 相等depth的也保留下来！
		next if $PDepth[1] < $ParentDepth[0];
		next if (defined $ParentDepth[1] && $PDepth[0] > $ParentDepth[1]);	# 比输入的最大值还大的过滤掉，相等的保留下来。
		next if (defined $ParentDepth[1] && $PDepth[1] > $ParentDepth[1]);
		foreach my $i (@filter_num) {
			if ($info[$i] eq "--"){
				push @new_info,"--";
			}else{
				my @dep_typ = split /\_/,$info[$i];		# 4_ll
				if ($dep_typ[0] < $OffDepth[0]) {
					push @new_info,"--";
				}elsif (defined $OffDepth[1] && $dep_typ[0] > $OffDepth[1]) {
					# $info[$i]="--";
					push @new_info,"--";
				}else{
					push @new_info,$dep_typ[1];
				}
			}
		}
	print Geno join("\t",$id,$type,@new_info),"\n";
}
close In;
#######################################################################################
# print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub Filter {
	my $depth_type = $_[0];
	my ($dep, $typ)=split /_/,$depth_type;
	my $tt="";
	if ($dep < $_[1]) {
		$tt = "--";
	}
	return $tt;
}

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
	perl $Script -i ./pop.primary.marker -o ./pop -Pdep 10 -Odep 1 -mtype ALL
Usage:
  Options:
	"i:s"=>\$fIn,	# pop.primary.marker
	"o:s"=>\$fOut,	# pop.primary.marker.xls
	"Pdep:s"=>\$ParentDepth,	# eg: 10_20 || _20 || 1 最好输入10_
	"Odep:s"=>\$OffDepth,		# eg: 1_20 || _20 ||1	最好输入1_
	"mtype:s"=>\$markertype,	# mtype = SNP||INDEL||ALL

USAGE
        print $usage;
        exit;
}

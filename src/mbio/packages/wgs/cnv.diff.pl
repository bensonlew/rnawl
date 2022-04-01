#!/usr/bin/perl -w
use 5.10.1;
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fIn1,$fIn2,$bool,$dOut,$Type,$Length);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"cnv1:s"=>\$fIn1,
	"cnv2:s"=>\$fIn2,
	"b|bool:s"=>\$bool,
	"t|type:s"=>\$Type,
	"l|leng:s"=>\$Length,
	"o:s"=>\$dOut,
			); 		#;/or &USAGE;
&USAGE unless ($fIn1 and $fIn2 and $dOut and $bool );
my ($Min,$Max)=split/-/,$Length if(defined $Length);
my %Type;
if(defined $Type){
	my @Type=split/,/,$Type;
	$Type{$_}=1 foreach (@Type);
}
$Min||=0;
$fIn1=DIR($fIn1);
$fIn2=DIR($fIn2);
$dOut=DIR($dOut);
my $name1=(split/[\/\.]/,$fIn1)[-4];
my $name2=(split/[\/\.]/,$fIn2)[-4];
open In,$fIn1;
my %cnv;
while (<In>) {
	chomp;
	next if ($_ eq "" ||/^$/ || /^#/);
	my ($chr1,$pos1,$pos2,$length,$type,$pvalue,$genenum,$genedetail)=split /\t/;
	$genedetail||="-";
	next if(defined $Max && $length >$Max);
	next if(defined $Type && !exists $Type{$type});
	$cnv{join("\t",$chr1,$pos1,$pos2,$length,$type)}=join("\t",$name1,$pvalue,$genenum,$genedetail) if ($length >=$Min);
}
close In;
open In,$fIn2;
open Out,">$dOut/$name1\_$name2.xls";
print Out join("\t","Chr1","Pos1","Pos2","Length","Type",$name1,$name2,"$name1\_pvalue","$name2\_pvalue","Gene number","Gene detail"),"\n";
while (<In>) {
	chomp;
	next if ($_ eq "" ||/^$/ || /^#/);
	my ($chr1,$pos1,$pos2,$length,$type,$pvalue,$genenum,$genedetail)=split /\t/;
	$genedetail||="-";
	next if(defined $Max && $length >$Max);
	next if(defined $Type && !exists $Type{$type});
	if($length >=$Min){
		if($bool =~ /TRUE/i && exists $cnv{join("\t",$chr1,$pos1,$pos2,$length,$type)}){
			my $sample1_pvalue=(split/\t/,$cnv{join("\t",$chr1,$pos1,$pos2,$length,$type)})[1];
			print Out join("\t",$chr1,$pos1,$pos2,$length,$type,"T","T",$sample1_pvalue,$pvalue,$genenum,$genedetail),"\n";
		}else{
			$cnv{join("\t",$chr1,$pos1,$pos2,$length,$type)}=join("\t",$name2,$pvalue,$genenum,$genedetail);
		}
	}
}
close In;
if($bool !~ /TRUE/i){
	foreach my $i (sort {$a cmp $b} keys %cnv){
		my ($sample,$pvalue,$genenum,$genedetail)=split/\t/,$cnv{$i};
		print Out "$i\t";
		if($sample eq $name1){
			print Out join("\t","T","F",$pvalue,"-"),"\t";
		}else{
			print Out join("\t","F","T","-",$pvalue),"\t";
		}
		print Out join("\t",$genenum,$genedetail),"\n";
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
	perl $Script -cnv1 Lands.cnv.anno -cnv2 A8-10.cnv.anno -bool TRUE/FALSE -o dir
Usage:
  Options:
  -cnv1	<file>	input cnv.anno
  -cnv2	<file>	input cnv.anno
  -b|bool <file>	default TEUE; or input FALSE
  -t|type	<str>	deletion|duplication
  -l|length	<str>	0: or :10000 or 3000:10000
  -o	<dir>	output dir
  -h         Help

USAGE
        print $usage;
        exit;
}

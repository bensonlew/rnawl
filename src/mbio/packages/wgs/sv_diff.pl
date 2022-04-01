#!/usr/bin/perl -w
use 5.10.1;
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fIn1,$fIn2,$bool,$fOut,$Type,$Length,$Depth1,$Depth2);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"sv1:s"=>\$fIn1,
	"sv2:s"=>\$fIn2,
	"d1|depth1:s"=>\$Depth1,
	"d2|depth2:s"=>\$Depth2,
	"b|bool:s"=>\$bool,
	"t|type:s"=>\$Type,
	"l|leng:s"=>\$Length,
	"o:s"=>\$fOut,
			); 		#;/or &USAGE;
&USAGE unless ($fIn1 and $fIn2 and $fOut and $bool);
my ($Min,$Max)=split/:/,$Length if(defined $Length);
my ($Dep1_min,$Dep1_max)=split/:/,$Depth1 if(defined $Depth1);
my ($Dep2_min,$Dep2_max)=split/:/,$Depth2 if(defined $Depth2);
$Min||=0;$Dep1_min||=0;$Dep2_min||=0;
my %Type;
if(defined $Type){
	my @Type=split/,/,$Type;
	$Type{$_}=1 foreach (@Type);
}
$fIn1=DIR($fIn1);
$fIn2=DIR($fIn2);
`touch $fOut` if(!-f $fOut );
$fOut=DIR($fOut);
my $name1=(split/[\/\.]/,$fIn1)[-4];
my $name2=(split/[\/\.]/,$fIn2)[-4];
open In,$fIn1;
my %sv;
while (<In>) {
	chomp;
	next if ($_ eq "" ||/^$/ || /^#/);
	my ($chr1,$pos1,$chr2,$pos2,$length,$type,$pvalue,$depth,$genenum,$genedetail)=split /\t/;
        $genedetail||="-";
	next if(defined $Max && abs($length)>$Max);
	next if(defined $Dep1_max && $depth>$Dep1_max);
	next if(defined $Type && !exists $Type{$type});
	$sv{join("\t",$chr1,$pos1,$chr2,$pos2,$length,$type)}=join("\t",$name1,$pvalue,$depth,$genenum,$genedetail) if (abs($length) >=$Min && $depth>=$Dep1_min);
}
close In;
open In,$fIn2;
open Out,">$fOut";
print Out join("\t","#Chr1","Pos1","Chr2","Pos2","Length","Type",$name1,$name2,"$name1\_pvalue","$name2\_pvalue","$name1\_depth","$name2\_depth","Gene number","Gene detail"),"\n";
while (<In>) {
	chomp;
	next if ($_ eq "" ||/^$/ || /^#/);
	my ($chr1,$pos1,$chr2,$pos2,$length,$type,$pvalue,$depth,$genenum,$genedetail)=split /\t/;
	$genedetail||="-";
	next if(defined $Max && abs($length) >$Max);
	next if(defined $Dep2_max && $depth>$Dep2_max);
	next if(defined $Type && !exists $Type{$type});
	if(abs($length) >=$Min and $depth >=$Dep2_min ){
		if($bool =~ /TRUE/i && exists $sv{join("\t",$chr1,$pos1,$chr2,$pos2,$length,$type)}){
			my (undef,$sample1_pvalue,$sample1_depth,undef)=split/\t/,$sv{join("\t",$chr1,$pos1,$chr2,$pos2,$length,$type)};
			print Out join("\t",$chr1,$pos1,$chr2,$pos2,$length,$type,"T","T",$sample1_pvalue,$pvalue,$sample1_depth,$depth,$genenum,$genedetail),"\n";
		}else{
			$sv{join("\t",$chr1,$pos1,$chr2,$pos2,$length,$type)}=join("\t",$name2,$pvalue,$depth,$genenum,$genedetail);
		}
	}
}
close In;
if($bool !~ /TRUE/i){
	foreach my $i (sort {$a cmp $b} keys %sv){
		my ($sample,$pvalue,$depth,$genenum,$genedetail)=split/\t/,$sv{$i};
		print Out "$i\t";
		if($sample eq $name1){
			print Out join("\t","T","F",$pvalue,"-",$depth,"-"),"\t";
		}else{
			print Out join("\t","F","T","-",$pvalue,"-",$depth),"\t";
		}
		print Out join("\t",$genenum,$genedetail),"\n";
	}
}
close Out;
#######################################################################################
# print STfOut "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
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
	perl $Script -sv1 Lands.sv.anno -sv2 A8-10.sv.anno -bool TRUE/FALSE -o dir
Usage:
  Options:
  -sv1	<file>	input sv.anno
  -sv2	<file>	input sv.anno
  -b	<file>	default TRUE; or input FALSE
  -o	<dir>	output dir
may no need:
  -d1	<str>	same as -l
  -d2	<str>	same as -l
  -t	<str>	DEL,INV,ITX,CTX,INS || DEL
  -l	<str>	0: or :10000 or 3000:10000
  -h         Help

USAGE
        print $usage;
        exit;
}

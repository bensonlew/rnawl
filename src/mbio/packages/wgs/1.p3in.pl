#!/usr/bin/env perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($diff,$dOut,$Range,$Tm1,$Tm2,$pair,$ref);
use Data::Dumper;
use List::Util qw(max);
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
		"help|?" =>\&USAGE,
		"d|diff:s"=>\$diff,
		"r|range:s"=>\$Range,
		"o|out:s"=>\$dOut,
		"T1|Tm1:s"=>\$Tm1,
		"T2|Tm2:s"=>\$Tm2,
		"p|pair:s"=>\$pair,
		"ref:s"=>\$ref,
                        ) or &USAGE;
&USAGE unless ($diff and $dOut and $Tm1 and $Tm2 and $pair and $Range and $ref);
my @range=split/[,-]/,$Range;
my $max=int(max(@range)/2);
$Range=~s/,/ /g;
############# `touch $dOut`if(!-f $dOut) || die "output file exists";
mkdir $dOut if(!-d $dOut);
$diff=DIR($diff);
$dOut=DIR($dOut);
open In,$diff;
open Out,">$dOut/variation.p3in";
##### file sort!!!!!!
my $Chr="";
my %Chrinfo;
my %CHR;
my $format="table";
my $nline=0;
while (<In>){
	chomp;
	next if ($_ eq ""|| /^$/);  # !!!! inputfile's header must include '^#';
	if (/^##/ && $nline == 0){
		$format="vcf";
		$nline=1;
		next;
	}
	next if (/^#/);
	$nline++;
	my ($chr,$pos,$ref,$alt);
	if ($format eq "vcf"){
		($chr,$pos,undef,$ref,$alt,undef)=split(/\s+/,$_);
	}else{
		($chr,$pos,$ref,$alt,undef)=split(/\s+/,$_);
	}
	$Chrinfo{$chr."_".$pos}=join("\t",$ref,$alt);
	push @{$CHR{$chr}},$pos;
}	
close In;
open Ref,"$ref";
$/=">";
while (<Ref>){
	chomp;
	next if ($_ eq "" || /^$/ );
	my($chr,@seq)=split/\n/,$_;
	my $seq=join("",@seq);
#	my($ref,$alt);  # add ID;
	if(exists $CHR{$chr}){
		foreach my $pos (sort @{$CHR{$chr}}){
			my($ref,$alt)=split /\t/,$Chrinfo{$chr."_".$pos};
			my $START_pos=$pos-$max-300;
			my $len=length($ref)-1+$max*2+300*2;
# eg:100-300 product;SNP;so $len=150*2+600=900;900 length template.
# left right +_300bp design primer.
# $max is  longest length's half  of primer product.
# end postion:$pos+length($ref)-1+$max+300 ($max=largest range num's half)
##$pos+length($ref)-1+20-($START_pos)
			my $TARGET=join(",",$pos-10-($START_pos),20+length($ref)-1);  # left center right.diatance equal.
			if($START_pos<0){
				$START_pos=0;
				$TARGET=join(",",$pos-2,18+length($ref)-1);
			}
			my $ID=join("_",$chr,$pos,$ref,$alt,$START_pos);
			my $template=substr($seq,$START_pos,$len); 
			print Out "PRIMER_SEQUENCE_ID=$ID
SEQUENCE_TEMPLATE=$template
SEQUENCE_TARGET=$TARGET
PRIMER_PRODUCT_SIZE_RANGE=$Range
PRIMER_MIN_TM=$Tm1
PRIMER_MAX_TM=$Tm2
PRIMER_NUM_RETURN=$pair
PRIMER_TASK=generic
PRIMER_MAX_END_STABILITY=500
=
";
		}
	}
}
close Out;
undef $/;
#######################################################################################
print "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
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

sub USAGE {#
        my $usage=<<"USAGE";
Description: 
Version:  $Script
Contact: qingmei.cui

Usage:
  Options:
-d	"d|diff"	<file>	diff.xls
-r	"r|range"	<file>	100-300,600-800 || 600-800 ; you can use "," separate 
-o 	"o|out "	<dir>	output dir/p3in
-T1	"T1|Tm1"	<float>	Tm1 Eg:57 || 57.5	
-T2	"T2|Tm2"	<float>	Tm2 Eg:63 || 63.5
-p	"p|pair"	<int>	primer pair numbers;Eg:3 ;p>=1 && p<=5
-ref	"ref"		<file>	input ref.fa

USAGE
        print $usage;
        exit;
}

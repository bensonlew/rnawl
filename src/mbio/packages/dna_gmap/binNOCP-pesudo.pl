#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fgeno,$dir,$Key,$fpos,$winsize,$stepsize);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$fgeno,
	"o:s"=>\$dir,
	"k:s"=>\$Key,
			) or &USAGE;
&USAGE unless ($fgeno and $dir and $Key);
mkdir $dir if (!-d $dir);
open In,$fgeno;
open Male,">$dir/$Key.male.matrix";
open Female,">$dir/$Key.female.matrix";
my @Head;
while (<In>) {
	#����ÿ��marker����Դ�жϣ�0�����޷�ȷ����1����ڶ�����-1�����1��Ⱦɫ��
	#ͬʱ������⽻�����δ�ж������࣬maybe not right
	chomp;
	next if ($_ eq "" || /^$/  );
	if (/^#/) {
		(undef,undef,@Head)=split(/\s+/,$_);
		print Male $_,"\n";
		print Female $_,"\n";
	}else{
		my ($id,$type,@geno)=split(/\s+/,$_);
		my @male;
		my @female;
		if ($type eq "aaxbb") {
			for (my $i=0;$i<@geno;$i++) {
				if ($geno[$i] eq "aa") {#nn��ll
					push @male,"aa 1";
					push @female,"aa 1";
				}elsif ($geno[$i] eq "bb") {#np��ll
					push @male,"bb -1";
					push @female,"bb -1";
				}elsif ($geno[$i] eq "ab") {#lm��nn
					push @male,"ab 1";
					push @female,"ab -1";
				}else{
					push @male,"-- 0";
					push @female,"-- 0";
				}
			}
			print Male "$id\t$type\t",join("\t",@male),"\n";
			print Female "$id\t$type\t",join("\t",@female),"\n";
		}
	}
}
close In;
close Male;
close Female;



#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub USAGE {#
        my $usage=<<"USAGE";
Contact:        long.huang\@majorbio.com;
Script:			$Script
Description: pesudo cross
	eg:
	perl $Script -i -o -k

Usage:
  Options:
  Options:
  -i	<file>	input genotype file name
			#MakrerID\tTYPE\tSample
  -o	<dir>	output dir
  -k	<str>	output file name 
  -h         Help
  
  -h         Help

USAGE
        print $usage;
        exit;
}
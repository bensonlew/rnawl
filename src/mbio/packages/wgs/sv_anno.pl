#!/usr/bin/env perl
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fSv,$fOut,$fStat,$fGff);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$fSv,
	"g:s"=>\$fGff,
	"o:s"=>\$fOut,
			) or &USAGE;
&USAGE unless ($fSv and $fOut  and $fGff );
open In,$fGff;
my %gene;
while (<In>) {
	chomp;
	next if ($_ eq "" ||/^$/ ||/^#/ );
	my ($chr,$soft,$type,$pos1,$pos2,undef,undef,undef,$info)=split(/\t/,$_);
	next if ($type ne "gene" && $type ne "mRNA");
	my $id;
	if ( $info =~ /Name=(\S+)\;/ || $info =~ /ID=(\S+)\;/ ) {
		$id=$1;	
	}
	$gene{$chr}{$pos1}{$pos2}=$id;
}
close In;
open In,$fSv;
open Out,">$fOut";
print Out "#chr\tpos1\tchr\tpos2\tlength\ttype\tpvalue\tdepth\tgene number\tgene\n";
my %stat;
my %length;
while (<In>) {
	chomp;
	next if ($_ eq "" ||/^$/ ||/^#/);
	my ($chr1,$pos1,undef,$chr2,$pos2,undef,$type,$length,$pvalue,$depth,undef)=split(/\t/,$_);
	my $info=join("\t",$chr1,$pos1,$chr2,$pos2,$length,$type,$pvalue);
	if ($chr1 ne $chr2) {
		my $postmp=$pos2;
		$pos2=$pos1+$length;
		my $chr=$chr1;
		my @gene;
		foreach my $pos3 (sort {$a<=>$b} keys %{$gene{$chr}}) {
			foreach my $pos4 (sort {$a<=>$b} keys %{$gene{$chr}{$pos3}}) {
				if (($pos3 >= $pos1 && $pos3 <= $pos2) || ($pos4 >= $pos1 && $pos4<=$pos2)) {
					push @gene,$gene{$chr}{$pos3}{$pos4};
				}
				last if ($pos3 > $pos2);
			}	
		}
		$chr=$chr2;
		$pos1=$postmp;
		$pos2=$postmp+$length;;
		foreach my $pos3 (sort {$a<=>$b} keys %{$gene{$chr}}) {
			foreach my $pos4 (sort {$a<=>$b} keys %{$gene{$chr}{$pos3}}) {
				if (($pos3 >= $pos1 && $pos3 <= $pos2) || ($pos4 >= $pos1 && $pos4<=$pos2)) {
					push @gene,$gene{$chr}{$pos3}{$pos4};
				}
				last if ($pos3 > $pos2);
			}	
		}
		 if (scalar @gene >0){
			print Out join("\t",$info,scalar @gene,join(":",@gene)),"\n";
		 }else{
			print Out join("\t",$info,scalar @gene,"-"),"\n";
		}
	}else{
		my $chr=$chr1;
		my @gene;
		foreach my $pos3 (sort {$a<=>$b} keys %{$gene{$chr}}) {
			foreach my $pos4 (sort {$a<=>$b} keys %{$gene{$chr}{$pos3}}) {
				if (($pos3 >= $pos1 && $pos3 <= $pos2) || ($pos4 >= $pos1 && $pos4<=$pos2)) {
					push @gene,$gene{$chr}{$pos3}{$pos4};
				}
				last if ($pos3 > $pos2);
			}	
		}
		 if (scalar @gene >0){
			print Out join("\t",$info,$depth,scalar @gene,join(":",@gene)),"\n";
		 }else{
			print Out join("\t",$info,$depth,scalar @gene,"-"),"\n";
		}
	}
}
close Out;
close In;
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub USAGE {#
        my $usage=<<"USAGE";
Contact:        long.huang\@majorbio.com;
Script:			$Script
Description:
	fq thanslate to fa format
	eg:
	perl $Script -i -o -k -c

Usage:
  Options:
  -i	<file>	input cnv name
  -g	<file>	input gff name
  -o	<file>	output merge  file
  -s	<file>	output stat file
  -h         Help

USAGE
        print $usage;
        exit;
}

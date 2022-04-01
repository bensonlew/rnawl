#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fIn1,$fIn2,$type,$fOut);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$fIn1,
	"k:s"=>\$fIn2,
	"type:s"=>\$type,
	"o:s"=>\$fOut
			) or &USAGE;
&USAGE unless ($fIn1 and $fIn2 and $fOut );
open In1,$fIn1;
open In2,$fIn2;
my @in1=split(/\//,$fIn1);
my @in2=split(/\//,$fIn2);
my ($name1,undef)=split(/\./,$in1[-1]);
my ($name2,undef)=split(/\./,$in2[-1]);
my (%hash1,%hash2,%hash3,%hash4);
while (<In1>){
	chomp;
	next if (/^$/ or "" or /#/);
	chomp;
	my @line=split(/\t/);
	my $chr_pos=join(":",$line[0],$line[4],$line[5]);
	my ($str,$num,undef)=split(/(\d+)/,$line[0]);
	$hash1{$chr_pos}=$_;
	$hash2{$str}{$num}{$line[4]}{$line[5]}=1;
	$hash4{$line[0]}=1;
}
close In1;
while(<In2>){
	chomp;
	next if (/^$/ or "" or /#/);
	chomp;
	my @line=split(/\t/);
	my $chr_pos=join(":",$line[0],$line[4],$line[5]);
	my ($str,$num,undef)=split(/(\d+)/,$line[0]);
	$hash2{$str}{$num}{$line[4]}{$line[5]}=1;
	$hash3{$chr_pos}=$_;
	$hash4{$line[0]}=1;
}
close In2;
$_=$type;
$type=~s/\s+//g;
if ( $type eq "same"){
	open Write,">$fOut/ssr.$type.result.xls";
	print Write"#Chr\tSSR.nr\tSSR type\tSSR\tSize\tStart\tEnd\t$name1\t$name2\tFORWARD PRIMER1 (5'-3')\tTm\tGC(%)\tLength(bp)\tREVERSE PRIMER1 (5'-3')\tTm\tGC(%)\tLength(bp)\tPRODUCT SIZE\n";
	my $ssr_nr=0;
	foreach my $key1 (sort {$a cmp $b} keys %hash2){
		foreach my $key2 (sort {$a <=> $b} keys %{$hash2{$key1}}){
			foreach my $key3 (keys %{$hash2{$key1}{$key2}}){
				$ssr_nr++;
				foreach my $key4 (sort {$a <=> $b} keys %{$hash2{$key1}{$key2}{$key3}}){
					if (exists $hash1{"$key1$key2:$key3:$key4"} and exists $hash3{"$key1$key2:$key3:$key4"}){
						my @one=split(/\t/,$hash1{"$key1$key2:$key3:$key4"},8);
						my @two=split(/\t/,$hash3{"$key1$key2:$key3:$key4"},8);
						if ($key1 eq "sca" or $key1 eq "chr"){
							print Write"ref_$one[0]\t$ssr_nr\t$one[2]\t$one[3]\t$one[4]\t$one[5]\t$one[6]\tT\tT\t$one[-1]\n";
						}
						else{
							print Write"$name1\_$one[0]\t$ssr_nr\t$one[2]\t$one[3]\t$one[4]\t$one[5]\t$one[6]\tT\tT\t$one[-1]\n";
						}
					}
				}
			}
		}
	}
	close Write;
}
else{
	open Write,">$fOut/ssr.$type.result.xls";
	print Write"#Chr\tSSR.nr\tSSR type\tSSR\tSize\tStart\tEnd\t$name1\t$name2\t$name1 FORWARD PRIMER1 (5'-3')\tTm\tGC(%)\tLength(bp)\t$name1 REVERSE PRIMER1 (5'-3')\tTm\tGC(%)\tLength(bp)\tPRODUCT SIZE\t$name2 FORWARD PRIMER1 (5'-3')\tTm\tGC(%)\tLength(bp)\t$name2 REVERSE PRIMER1 (5'-3')\tTm\tGC(%)\tLength(bp)\tPRODUCT SIZE\n";
	my $ssr_nr=0;
	foreach my $key1 (sort {$a cmp $b} keys %hash2){
		foreach my $key2 (sort {$a <=> $b} keys %{$hash2{$key1}}){
			foreach my $key3 (keys %{$hash2{$key1}{$key2}}){
				$ssr_nr++;
				foreach my $key4 (sort {$a <=> $b} keys %{$hash2{$key1}{$key2}{$key3}}){
					if(exists $hash1{"$key1$key2:$key3:$key4"} and !exists $hash3{"$key1$key2:$key3:$key4"}){
						my @one=split(/\s+/,$hash1{"$key1$key2:$key3:$key4"},8);
						if ($key1 eq "sca" or $key1 eq "chr"){
							print Write"ref_$one[0]\t$ssr_nr\t$one[2]\t$one[3]\t$one[4]\t$one[5]\t$one[6]\tT\tF\t$one[-1]\t--\t--\t--\t--\t--\t--\t--\t--\t--\n";
						}
						else{
							print Write"$name1\_$one[0]\t$ssr_nr\t$one[2]\t$one[3]\t$one[4]\t$one[5]\t$one[6]\tT\tF\t$one[-1]\t--\t--\t--\t--\t--\t--\t--\t--\t--\n";
						}
						next;
					}
					if(!exists $hash1{"$key1$key2:$key3:$key4"} and exists $hash3{"$key1$key2:$key3:$key4"}){
						my @two=split(/\s+/,$hash3{"$key1$key2:$key3:$key4"},8);
						if ($key1 eq "sca" or $key1 eq "chr"){
							print Write"ref_$two[0]\t$ssr_nr\t$two[2]\t$two[3]\t$two[4]\t$two[5]\t$two[6]\tF\tT\t--\t--\t--\t--\t--\t--\t--\t--\t--\t$two[-1]\n";
						}
						else{
							print Write"$name1\_$two[0]\t$ssr_nr\t$two[2]\t$two[3]\t$two[4]\t$two[5]\t$two[6]\tF\tT\t--\t--\t--\t--\t--\t--\t--\t--\t--\t$two[-1]\n";
						}
						next;
					}
				}
			}
		}
	}
	close Write;
}

#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub USAGE {#
        my $usage=<<"USAGE";
Contact:        chongqing.shi\@majorbio.com;
Script:			$Script
Description:
	fq thanslate to fa format
	eg:
	perl $Script -i -o -k -c

Usage:
  Options:
  -i	<file>	misa.out1 file name
  -k	<file>	misa.out1 file name
  -type	<flag>	same or diff
  -o	<file>	out file name
  -h         Help

USAGE
        print $usage;
        exit;
}

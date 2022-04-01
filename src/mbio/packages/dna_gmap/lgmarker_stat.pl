#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($lg,$info,$out);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"lg:s"=>\$lg,
	"out:s"=>\$out,
	"info:s"=>\$info,
			) or &USAGE;
&USAGE unless ($info and $lg and $out);

my %stat;
open IN,$info;
while (<IN>){
	chomp;
	next if ($_ eq ""|| /^$/|| /^#/);
	my ($marker,$gtype,$vtype,undef)=split(/\s+/,$_,4);#MarkerID       Gtype   Vtype   Pdep    Mdep    16S456_1        16S456_2
	$stat{$marker}=$vtype;
}
close IN;

my %for;
my $formation;
my $bin||=0;;
open OUT,">$out";
print OUT "LG ID\tMarker Number\tSNP Number\tInDel Number\n";
open In,$lg;
$/=">";
while(<In>){
	chomp;
	next if ($_ eq ""|| /^#/ || /^$/);
	my($lgid,$lgnum,@marker)=split(/\s+/,$_);
	my $type=$marker[-1];
	$bin=1 if($type=~/(\w+)(-)(\d+)/);
	$formation=join(",",@marker);
	$for{$lgid}=$formation;
}
close In;
# print Dumper \%for;
#print $bin,"\n",$formation,"\n";
foreach my$lgid(sort{$a cmp $b}keys%for){
	my @indel=();
	my @snp=();
	#my ($marknum,$markers)=split(/\t/,$for{$lgid});
	my @mark=split(/\,/,$for{$lgid});
	if($bin eq "0"){
		foreach my $marker(keys%stat){
			foreach(@mark){
				if($_ eq $marker){
					push @indel,$marker if($stat{$marker} eq "INDEL");
					push @snp,$marker if($stat{$marker} eq "SNP");
				}
			}
		}
		#print OUT "$lgid\t",scalar@mark,"\t",scalar@snp,"\t",scalar@indel,"\n";
	}else{
		foreach my $marker(keys%stat){
			my ($chrid,$pos)=split(/\_/,$marker);
			foreach(@mark){
                next if($_ !~ /\_/);
				my($chr,$win)=split(/\_/,$_);
				my($start,$end)=split(/\-/,$win);
				if($chr eq $chrid){
					if($pos > $start && $pos <=$end){
						push @indel,$marker if($stat{$marker} eq "INDEL");
						push @snp,$marker if($stat{$marker} eq "SNP");
					}
				}
			}
		}
	}
	my $marknum=scalar@indel + scalar@snp;
	print OUT "$lgid\t",$marknum,"\t",scalar@snp,"\t",scalar@indel,"\n";
}
close OUT;
#######################################################################################
#print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub USAGE {#
        my $usage=<<"USAGE";
Contact:        caixia.tian\@majorbio.com;
Script:			$Script
Description:
	fq thanslate to fa format
	eg:
	perl $Script -i -o -k -c

Usage:
  Options:
  -lg	<file>	Total.lg file
  -info	<file>  detail.info file
  -out	<file>	output result file
  -h         Help

USAGE
        print $usage;
        exit;
}

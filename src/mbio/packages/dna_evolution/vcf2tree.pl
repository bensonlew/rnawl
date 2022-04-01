#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fIn,$fOut);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$fIn,
	"o:s"=>\$fOut
			) or &USAGE;
&USAGE unless ($fIn and $fOut);
open In,$fIn;
my %seq;
my @Indi;
my %BASE=(
	"AA"=>"A","GG"=>"G","CC"=>"C","TT"=>"T",
	"AT"=>"W","AG"=>"R","AC"=>"M",
	"CG"=>"S","CT"=>"Y",
	"GT"=>"K"
);
while (<In>) {
	chomp;
	next if ($_ eq "" || /^$/ ||/^##/);
	if (/^#/) {
		(undef,undef,undef,undef,undef,undef,undef,undef,undef,@Indi)=split(/\t/,$_);
		foreach my $id (@Indi) {
			$seq{$id}="";
		}
	}else{
		my ($CHROM,$POS,$ID,$REF,$ALT,$QUAL,$FILTER,$INFO,$FORMAT,@indi)=split(/\t/,$_);
		my @ale=split(/\,/,join(",",$REF,$ALT));
		for (my $i=0;$i<@indi;$i++) {
			my $geno=(split(/\:/,$indi[$i]))[0];
			my ($g1,$g2)=split(/\//,$geno);
			if ($g1 eq ".") {
				$seq{$Indi[$i]}.="?";
			}else{
				if (!exists $BASE{join("",sort($ale[$g1],$ale[$g2]))}){
					$seq{$Indi[$i]}.="?";
					next;
				};
				if (!exists $seq{$Indi[$i]}) {
					print $Indi[$i],"\n";die;
				}
				$seq{$Indi[$i]}.=$BASE{join("",sort($ale[$g1],$ale[$g2]))};
			}
		}
	}
}
close In;
open Out,">$fOut.phylip";
open FA,">$fOut.fasta";
my $head=0;
foreach my $id (sort keys %seq) {
	if ($head == 0) {
		print Out scalar keys %seq," ",length($seq{$id}),"\n";
		$head =1;
	}
	print Out $id," ",$seq{$id},"\n";
	print FA ">$id\n$seq{$id}\n";
}
close Out;
close FA;
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
  -i	<file>	input file name
  -o	<file>	split windows sh
  -h         Help

USAGE
        print $usage;
        exit;
}

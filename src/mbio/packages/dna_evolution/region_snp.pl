#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fin,$fout,$min);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	#"i:s"=>\$fIn,
	#"i2:s"=>\$fIn,
	"input:s"=>\$fin,
	"output:s"=>\$fout,
	"select:s"=>\$min,
			) or &USAGE;
&USAGE unless ($fout);
open IN,$min;
my %region;
while (<IN>) {
	chomp;
	next if ($_ eq "" || /^$/ || /pos/);
	s/"//g;
	my ($chr,$pos1,$pos2);
	if ($min=~/pi_tajimaD_fst.select/) {
		($chr,$pos1,undef)=split/\s+/;
		$pos2 = $pos1 + 10000 - 1;
	}else{
		($chr,$pos1,$pos2,undef)=split(/\s+/,$_);
	}
	my $regioned=0;
	foreach my $region (sort keys %{$region{$chr}}) {
		my ($pos3,$pos4)=split(/\s+/,$region);
		#print $pos1,"\t",$pos2,"\t",$pos3,"\t",$pos4,"\n";
		if ($pos1 >= $pos3 && $pos1 <= $pos4) {
			my ($p1,$p2,$p3,$p4)=sort {$a<=>$b} ($pos1,$pos2,$pos3,$pos4);
			my $newregion=join("\t",$p1,$p4);
			delete $region{$chr}{$region};
			$region{$chr}{$newregion}++;
			$regioned=1;
		}
	}
	if ($regioned == 0) {
		$region{$chr}{join("\t",$pos1,$pos2)}++;
	}
}
close IN;
#print Dumper \%region;die;
open IN,$fin;
open OUT,">$fout";
while (<IN>) {
	chomp;
	my ($chr,$pos,$id,$info)=split/\s+/,$_,4;
	foreach my $region (sort keys %{$region{$chr}}) {
		my ($pos1,$pos2)=split(/\s+/,$region);
		if ($pos >= $pos1 && $pos <= $pos2) {
			print OUT "$id\n";
		}
	}
}
close IN;
close OUT;
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub USAGE {#
        my $usage=<<"USAGE";
Contact:        minghao.zhang\@majorbio.com;
Script:			$Script
Description:
	fq thanslate to fa format
	eg:
	perl $Script -i -o -k -c

Usage:
  Options:
	"input:s"=>\$fin,  pop.recode.vcf                                                            
    "output:s"=>\$fout,  output file
	"select:s"=>\$selectin,	*.select
  -h         Help

USAGE
        print $usage;
        exit;
}

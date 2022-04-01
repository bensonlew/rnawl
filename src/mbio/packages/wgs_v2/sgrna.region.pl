#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fin,$fout,$sgr);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$fin,
	"sgr:s"=>\$sgr,
	"o:s"=>\$fout,
			) or &USAGE;
&USAGE unless ($fout);
open In,$sgr;
#$/ = ">";
my $sgrna;
while (<In>) {
	chomp;
	$sgrna=$_;
	#next if ($_ eq "" || /^$/);
	#(undef,$sgrna) = split(/\n/,$_,2);
}
close In;

open IN,"gunzip -c $fin |";
my $n;
open OUT,">$fout";
while (<IN>) {
	chomp;
	$n++;
	my (undef,undef,undef,$chr,$pos1,$pos2,undef)=split/\s+/;
	print OUT "$chr\t$pos1\t$pos2\tR$n\t$sgrna\tNA\tNA\n";
}
close IN;
close OUT;
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub USAGE {#
        my $usage=<<"USAGE";
Contact:        meng.luo\@majorbio.com;
Script:			$Script
Description:
	fq thanslate to fa format
	eg:
	perl $Script -i -o -k -c

Usage:
  Options:
  -i risearch result
  -sgr sgrna sequence (no PAM sequence)
  -o output file name 
  -h         Help

USAGE
        print $usage;
        exit;
}

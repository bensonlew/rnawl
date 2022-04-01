#/usr/bin/perl -w
use strict;
use warnings;
my ($fa,$position,$connec,$symbol,$out);
use Getopt::Long;
use Data::Dumper;
use File::Basename qw(basename dirname);
GetOptions(
	"help|?" =>\&USAGE,
	"fa:s"=>\$fa,
	"position:s"=>\$position,
	"connect:s" =>\$connec,
	"symbol:s" =>\$symbol,
	"out:s"=>\$out,
			) or &USAGE;
&USAGE unless ($fa and $symbol);

$position||="before";
$connec ||="_";
$out ||="NewID.fa";

open FA, "<$fa"||die "no such files";
open OUT,">$out";
$/ = ">";
my $newID;
while (<FA>){
	chomp;
	next if ($_ eq ""||/^$/);
	my($seq,@seq) = split(/\n/,$_);
	my($seqid,@id) = split(/[\s\t]/,$seq);
	if ($position eq "before" ){
		$newID = $symbol.$connec.$seqid
	}else{
		$newID = $seqid.$connec.$symbol
	}
#	print $seq;
	print OUT ">",$newID," ",join(" ",@id),"\n",@seq,"\n";
}
$/ = "\n";
close FA;
close OUT;


sub USAGE {
        my $usage=<<"USAGE";
Contact:        licui.li\@majorbio.com;
Script:		perl faChangeHeader.pl -fa <FA> -symbol \_ 
Description: change FA sequence header
Usage:
  Required:
	-fa	<file>	input FA file
	-symbol	<str>	the added character
  Optional:
	-connect	<>	any connector,default "_"
	-position	<A or B>	the position of symbol added,A with After and B with Before
	-out	<file>	outputfile
	-h         Help

USAGE
        print $usage;
        exit;
}

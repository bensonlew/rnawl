#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fGenome,$fOut,$Gff,$match);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$fGenome,
	"o:s"=>\$fOut,
	"g:s"=>\$Gff,
	"f:s"=>\$match
	) or &USAGE;
&USAGE unless ($fGenome and $fOut and $Gff);
my %change;
if ($match) {
	open In,$match;
	while (<In>) {
		chomp;
		next if ($_ eq ""||/^$/);
		my ($id,$change)=split(/\s+/,$_);
		$change{$id}=$change;
	}
	close In;
}
open In,$fGenome;
open Out,">$fOut.fa";
$/=">";
my $n=0;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/ ||/^#/);
	my ($info,@seq)=split(/\n/,$_);
	my $id=(split(/\s+/,$info))[0];
	if (!exists $change{$id}) {
		$n++;
		$change{$id}="sca$n";
	}
	print Out ">$change{$id}\n";
	print Out join("\n",@seq),"\n";
}
close In;
close Out;
open In,$Gff;
open Out,">$fOut.gff";
$/="\n";
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/||/^#/) ;
	my ($id,@info)=split(/\s+/,$_);
	die "gff && ref file doesn't match\n" if (!exists $change{$id});
	print Out join("\t",$change{$id},@info),"\n";
}
close In;

close Out;

open Out,">$fOut.changelog";
foreach my $id (sort keys %change) {
	print Out join("\t",$id,$change{$id}),"\n";
}
close Out;
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub USAGE {#
        my $usage=<<"USAGE";
Contact:        long.huang\@majorbio.com;
Script:			$Script
Description:
	reformat genome,rename scaffold name at genome fa file and gff file
eg:
	perl $Script -i Genome.fa -g Genome.gff -k keyname -o dir 

Usage:
  Options:
  -i	<file>	input genome name,fasta format,
  -g	<file>	input genome gff file,
  -o	<str>	output file prefix
  -f	<file>	chromosome change file

  -h         Help

USAGE
        print $usage;
        exit;
}

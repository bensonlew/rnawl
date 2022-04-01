#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fIn,$fOut,$tophit,$topmatch,$eval,$fa);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"d:s"=>\$fIn,
	"o:s"=>\$fOut,
	"i:s"=>\$fa,
	"hit:s"=>\$tophit,
	"match:s"=>\$topmatch,
	"evalue:s"=>\$eval,
			) or &USAGE;
&USAGE unless ($fIn and $fOut);
open In,$fa;
my %gene;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/);
	my ($fafile,undef)=split(/\s+/,$_);
	$/=">";
	open Seq,$fafile;
	while (<Seq>) {
		chomp;
		next if ($_ eq ""||/^$/);
		my ($id,undef)=split(/\s+/,$_);
		$gene{$id}=1;
	}
	close Seq;
	$/="\n";
}
close In;
my %EggNOG=(
	"J"=>"Translation, ribosomal structure and biogenesis",
	"A"=>"RNA processing and modification",
	"K"=>"Transcription",
	"L"=>"Replication, recombination and repair",
	"B"=>"Chromatin structure and dynamics",
	"D"=>"Cell cycle control, cell division, chromosome partitioning",
	"Y"=>"Nuclear structure",
	"V"=>"Defense mechanisms",
	"T"=>"Signal transduction mechanisms",
	"M"=>"Cell wall/membrane/envelope biogenesis",
	"N"=>"Cell motility",
	"Z"=>"Cytoskeleton",
	"W"=>"Extracellular structures",
	"U"=>"Intracellular trafficking, secretion, and vesicular transport",
	"O"=>"Posttranslational modification, protein turnover, chaperones",
	"C"=>"Energy production and conversion",
	"G"=>"Carbohydrate transport and metabolism",
	"E"=>"Amino acid transport and metabolism",
	"F"=>"Nucleotide transport and metabolism",
	"H"=>"Coenzyme transport and metabolism",
	"I"=>"Lipid transport and metabolism",
	"P"=>"Inorganic ion transport and metabolism",
	"Q"=>"Secondary metabolites biosynthesis, transport and catabolism",
	"R"=>"General function prediction only",
	"S"=>"Function unknown",
);

my %anno;
my @blast=glob("$fIn/*.emapper.annotations");
foreach my $blast (@blast) {
	open In,$blast;
	while (<In>) {
		chomp;
		next if ($_ eq ""||/^$/||/^#/);
		my @info=split(/\t/,$_);
		if (scalar @info < 13) {
			$anno{$info[0]}="--\t--";
		}else{
			my @type=split(/\,/,$info[11]);
			my @detail;
			foreach my $type (@type) {
				$type=~s/\s//g;
				push @detail,$EggNOG{$type};
			}
			my $NOG=(split(/\,/,$info[9]))[0];
			my $NOGExpress=$info[12];
			$anno{$info[0]}=join("\t",join(",",@type).":".join(";",@detail),"$NOG:$NOGExpress");
		}
	}
	close In;
}

open Out,">$fOut";
print Out "#GeneID\tEGGNOG\tEGGNOG_ANNO\n";
foreach my $query (sort keys %gene) {
	$anno{$query}||="--\t--";
	print Out $query,"\t",$anno{$query},"\n";
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
	So this program is written, to get the information and list them in lines saved in a file from the blast m0 format outfile
	eg:
	perl $Script -i $fIn -o $fOut

Usage:
  Options:
	-i	<file>	input file name
	-o	<file>	output file name
	-tophit	<num>	input set how many subjects for a query to be displayed,default 1
	-topmatch	<num>	input set suits(results of one subject match one query) to be displayed.match number,defualt 1
	-evalue	<num>	input max evalue,default 10e-5
  -h         Help

USAGE
        print $usage;
        exit;
}

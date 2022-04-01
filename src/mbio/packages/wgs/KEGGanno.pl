#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($dIn,$fOut,$fa);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"d:s"=>\$dIn,
	"o:s"=>\$fOut,
	"i:s"=>\$fa,
			) or &USAGE;
&USAGE unless ($dIn and $fOut and $fa);
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
$/="\n";
my %anno;
my @kobas=glob("$dIn/*.kobas");
open Out,">$fOut";
print Out "#query\tKoID\tKoanno\n";
foreach my $kobas (@kobas) {
	$/="///";
	open In,$kobas;
	while (<In>) {
		chomp;
		next if ($_ eq "" || /^$/ || /^#/);
		my @line=split(/\n/,$_);
		next if (scalar @line <4);
		my (undef,$queid,$KID,$pathway)=split(/\n/,$_,4);
		$queid=(split(/\t/,$queid))[-1];
		$KID=(split(/\t/,$KID,2))[-1];
		my @Path=split(/\n/,$pathway);
		my @koid;
		my @anno;
		foreach my $path (@Path) {
				my @info=split(/\t/,$path);
				push @koid,$info[3];
				push @anno,$info[1];
		}
		$anno{$queid}=join("\t",join(",",@koid),join(":",@anno));
	}
	close In;
}
foreach my $id (sort keys %gene) {
	if (!exists $anno{$id}) {
		print Out $id,"\t","--\t--\n";
	}else{
		print Out $id,"\t",$anno{$id},"\n";
	}
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
	eg:
	perl $Script -i -o 

Usage:
  Options:
  -i	<dir>	input *.kobas file name
  -o	<file>	output file name
  -h         Help

USAGE
        print $usage;
        exit;
}

#!/usr/bin/env perl
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fCnv,$fOut,$fStat,$fGff);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$fCnv,
	"g:s"=>\$fGff,
	"o:s"=>\$fOut,
			) or &USAGE;
&USAGE unless ($fCnv and $fOut and $fGff );
open In,$fGff;
my %gene;
while (<In>) {
	chomp;
	next if ($_ eq "" ||/^$/ ||/^#/ );
	my ($chr,$soft,$type,$pos1,$pos2,undef,undef,undef,$info)=split(/\s+/,$_);
	next if ($type ne "gene" && $type ne "mRNA");
	my $id;
	if ( $info =~ /Name=([^;]*);/ || $info =~ /ID=([^;]*);/ ) {
		$id=$1;	
	}
	$gene{$chr}{$pos1}{$pos2}=$id;
}
close In;
open In,$fCnv;
open Out,">$fOut";
print Out "#chr\tpos1\tpos2\tlength\ttype\tpvalue\tgene number\tgene\n";
my %stat;
my %length;
while (<In>) {
	chomp;
	next if ($_ eq "" ||/^$/ ||/^$/);
	my ($type,$pos,$length,undef,$pvalue,undef)=split(/\t/,$_);
	my ($chr,$position)=split(/\:/,$pos);
	my ($pos1,$pos2)=split(/\-/,$position);
	my $info=join("\t",$chr,$pos1,$pos2,$length,$type,$pvalue);
	my @gene;
	foreach my $pos3 (sort {$a<=>$b} keys %{$gene{$chr}}) {
		foreach my $pos4 (sort {$a<=>$b} keys %{$gene{$chr}{$pos3}}) {
			if (($pos3 > $pos1 && $pos3 < $pos2) || ($pos4 > $pos1 && $pos4<$pos2)) {
				push @gene,$gene{$chr}{$pos3}{$pos4};
			}
			last if ($pos3 > $pos2);
		}	
	}
	print Out join("\t",$info,scalar @gene,join(":",@gene)),"\n";
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

#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($ref,$out,$gff);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$ref,
	"o:s"=>\$out,
	"g:s"=>\$gff,
	) or &USAGE;
&USAGE unless ($ref and $out and $gff);
open In,$ref;
$/=">";
my %seq;
while (<In>) {
	chomp;
	next if ($_ eq "" ||/^$/);
	my ($id,@seq)=split(/\n/,$_);
	my $seq=join("",@seq);
	$seq{$id}=$seq;
}
close In;
open In,$gff;
open Out,">$out";
$/="\n";
my $flag="gene";
my %out;
my @line;
my $n=0;
my $m=0;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/||/^#/);
	my (undef,undef,$types,undef,undef,undef,undef,undef,undef)=split(/\t/,$_);
	next if ($types eq "region");
	$n++;	
	if ($n == 1 && $types ne "gene") {
		$flag="mRNA";
	}
	if ($types eq $flag) {
		if (scalar @line == 0) {
			push @line,$_;
		}else{
			$m++;
			my $Entrez||="--";
			my $Geneba||="--";
			my $Chromo||="--";
			my $START||="--";
			my $END||="";
			my %Protei;
			my @Protei;
			my $GeneNa||="--";
			my %Transc;
			my @Transc;
			my @Geneba;
			my %Geneba;
			my $id;
			foreach my $line (@line) {
				my ($chr,$source,$type,$start,$end,undef,undef,undef,@info)=split(/\t/,$line);
				my $info=join(" ",@info);
				if ($type eq $flag) {
					if ($info =~ /Name=([^;,]*)/) {
						$GeneNa=$1;
					}
					$START=$start;
					$END=$end;
					$Chromo=$chr;
					if ($info=~/GeneID:([^;,]*)/) {
						$Entrez=$1;
					}
					if ($info=~/ID=([^;,]*)/) {
						$id=$1;
					}
				}
				if ($type ne "CDS" && $type ne "exon") {
					if ($info=~/Genbank:([^;,]*)/) {
						if (!exists $Geneba{$1}) {
							push @Geneba,$1;
						}
						$Geneba{$1}=1;
					}
				}
				if ($type eq "exon") {
					if ($info=~/;transcript_id=([^;,]*)/) {
						if (!exists $Transc{$1}) {
							push @Transc,$1;
						}
						$Transc{$1}=1;
					}
				}

				if ($type eq "CDS") {
					if ($info=~/;protein_id=([^;,]*)/) {
						if (!exists $Protei{$1}) {
							push @Protei,$1;
						}
						$Protei{$1}=1;
					}
				}

			}
			my $outid=join(":",$id,join("|",$GeneNa,$Entrez,join(";",@Geneba),join(";",@Transc),join(";",@Protei)),$Chromo,$START,$END);
			my @outid=split(/\s+/,$outid);
			$outid=join("\_",@outid);
			print Out ">$outid\n".substr($seq{$Chromo},$START,$END-$START+1),"\n";
			@line=();
			push @line,$_;
		}
	}else{
		push @line,$_;
	}
}
close In;
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

  -h         Help

USAGE
        print $usage;
        exit;
}

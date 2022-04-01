#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
my ($ref,$out,$gff);
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$ref,
	"o:s"=>\$out,
	"g:s"=>\$gff,
	) or &USAGE;
&USAGE unless ($ref and $out and $gff);
###########################################################
$ref = ABSOLUTE_DIR($ref);
$gff = ABSOLUTE_DIR($gff);

open IN,$ref;
my %seq;       #%seq, id-->sequence
$/='>';
while (<IN>) {
	chomp;
	next if ($_ eq "" ||/^$/||/^#/);
	my ($title,@seq) = split(/\n/,$_);
	my ($id,undef) = split /\s+/,$title;
	$id=~s/\:/\./g;         ###
	next if ($id=~/cds/);
	my $seq = join("",@seq);
	$seq{$id} = $seq;
}
close IN;

open IN,$gff;
open OUT,">$out";
$/ = "\n";
my %genelist;  #%genelist, id-->geneinfo
my %rnaparent; #%rnalist, rnaid-->rnaparent
my %location;  #%location, id-->location
my %Genebanks; #%Genebanks, parentid-->genebank1;genebank2
my %Transcs;   #%Transcs, parentid-->trans1;trans2
my %Proteins;  #%Proteins, parentid-->pro1;pro2
while (<IN>) {
	chomp;
	next if ($_ eq ""||/^$/||/^#/);
	my ($chr,$source,$type,$start,$end,undef,undef,undef,$info) = split(/\t/,$_);
	my $id;       #ID=
	my $Parent;   #Parent=
	my $GeneName; #Name=
	my $Entrez;   #GeneID:

	if ($info=~/ID=([^;,]*)/) {
		$id = $1;
		$id=~s/\:/\./g;
	}else{
		next;
	}

	if ($info=~/Parent=([^;,]*)/) {
		$Parent = $1;
		if ($info=~/Genbank:([^;,]*)/) {
			$Genebanks{$Parent}{$1}++;
		}
		if ($info=~/;transcript_id=([^;,]*)/) {
			$Transcs{$Parent}{$1}++;
		}
		if ($info=~/;protein_id=([^;,]*)/) {
			$Proteins{$Parent}{$1}++;
		}
	}
	if ($id=~/^gene/) {
		if ($info =~/Name=([^;,]*)/) {
			$GeneName = $1;
		}else{
			$GeneName = '--';
		}
		if ($info=~/GeneID:([^;,]*)/) {
			$Entrez = $1;
		}else{
			$Entrez = '--';
		}
		$location{$id} = "$chr:$start:$end";
		$genelist{$id} = "$GeneName|$Entrez";
	}elsif($id=~/^rna/){
		$rnaparent{$id} = $Parent;
		$location{$id} = "$chr:$start:$end";
	}
	if (exists $seq{$id}){
		if ($info =~/Name=([^;,]*)/) {
        	$GeneName = $1;
        }else{
        	$GeneName = '--';
        }
        if ($info=~/GeneID:([^;,]*)/) {
        	$Entrez = $1;
        }else{
        	$Entrez = '--';
        }
        $location{$id} = "$chr:$start:$end";
        $genelist{$id} = "$GeneName|$Entrez";
	}
}

foreach(sort keys %seq){
	if ($_=~/rna-/) {
		print OUT ">$_:$genelist{$rnaparent{$_}}";
	}else{
		print OUT ">$_:$genelist{$_}";
	}
	my $genebank = join ";",keys %{ $Genebanks{$_} };
	my $trans = join ";",keys %{ $Transcs{$_} };
	my $protein = join ";",keys %{ $Proteins{$_} };
	$genebank = '--' if ($genebank eq ""||$genebank=~/^$/);
	$trans = '--' if ($trans eq ""||$trans=~/^$/);
	$protein = '--' if ($protein eq ""||$protein=~/^$/);
	print OUT "|$genebank|$trans|$protein";
	if ($_=~/rna-/) {
		print OUT ":$location{$rnaparent{$_}}";
	}else{
		print OUT ":$location{$_}";
	}
	print OUT "\n$seq{$_}\n";
}
close IN;
close OUT;
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub ABSOLUTE_DIR #$pavfile=ABSOLUTE_DIR($pavfile);
{
	my $cur_dir=`pwd`;chomp($cur_dir);
	my ($in)=@_;
	my $return="";
	if(-f $in){
		my $dir=dirname($in);
		my $file=basename($in);
		chdir $dir;$dir=`pwd`;chomp $dir;
		$return="$dir/$file";
	}elsif(-d $in){
		chdir $in;$return=`pwd`;chomp $return;
	}else{
		warn "Warning! just for existing file and dir \n$in";
		exit;
	}
	chdir $cur_dir;
	return $return;
}
sub USAGE {#
        my $usage=<<"USAGE";
Contact:        tong.wang\@majorbio.com;
Script:		$Script
Description:
	enrich pro.fa from gffread by gff file
eg:
	perl $Script -i ref.temp.pro.fa -g Genome.gff -o ref.pro.fa

Usage:
  Options:
  -i	<file>	input ref.temp.pro.fa
  -g	<file>	input genome gff file
  -o	<file>	output file

  -h         Help

USAGE
        print $usage;
        exit;
}

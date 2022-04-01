#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($falist,$in,$out);        #1
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"falist:s"=>\$falist,
	"anno:s"=>\$in,     #2
	"o:s"=>\$out,   #3
			) or &USAGE;
&USAGE unless ($falist and $in and $out); #4
#######################################################################################
$falist = ABSOLUTE_DIR($falist);
$in = ABSOLUTE_DIR($in);

#get gene name list; geneID-->genetitle
my %gene;
open IN,$falist;
while (<IN>){
	chomp;
	next if ($_ eq ""||/^$/);
	s/\s//g;
	$/=">";
	open Seq,$_;
	while (<Seq>) {
		chomp;
		next if ($_ eq ""||/^$/);
		my ($id,undef) = split /\s+/,$_;
		my ($name,undef) = split /\:/,$id;
		$gene{$name} = $id;
	}
	close Seq;
	$/="\n";
}
close IN;

#merge interProscan anno result
#%interProanno; geneID-->gene_annos
#%interProacc; geneID-->gene_accessions
my (%interProanno,%interProacc);
my @outanno = glob("$in/*.interPro.interproscan"); #interpro result suffix
foreach my $peranno (@outanno){
	open IN,$peranno;
	my $temp='';
	while (<IN>){
		chomp;
		next if ($_ eq ""||/^$/||/^#/);
		my ($query,$MD5,$seqlength,$analysis,$signatureAccession,$signatureDescription,$start,$stop,$score,$status,$date,$ipaccession,$ipanno,$Go,$Path) = split /\t/,$_;
		my ($id,undef) = split /\:/,$query;
		next if ($analysis=~/Pfam/);
		#maybe blank,  default '--'
		$signatureDescription = '--' if (!defined $signatureDescription||$signatureDescription eq "");
		$ipaccession = '--' if (!defined $ipaccession||$ipaccession eq "");
		$ipanno = '--' if (!defined $ipanno||$ipanno eq "");
		$Go = '--' if (!defined $Go||$Go eq "");
		$Path = '--' if (!defined $Path||$Path eq "");
		if (!exists $interProanno{$id}){
			$interProacc{$id} = $ipaccession;
			$interProanno{$id} = "$ipanno|$analysis|$signatureAccession|$signatureDescription";
		}elsif("$id,$analysis,$signatureAccession" ne $temp){
			$interProacc{$id} .= ",$ipaccession";
			$interProanno{$id} .= ":$ipanno|$analysis|$signatureAccession|$signatureDescription";
		}
		$temp = "$id,$analysis,$signatureAccession";
	}
	close IN;
}

#print result
open OUT,">$out";
print OUT "#GeneID\tInterProAccession\tAnno\n";
foreach (sort keys %gene){
	if (!exists $interProacc{$_}){
		print OUT "$gene{$_}\t--\t--\n";
	}else{
		print OUT "$gene{$_}\t$interProacc{$_}\t$interProanno{$_}\n";
	}
}
close OUT;
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

sub ABSOLUTE_DIR #$pavfile=&ABSOLUTE_DIR($pavfile);
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
		warn "Warning just for file and dir \n$in";
		exit;
	}
	chdir $cur_dir;
	return $return;
}

sub USAGE {           #5
        my $usage=<<"USAGE";
Contact:	tong.wang\@majorbio.com
Version:	$version
Script:		$Script
Description:	Merge pfam result
Usage:

  -falist	<file>	input falist file
  -anno		<dir>	input pfam result dirname 
  -o		<file>	output file 
  -h		Help

USAGE
        print $usage;
        exit;
}

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

#merge Pfam anno result
#%pfamanno; geneID-->gene_annos
#%pfamacc; geneID-->gene_accessions
my (%pfamanno,%pfamacc);
my @outanno = glob("$in/*.pfam.anno"); #pfam result suffix
foreach my $peranno (@outanno){
	open IN,$peranno;
	while (<IN>){
		chomp;
		next if ($_ eq ""||/^$/||/^#/);
		my ($seq_id,$align_start,$align_end,$enve_start,$enve_end,$hmm_acc,$hmm_name,$type,$hmm_start,$hmm_end,$hmm_length,$bit_score,$E_value,$significance,$clan,$predicted_active_site,$strand,$nt_start,$nt_end) = split /\s+/,$_;
		my ($id,undef) = split /\:/,$seq_id;
		if (!exists $pfamanno{$id}){
			$pfamacc{$id} = $hmm_acc;
			$pfamanno{$id} = "$hmm_name|$type|$clan|$predicted_active_site";
		}else{
			$pfamacc{$id} .= ";$hmm_acc";
			$pfamanno{$id} .= ";$hmm_name|$type|$clan|$predicted_active_site";
		}
	}
	close IN;
}

#print result
open OUT,">$out";
print OUT "#GeneID\tPfamAccession\tPfamAnno\n";
foreach (sort keys %gene){
	if (!exists $pfamacc{$_}){
		print OUT "$gene{$_}\t--\t--\n";
	}else{
		print OUT "$gene{$_}\t$pfamacc{$_}\t$pfamanno{$_}\n";
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

#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($falist,$in,$out);        #1
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.1.0";
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

#get seq name list;
#%seq; $seq{seqID} = seqtitle
#title --> geneID:rnaID:GeneName|Entrez|Genebank|Transcript_id|Protein_id:loc:start:end
my %seq;
#my $key;
open IN,$falist;
while (<IN>){
	chomp;
	next if ($_ eq ""||/^$/);
	s/\s//g;
	$/=">";
	open SEQ,$_;
	while (<SEQ>) {
		chomp;
		next if (/^\s*$/);
		my ($title,undef) = split /\n/,$_;
		my ($id, undef) = split /\:/,$title; ##keyword rnaID
		$seq{$id} = $title;
	}
	close SEQ;
	$/="\n";
}
close IN;
#foreach $key (keys %seq)
#{
#        print "$key=>$seq{$key}\n";
#}
#merge interProscan anno result
#%interProanno; $interProanno{seqID} = seq_annos
#%interProacc;  $interProacc{seqID}  = seq_accessions
my (%interProanno,%interProacc,%temp);
my @outanno = glob("$in/*.interPro.interproscan"); #interpro result suffix
foreach my $peranno (@outanno){
	open IN,$peranno;
	while (<IN>){
		chomp;
		next if ($_ eq ""||/^$/||/^#/);
		my ($query,$MD5,$seqlength,$analysis,$signatureAccession,$signatureDescription,$start,$stop,$score,$status,$date,$ipaccession,$ipanno,$Go,$Path) = split /\t/,$_;
		my ($id, undef) = split /\:/,$query; ##keyword rnaID
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
		}elsif(!exists $temp{"$id,$analysis,$signatureAccession"}){
			$interProacc{$id} .= ",$ipaccession";
			$interProanno{$id} .= ":$ipanno|$analysis|$signatureAccession|$signatureDescription";
		}
		$temp{"$id,$analysis,$signatureAccession"}++;
	}
	close IN;
}

#print result
open OUT,">$out";
print OUT "#SeqID\tInterProAccession\tInterProAnno\n";
foreach (sort keys %seq){
	if (!exists $interProacc{$_}){
		print OUT "$seq{$_}\t--\t--\n";
	}else{
		print OUT "$seq{$_}\t$interProacc{$_}\t$interProanno{$_}\n";
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
		warn "Warning! just for existing file and dir.\n$in";
		exit;
	}
	chdir $cur_dir;
	return $return;
}

sub USAGE {           #5
        my $usage=<<"USAGE";
Contact:	tong.wang\@majorbio.com
Version:	$version
Time:		20181128
Script:		$Script
Description:	Merge InterProScan result
Usage:

  -falist	<file>	input falist file
  -anno		<dir>	input InterProScan result dirname; file = *.InterPro.interproscan
  -o		<file>	output file 
  -h		Help

USAGE
        print $usage;
        exit;
}

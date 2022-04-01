#!/mnt/bin/env perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($ref,$gff,$gene,$out);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
        "help|?" =>\&USAGE,
        "ref:s"=>\$ref,
        "gff:s"=>\$gff,
		"out:s"=>\$out,
		"gene:s"=>\$gene,
                        ) or &USAGE;
&USAGE unless ($ref and $gff and $gene);
my %stat;
open Ref,$ref;
if($ref=~/gz$/){
        close Ref;
        open Ref,"gunzip -c $ref|";
}
$/=">";
while(<Ref>){
        chomp;
        next if ($_ eq ""|| /^$/);
        my($chr,@seq)=split(/\s+/,$_);
        my $seq=join("",@seq);
        $stat{$chr}=$seq;
}
close Ref;

my(%gene,%tran,%info,$gene_type);
open Gff,$gff;
if($gff=~/gz$/){
	close Gff;
	open Gff,"gunzip -c $gff|";
}
$/="\n";
while (<Gff>){
	chomp;
	next if ($_ eq "" || /^$/);#chr1	RefSeq	gene	3631	5899	.	+	.	ID=gene0;Dbxref=Araport:AT1G01010,TAIR:AT1G01010,GeneID:839580;Name=NAC001;gbkey=Gene;gene=NAC001;gene_biotype=protein_coding;gene_synonym=ANAC001,NAC	domain	containing	protein	1,T25K16.1,T25K16_1;locus_tag=AT1G01010
	my $line=$_ ;
	my($chr,$source,$type,$start,$end,$score,$strand,$phase,$info,@unkon)=split(/\t/,$line);

  if($info=~/ID=($gene)\;(.*)/){
    if($type eq "gene"){
      $gene_type="gene";  # 确定传进来的gene是不是一个gene_id
    }else{
      $gene_type="other";
    }
    last;
  }
}
close Gff;
print $gene_type, "\n";

open Gff,$gff;
if($gff=~/gz$/){
	close Gff;
	open Gff,"gunzip -c $gff|";
}
while (<Gff>){
	chomp;
	next if ($_ eq "" || /^$/);#chr1	RefSeq	gene	3631	5899	.	+	.	ID=gene0;Dbxref=Araport:AT1G01010,TAIR:AT1G01010,GeneID:839580;Name=NAC001;gbkey=Gene;gene=NAC001;gene_biotype=protein_coding;gene_synonym=ANAC001,NAC	domain	containing	protein	1,T25K16.1,T25K16_1;locus_tag=AT1G01010
	my $line=$_ ;
	my($chr,$source,$type,$start,$end,$score,$strand,$phase,$info,@unkon)=split(/\t/,$line);
	my ($geneid,$transr,$transe,$exonid);

  if($gene_type eq "gene"){
    if($type eq "mRNA"){
  		# $geneid=$2 if($info=~/(.*);Parent=(\w+)\;(.*)/);
      # $transr=$1 if($info=~/ID=(\w+)\;(.*)/);
      $geneid=$2 if($info=~/(.*);Parent=(.*?)\;(.*)/);
  		$transr=$1 if($info=~/ID=(.*?)\;(.*)/);
  		if($geneid eq $gene){
  			$gene{$transr}=$geneid;
  		}
  		next;
  	}
  	if($type eq "exon"){
  		# $exonid=$1 if($info=~/ID=(\w+)\;(.*)/);
  		# $transe=$2 if($info=~/(.*);Parent=(\w+)\;(.*)/);
      $exonid=$1 if($info=~/ID=(.*?)\;(.*)/);
      $transe=$2 if($info=~/(.*);Parent=(.*)/);
      my @item = split ";", $transe;
      $transe = $item[0];
  		$tran{$exonid}=$transe;
  		$info{$exonid}=join("\,",$exonid,$start,$end,$chr,$transe,$strand,$type);
      next;
  	}
  	if($type eq "CDS"){
  		# $transe=$2 if($info=~/(.*);Parent=(\w+)\;(.*)/);
  		# $exonid=$1 if($info=~/ID=(\w+)\;(.*)/);
      $exonid=$1 if($info=~/ID=(.*?)\;(.*)/);
      $transe=$2 if($info=~/(.*);Parent=(.*)/);
      my @item = split ";", $transe;
      $transe = $item[0];
  		$tran{$exonid}=$transe;
  		$info{$exonid}=join("\,",$exonid,$start,$end,$chr,$transe,$strand,$type);
      next;
  	}
  }else{
    # if($info=~/(.*)=($gene)\;(.*)/){
    if($info=~/(.*)=($gene)(.*)/){
      if($type eq "mRNA"){
    		# $transr=$1 if($info=~/ID=(\w+)\;(.*)/);
        $transr=$1 if($info=~/ID=(.*?)\;(.*)/);
    		$gene{$transr}=$gene;
    		next;
    	}
    	if($type eq "exon"){
    		# $exonid=$1 if($info=~/ID=(\w+)\;(.*)/);
    		# $transe=$2 if($info=~/(.*);Parent=(\w+)\;(.*)/);
        $exonid=$1 if($info=~/ID=(.*?)\;(.*)/);
        $transe=$2 if($info=~/(.*);Parent=(.*)/);
        my @item = split ";", $transe;
        $transe = $item[0];
    		$tran{$exonid}=$transe;
    		$info{$exonid}=join("\,",$exonid,$start,$end,$chr,$transe,$strand,$type);
        next;
    	}
    	if($type eq "CDS"){
    		# $transe=$2 if($info=~/(.*);Parent=(\w+)\;(.*)/);
    		# $exonid=$1 if($info=~/ID=(\w+)\;(.*)/);
    		$exonid=$1 if($info=~/ID=(.*?)\;(.*)/);
        $transe=$2 if($info=~/(.*);Parent=(.*)/);
        my @item = split ";", $transe;
        $transe = $item[0];
    		$tran{$exonid}=$transe;
    		$info{$exonid}=join("\,",$exonid,$start,$end,$chr,$transe,$strand,$type);
        next;
    	}
    }
  }
}
close Gff;

open Out,">$out/$gene.exon.fa";
foreach my $transr(sort keys %gene){
	if($gene{$transr} eq $gene){
		foreach my $exonid(sort keys %info){
			my($exon,$start,$end,$chr,$transe,$strand,$type)=split("\,",$info{$exonid});
			if($transe eq $transr){
				my $length= $end - $start;
				my $exonseq=substr($stat{$chr},$start - 1,$length);
				if($strand=~/\+/){
					print Out ">$transe\t$exon\t$start\t$end\t$type\n$exonseq\n";
				}else{
					$exonseq=~s/a|A/T/g;
					$exonseq=~s/t|T/A/g;
					$exonseq=~s/g|G/C/g;
					$exonseq=~s/c|C/G/g;
					my @a;
					for(my $i = 0; $i < length($exonseq); $i++){
						$a[$i] = substr($exonseq,$i,1);
					}
					@a=reverse@a;
					print Out ">$chr\t$transe\t$exon\t$start\t$end\t$type\n";
					print Out join("",@a),"\n";
				}
			}
		}
	}
}
close Out;
#######################################################################################
#print STDOUT "\nDone. Total elapsed time : "d,time()-$BEGIN_TIME,"s\n";
########################################################################################
sub USAGE {#
        my $usage=<<"USAGE";
Contact:        caixia.tian\@majorbio.com;
Script:                 $Script
Description:
        get gene's all exon seq
        eg:
        perl $Script -ref -gff -out -gene

Usage:
  Options:
	-ref	<file>  input ref.fa
	-gff	<file>  input ref.gff
	-out	<file>	output file name
	-gene	<str>	geneid
	-h			Help

USAGE
        print $usage;
        exit;
}

#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($bam,$id,$fai,$out);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
    "bam:s"=>\$bam,
    "id:s"=>\$id,
    "fai:s"=>\$fai,
	"out:s"=>\$out,
			) or &USAGE;
&USAGE unless ($bam and $fai and $id and $out);
mkdir $out if (!-d $out);
mkdir "$out/combine" if (!-d "$out/combine");
open IN,"$fai";
my %hash;
while(<IN>){
    chomp;
    next if ($_ eq ""||/^$/);
    my ($id,$len,@others)=split(/\s+/,$_);
    $hash{$id}=$len;
    }
close IN;
open IN,"samtools view $bam|";
#open IN,"$bam";
#my %region;
my %fq;
my %stat;
my %filehand;
while(<IN>){
    $_=~s/[\n\r]//g;
    my ($fqid,$flags,$chr1,$pos1,$a,$b,$chr2,$pos2,$c,$seq,$qual,@others)=split/\t/,$_;
    next if $seq eq "*";
        my $rpos1=join("-",$chr1,$pos1-500,$pos1+length($seq)+500);
        if($chr2 eq "="){$chr2=$chr1;}
		my $rpos2=join("-",$chr2,$pos2-500,$pos2+length($seq)+500);
		#my $rpos=join(".",(split/\-/,$rpos1));
        my $rpos=$rpos1;
		#$rpos=join(".",(split/\-/,$rpos2)) if($rpos1 =~ /Chrinsert/);
        #next if $rpos1 =~ /Chrinsert/ and $rpos1 =~/=/;
        $rpos=$rpos2 if($rpos1 =~ /Chrinsert/ and $rpos1 !~/=/);

        my $flag=0;
#        my %filehand;
	if (/Chrinsert/) {
		if ($rpos1 =~ /Chrinsert/ && $rpos2 =~ /Chrinsert/){
			if (!exists $filehand{"both"}) {
				open $filehand{"both"},">>$out/$id.onlyinsert.fq";
				print {$filehand{"both"}} "\@$fqid\n$seq\n+\n$qual\n";
			}
		}else{
			foreach my $region (sort keys %filehand) {
				my ($c1,$p1,$p2)=split(/\-/,$region);
				my ($c2,$p3,$p4)=split(/\-/,$rpos);
				next if ($c1 ne $c2);
				if (!($p1 > $p4) && !($p2 < $p3)) {
					$flag=1;
					$stat{$region}{readnum}++;
					$stat{$region}{basenum}+=length($seq);
					$stat{$region}{q30}+=q30($qual);
					print {$filehand{$region}} "\@$fqid\n$seq\n+\n$qual\n";
				}
			}
			if ($flag == 0) {
				open $filehand{$rpos},">>$out/$id.$rpos.fq";
				$stat{$rpos}{readnum}++;
				$stat{$rpos}{basenum}+=length($seq);
				$stat{$rpos}{q30}+=q30($qual);
				print {$filehand{$rpos}} "\@$fqid\n$seq\n+\n$qual\n";
			}
		}
	}else{
		foreach my $region (sort keys %filehand) {
			my ($c1,$p1,$p2)=split(/\-/,$region);
			my ($c2,$p3,$p4)=split(/\-/,$rpos);
			next if ($c1 ne $c2);
            my $flag=0;
			if (!($p1 > $p4) && !($p2 < $p3)) {
				$flag=1;
				$stat{$region}{readnum}++;
				$stat{$region}{basenum}+=length($seq);
				$stat{$region}{q30}+=q30($qual);
				print {$filehand{$region}} "\@$fqid\n$seq\n+\n$qual\n";
			}
		}
	}

}
close IN;
my %merge;
foreach my $region (sort keys %filehand) {
	my ($c1,$p1,$p2)=split(/-/,$region);
	my $flag=0;
	foreach my $merge (sort keys %merge) {
		my ($c2,$p3,$p4)=split(/-/,$merge);
		next if ($c1 ne $c2);
		if (!($p1 > $p4) && !($p2 < $p3)) {
			$flag=1;
			my $news=min($p1,$p2,$p3,$p4);
			my $newe=max($p1,$p2,$p3,$p4);
			push @{$merge{join("-",$c1,$news,$newe)}},@{$merge{$region}};
			push @{$merge{join("-",$c1,$news,$newe)}},@{$merge{$region}};
			delete $merge{$region};
			delete $merge{$region};
		}
	}
	if ($flag ==1 ) {
		push @{$merge{$region}},"$out/$id.$region.fq";
	}
}
foreach my $region (sort keys %merge) {
	my $str=join(" ",@{$merge{$region}});
	`cat $str > $out/$region.fq && rm -f $str`;
}
open ST,">$out/stat.xls";
foreach my $regions (%stat){
    print ST "$id\t$regions\t$stat{$regions}{readnum}\t$stat{$regions}{basenum}\t$stat{$regions}{q30}\n";
}
sub q30{
    my @qstr=split//,$_;
    my $q20;
    my $q30;
    my $num=scalar @qstr;
    for(my $i=0;$i<=$#qstr;$i++){
        my $qual = ord($qstr[$i])-33;
         if($qual >=30){
             $q30++;
             $q20++;
             }
         if($qual >=20){
             $q20++;
             }
        }
        #my $q30s=sprintf("%.4f",$q30/$num);
        my $q30s=$q30/$num;
     return $q30s;
    }
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub USAGE {#
        my $usage=<<"USAGE";
Contact:        nanshan.yang\@majorbio.com;
Script:			$Script
Description:
	
	eg:
	perl $Script -g -i -o 

Usage:
  Options:
  -bam    <file>  input bam file
  -id   bam id 
  -fai  fai file
  -out	<dir>	out dir     
  -h         Help

USAGE
        print $usage;
        exit;
}

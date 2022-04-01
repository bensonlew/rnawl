#!/usr/bin/perl -w
use strict;
use warnings;
use Carp;


use Getopt::Long qw(:config no_ignore_case bundling);

my %opts;
GetOptions(\%opts,"i=s","o=s","r=s","m=s", "n=s", "nf=s");


my $usage = <<__EOUSAGE__;

#################################################################################### 
# Usage:
#	perl $0 --i ref.bed --r ref.des
# Required:
#
#  --i		bed file
#  --r          description file
#  --m          gene2trans2bed file
#  --n          tran2name file
#  --nf         where gene name from  
#  --o	        biomart file
#
####################################################################################


__EOUSAGE__

    ;


die $usage  if( !$opts{i} || !$opts{m});

my  $output= $opts{i}."biomart";
my $source = "unknown";
if (exists $opts{nf}){
    $source = $opts{nf};
}

my %des;
if($opts{r}){
    open (DES, $opts{r}) || die "Can not open $opts{r}\n";
    while(<DES>){
        chomp;
        my ($gene, $name) = split(/\t/, $_);
        $name=~s/[\r\n]//g;
        $des{$gene} = $name;
    }
    close DES;
}



my %t2g =();
my %t2p =();

open (G2I, $opts{m}) || die "Can not open $opts{m}\n";
while(<G2I>){
    chomp;
    my ($gene, $trans, $pep) = split(/\t/, $_);
    $t2g{$trans} = $gene;
    $t2p{$trans} = $pep;
}
close G2I;

my %t2n =();
if($opts{n}){
    open (T2N, $opts{n}) || die "Can not open $opts{n}\n";
    while(<T2N>){
        chomp;
        my ($tran, $name) = split(/\t/, $_);
        $name=~s/[\r\n]//g;
        $t2n{$tran} = $name;
    }
    close T2N;
}


open (INFILE,$opts{i}) || die "Can not open $opts{i}\n";
open (OUTFILE, "> $opts{o}") || die "Can not open $opts{o}\n";

while (<INFILE>) {
    chomp;
    my @line = split(/\t/,$_);
    my $strand = 0;
    my $start = $line[1];
    if ($line[5] eq "+"){
        $strand = 1;
    }elsif($line[5] eq "-"){
        $strand = -1;
        $start = $line[2];
    }
    my $len = 0;
    my @exon_len = split(",", $line[10]);
    foreach(@exon_len){
        $len += $_;
    }
    my $gene = $line[3];
    my $pep = "";
    if(exists $t2g{$line[3]}){
        $gene = $t2g{$line[3]};
    }
    if(exists $t2p{$line[3]}){
        $pep = $t2p{$line[3]};
    }
    print OUTFILE $gene."\t".$line[3]."\t";
    if(exists $opts{n}){
        print OUTFILE $t2n{$line[3]}."\t".$source."\t";
    }else{
        # print "\t\t";
    }

    print OUTFILE $pep."\t";

    if(exists $des{$line[3]}){
        print OUTFILE $des{$line[3]}."\t";
    }else{
        print OUTFILE "\t";
    }
    print OUTFILE $line[0]."\t".$line[1]."\t".$line[2]."\t".$strand."\t".$line[1]."\t".$line[2]."\t".$start."\t".$len."\tprotein_coding\tprotein_coding\n";
}

close INFILE;
close OUTFILE;

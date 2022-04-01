#!/usr/bin/perl -w
use strict;
use warnings;
use Carp;


use Getopt::Long qw(:config no_ignore_case bundling);

my %opts;
GetOptions(\%opts,"i=s","o=s","r=s");


my $usage = <<__EOUSAGE__;

#################################################################################### 
# Usage:
#	perl $0 --i ref.bed --r ref.des
# Required:
#
#  --i		bed file
#  --r          description file
#  --o	        biomart file
#
####################################################################################


__EOUSAGE__

    ;


die $usage  if( !$opts{i} || !$opts{r} );

my  $output= $opts{i}."biomart";
my %des;
open (DES, $opts{r}) || die "Can not open $opts{r}\n";
while(<DES>){
    chomp;
    my ($gene, $name) = split(/\t/, $_);
    $des{$gene} = $name;
}
close DES;

open (INFILE,$opts{i}) || die "Can not open $opts{i}\n";

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
    print $line[3]."\t".$line[3]."\t".$line[3]."\t";
    if(exists $des{$line[3]}){
            print $des{$line[3]}."\t";
    }else{
        print "\t";
    }
    print $line[0]."\t".$line[1]."\t".$line[2]."\t".$strand."\t".$line[1]."\t".$line[2]."\t".$start."\t".$len."\tprotein_coding\tprotein_coding\n";
}

#! /usr/bin/perl
use strict
use warnings

if (@ARGV != 5){
   print "perl cal_coverage fasta1 fasta2 fastas ref_fasta coverage;
   exit;
}

my ($fasta_1, $fasta_2, $fasta_s, $fasta, $coverage)=@ARGV;
open(FASTA1, "<$fastq_1") or die;
open(FASTA2, "<$fastq_2") or die;
open(FASTAS, "<$fastq_s") or die;
open(FASTA, "<$fasta") or die;
open(COVERAGE, ">$coverage") or die;
print COVERAGE "Coverage"
while(<FASTA1>){
    chomp;
    if (/>/){
        my $gene = $_;
    }else{
        my $gene1_len;
        my $len = length($_)
    };
    $gene1_len += $len ;
}
while(<FASTA2>){
    chomp;
    if (/>/){
        my $gene = $_;
    }else{
        my $gene2_len;
        my $len = length($_)
    }
    $gene2_len += $len
}
while(<FASTAS>){
    chomp;
    if (/>/){
        my $gene = $_;
    }else{
        my $genes_len;
        my $len = length($_)
    }
    $genes_len += $len
}
while(<FASTA>){
    chomp;
    if (/>/){
        my $gene = $_;
    }else{
        my $ref_len;
        my $len = length($_)
    }
    $ref_len += $len
}

sub depth_calcu
{
    my $total_len;
        $total_len = $gene1_len + $gene2_len + $genes_len
    my $dep = $total_len/$ref_len
    print COVERAGE $dep, "\n";
}


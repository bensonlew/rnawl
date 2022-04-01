#!/usr/bin/perl -w
use strict;
use warnings;
use Bio::SeqIO;
use List::Util qw/max min/;

die "perl $0 <ass><gtf><cut><out><type>\n" unless  @ARGV==5;

my ($ass,$gene,$cut,$out,$type)=@ARGV;
my ($seq_name,$s,$e,$as_name,$last,$string,$end,$start,$seq,$ass_name,$location);
my %seq;

my $in  = Bio::SeqIO->new(-file => "$ass",
                        -format => 'Fasta');

while ( my $sq = $in->next_seq() ) {
    $seq{$sq->id} = uc($sq->seq);
}

open IN,$gene||die;
open OUT, "> $out" || die "can not open output files for write";
my $min = 0;
my $max = 0;
my %gene = ();

if($type eq "rock_index"){
    chomp;
    while(<IN>){
        my @line = split(/\t/, $_);
        my $id = $line[6];
        my $strand = $line[2];
        my @pos = split /\.\./, $line[1];
        my $start = $pos[0];
        my $end = $pos[1];
        my $location = $line[0];
        if($strand eq "-"){
            $start = $pos[1];
            $end = $pos[0];
        }
        if ($location=~/^p/){
            $as_name = $location =~ s/p/Plasmid/r;
        }elsif ($location=~/^Chr/) {
            $as_name = $location =~ s/Chr/Chromosome/r;
        }else{
            $as_name = $location;
        }
        if (exists($seq{$location})){
             $as_name = $location;
        }else{
        }

        $last = $max;
        if ($strand eq "-"){
            $string = "-";
            if($start >= length($seq{$as_name})){
                next;
            }
            $start = min ($start + $cut,length($seq{$as_name}));
            $end = $end;
            $seq = reverse(substr($seq{$as_name},$end-1,$start-$end+1));
            $seq=~tr/AGCT/TCGA/;
        }else{
            $string = "+";
            if($start <= 1){
                next;
            }
            $start = max (1, $start - $cut);
            $end = $end;
            $seq = substr($seq{$as_name},$start-1,$end-$start+1);
        }
        $gene{$id} = 1;
        print OUT ">$id\t$string\t$start\t$end\t|$location\n$seq\n";
    }

}else{
    while(<IN>) {
        chomp;
        my @line = split(/\t/, $_);
        if (! ($line[8]) ){
            next;
        }
        if ( $line[8] =~ /gene_id "([^"]*)"/){
            my $id = $1;
            if (exists $gene{$id} && $gene{$id} == 1){
                next;
            }
            my $strand = $line[6];

            my $start = $line[3];
            my $end = $line[4];
            my $location = $line[0];
            if($strand eq "-"){
                $start = $line[4];
                $end = $line[3];
            }

            if ($location=~/^p/){
                $as_name = $location =~ s/p/Plasmid/r;
            }elsif ($location=~/^Chr/) {
                $as_name = $location =~ s/Chr/Chromosome/r;
            }else{
                $as_name = $location;
            }
            $last = $max;
            if ($strand eq "-"){
                $string = "-";
                if($start >= length($seq{$as_name})){
                    next;
                }
                $start = min ($start + $cut,length($seq{$as_name}));
                $end = $end;
                $seq = reverse(substr($seq{$as_name},$end-1,$start-$end+1));
                $seq=~tr/AGCT/TCGA/;
            }else{
                $string = "+";
                if($start <= 1){
                    next;
                }
                $start = max (1, $start - $cut);
                $end = $end;
                $seq = substr($seq{$as_name},$start-1,$end-$start+1);
            }
            $gene{$id} = 1;
            print OUT ">$id\t$string\t$start\t$end\t|$location\n$seq\n";
        }
    }


}


close IN;
close OUT;

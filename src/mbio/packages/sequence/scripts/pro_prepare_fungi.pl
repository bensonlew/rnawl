#!/usr/bin/perl -w
#guanqin.zou
use strict;
use warnings;
use List::Util qw/max min/;

die "perl $0 <ass><gene><cut><out>\n" unless  @ARGV==4;

my ($ass,$gene,$cut,$out)=@ARGV;
my ($seq_name,$s,$e,$as_name,$last,$string,$end,$start,$seq,$ass_name,$location);
my %seq;
open AS,$ass||die;
while(<AS>) {
  chomp;
  if ($_ =~ /^>(\S+)/) {
      $ass_name =$1;
  }else{
      $seq{$ass_name} = $_;
  }
}
close AS;
open IN,$gene||die;
open OUT, "> $out" || die "can not open output files for write";
my $min = 0;
my $max = 0;
while(<IN>) {
  chomp;
  if ($_ =~ /^>/) {
      my @line = split /\s+/, $_;
      $seq_name = $line[0];
      #$s = $line[4];
      $s = $line[1];
      #$e = $line[5];
      $e = $line[2];
      my @all = split /_/,$line[3];
      $location = $all[0];
      if ($location=~/^p/){
          $as_name = $location =~s/p/Plasmid/r;
      }elsif ($location=~/^Chr/) {
          $as_name = $location =~s/Chr/Chromosome/r;
      }else{
          $as_name = $location;
      }
      $last = $max;
      #if ($s > $e){
      if ($line[6] eq "-"){

          $string = "-";
          #$start = min ($s + $cut,length($seq{$as_name}));
          $start = min ($e + $cut,length($seq{$as_name}));
          #$end = $e;
          $end = $s;
          $seq = reverse(substr($seq{$as_name},$end-1,$start-$end+1));
          $seq=~tr/AGCT/TCGA/;
      }else{
          $string = "+";
          $start = max (1,$s - $cut);
          $end = $e;
          $seq = substr($seq{$as_name},$start-1,$end-$start+1);
      }
      $max = max ($s,$e);
      $min = min ($s,$e);
      if($min - $last >= $cut){
          print OUT "$seq_name\t$string\t$start\t$end\t|$location\n$seq\n";
      }
      }
  }
  close IN;
  close OUT;
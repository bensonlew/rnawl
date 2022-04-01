#! /usr/bin/perl
use strict;
use Bio::SeqIO;

if(@ARGV != 3)  {
    print "perl $0 Input gff sample\n";
    exit;
}
my ($input,$gff,$sample) = @ARGV;
open (IN,$gff) || die $!;
my %ass;
while(<IN>){
  chomp;
  my @temp =split /\t/;
  my @as =split /;/,$temp[7];
#if($as[0] eq "Name=16S_rRNA"){
if($as[0] =~ /16S/i){
  my @arry=split /\_/,$temp[1];
  if($arry[0] eq 'Chromosome'){
     $ass{'Chromosome'}=$temp[5] . "\t" . $temp[2] . "\t" . $temp[3];
     last;
}elsif($arry[0] =~ /Chromosome[1-9]*/){
    $ass{$arry[0]}=$temp[5] . "\t" . $temp[2] . "\t" . $temp[3];
    last;
}else{
    $ass{$arry[0]} = $temp[5] . "\t" . $temp[2] . "\t" . $temp[3];  ##
    }
    last;
}
}
close IN;
my $in = Bio::SeqIO->new(-file => "$input", -format => "fasta");
my $out = Bio::SeqIO->new(-file => ">16s.fasta", -format => "fasta");
while(my $seq = $in->next_seq())  {
    my $id = $seq->id;
    if(exists $ass{$id}){
       my ($stand,$start,$end)=split /\t/,$ass{$id};
       my $str;
       if($stand eq '-'){
           if($end < $start){
           $str = $seq->subseq($end, $start);
           }else{
           $str = $seq->subseq($start, $end);
           }
           $str=reverse($str);
           $str=~ tr/atcgATCG/tagcTAGC/;
       }else{
           $str = $seq->subseq($start, $end);
       }
       my $seq1 = Bio::Seq->new(-id =>$sample, -seq => $str);
     $out->write_seq($seq1);
}
}
#!/usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;

if(@ARGV != 5) {
    print "perl tiqu_16s.pl sample type Input summary out\n";
    exit;
}
my ($sample,$type,$input,$summary,$output) = @ARGV;
open (IN,$summary) || die $!;
my %hash;
while(<IN>){
   chomp;
   next if (/^Location/);
   my @temp =split /\t/;
   my $des;
   my $id;
   if($type eq "prephage"){
       $id =$temp[1];
       if ($temp[4] eq "+"){
           $des = $temp[0] . "\t" . $temp[2] . "\t" . $temp[3];
       }elsif($temp[4] eq "-"){
           $des = $temp[0] . "\t" . $temp[3] . "\t" . $temp[2];
       }
   }elsif($type eq "island"){
       $id =$temp[1];
       $des = $temp[0] . "\t" . $temp[2] . "\t" . $temp[3];
   }
   $hash{$id} = $des;
}
close IN;
my $in = Bio::SeqIO->new(-file => "$input", -format => "fasta");
my $out = Bio::SeqIO->new(-file => ">>$output", -format => "fasta");
while(my $seq = $in->next_seq())  {
    my $id =$seq->id;
    my $len = $seq->length;
    print "$id\t$len\n";
    foreach (sort keys %hash){
        my @temp =split /\t/,$hash{$_};
        if($id eq $temp[0])  {
           print "$temp[1]\t$temp[2]\n";
           my $str = $seq->subseq($temp[1], $temp[2]);
           my $seqname = $_ . "__" .$sample;
           my $seq1 = Bio::Seq->new(-id => $seqname, -seq => $str);
           $out->write_seq($seq1);
    }
  }
}

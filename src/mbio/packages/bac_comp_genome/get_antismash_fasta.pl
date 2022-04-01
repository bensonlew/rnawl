#!/usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;

if(@ARGV != 3) {
    print "perl $0 Input summary out\n";
    exit;
}
my ($input,$summary,$output) = @ARGV;
my $file = (split /\//,$summary)[-1];
my $sample =(split /.antismash_anno/,$file)[0];
open (IN,$summary) || die $!;
my %hash;
while(<IN>){
   chomp;
   next if (/^Cluster ID/);
   my @temp =split /\t/;
   my $id =$temp[0];
   my $des = $temp[1] . "\t" . $temp[3] . "\t" . $temp[4];
   $hash{$id} = $des;
}
close IN;
my $in = Bio::SeqIO->new(-file => "$input", -format => "fasta");
my $out = Bio::SeqIO->new(-file => ">>$output", -format => "fasta");
while(my $seq = $in->next_seq())  {
    my $id =$seq->id;
    foreach (sort keys %hash){
        my @temp =split /\t/,$hash{$_};
        if($id eq $temp[0])  {
           my $str = $seq->subseq($temp[1], $temp[2]);
           my $seqname = $_ . "__" .$sample;
           my $seq1 = Bio::Seq->new(-id => $seqname, -seq => $str);
           $out->write_seq($seq1);
    }
  }
}

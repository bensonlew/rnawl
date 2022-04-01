#!/usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;

if(@ARGV != 3)  {
    print "perl $0 Input bin_list prefix\n";
    exit;
}
my ($input,$binlist,$name) = @ARGV;
my $in = Bio::SeqIO->new(-file => "$input", -format => "fasta");
my $out = Bio::SeqIO->new(-file => ">>$name.fa", -format => "fasta");
my @bin_list = split /;/,$binlist;

while(my $seq = $in->next_seq()){
   my $id =$seq->id;
   for my $de (@bin_list){
       if ($de =~/(.*)\((.*)\)/){
           my ($start,$end)=split /,/,$2;
           if ($id eq $1){
              my $seq1 = $seq->subseq($start, $end);
              my $seq2 = Bio::Seq->new(-id => $id, -seq => $seq1);
              $out->write_seq($seq2);
           }
       }else{
           if ($id eq $de){
               my $seq1 = $seq->seq;
               my $seq2 = Bio::Seq->new(-id => $id, -seq => $seq1);
               $out->write_seq($seq2);
           }
       }
   }
}
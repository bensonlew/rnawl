#!/usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;

if(@ARGV != 2)  {
    print "perl $0 Input prefix\n";
    exit;
}
my ($input,$name) = @ARGV;
my $in = Bio::SeqIO->new(-file => "$input", -format => "fasta");
my $out = Bio::SeqIO->new(-file => ">>$name.large.scaf.fa", -format => "fasta");
my $out2 = Bio::SeqIO->new(-file => ">>$name.scaf.fa", -format => "fasta");
while(my $seq = $in->next_seq())  {
   my $id =$seq->id;
   my $str = $seq->seq;
   my $len =$seq->length;
   if($len >=655360){
       my $seq1 = Bio::Seq->new(-id => $id, -seq => $str);
       $out->write_seq($seq1);
   }else{
       my $seq2 = Bio::Seq->new(-id => $id, -seq => $str);
       $out2->write_seq($seq2);
   }
}
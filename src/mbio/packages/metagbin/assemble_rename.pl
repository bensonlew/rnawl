#!/usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;

if(@ARGV != 2)  {
    print "perl $0 Input sample_name\n";
    exit;
}
my ($input,$name) = @ARGV;
my $in = Bio::SeqIO->new(-file => "$input", -format => "fasta");
my $out = Bio::SeqIO->new(-file => ">>$name.scaf.fa", -format => "fasta");
while(my $seq = $in->next_seq())  {
   my $id =$name . "_" . $seq->id;
   my $str = $seq->seq;
   my $seq1 = Bio::Seq->new(-id => $id, -seq => $str);
   $out->write_seq($seq1);
}
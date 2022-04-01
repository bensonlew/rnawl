#! /usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;
##author: hao.gao

if(@ARGV!=2) {
    print STDERR "perl $0 seq output\n";
    exit;
}
my ($seq,$output)=@ARGV;

my $in = Bio::SeqIO->new(-file => $seq, -format => "fasta");

while(my $seq = $in->next_seq)  {
    if($seq->length >=300 )  {
       my $id=$seq->id;
       my $seqs =$seq->seq;
       my $seq1 =Bio::Seq->new(-id =>$id,-seq =>$seqs);
		  my $out=Bio::SeqIO->new(-file =>">>$output",-format =>"fasta");
		  $out->write_seq($seq1);
    }
}

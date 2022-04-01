#! /usr/bin/perl -W
use Bio::SeqIO;

my ($input,$output)=@ARGV;

my $in = Bio::SeqIO->new(-file => "$input", -format => "fasta");
my $out = Bio::SeqIO->new(-file => ">$output", -format => "fasta");
my  $seqs;
while(my $seq = $in->next_seq())  {
    my $sq =$seq->seq;
    $sq=~s/[nN]*//ig;
    $seqs .=$sq;
}
 my $seq1 = Bio::Seq->new(-id => 'test', -seq => $seqs);##这里test不能改
 $out->write_seq($seq1);

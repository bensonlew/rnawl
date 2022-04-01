#!perl -w
use strict;
use Bio::SeqIO;

if(@ARGV != 2) {
    print "perl $0 input output\n";
    exit;
}
my ($ref,$output) = @ARGV;
open (OUT,">$output") || die $!;
my $in = Bio::SeqIO->new(-file => "$ref", -format => "fasta");
while(my $seq = $in->next_seq()) {
    my $id = $seq->id;
    my  $sqes =$seq->seq;
   my $length =$seq->length;
    if($length >=1000){
      print OUT ">$id\n$sqes\n";
    }
   }

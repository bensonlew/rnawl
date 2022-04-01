#! /usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;
##author: hao.gao

if(@ARGV!=2) {
    print STDERR "perl $0 dir output\n";
    exit;
}
my ($dir,$output)=@ARGV;
open(OUT,">$output") or die;
my @file =glob "$dir/*.scafSeq";
my @arry;
my %hash;
for my $file (@file){
my $seq_number;
my $in = Bio::SeqIO->new(-file => $file, -format => "fasta");
my $total_sequence_length = 0;
my @x;
while(my $seq = $in->next_seq)  {
    if($seq->length >=300 )  {
        $seq_number++;
        push @x, $seq->length;
        $total_sequence_length += $seq->seq =~ tr/atcgnATCGN/atcgnATCGN/;
    }
}
@x = sort {$b <=> $a} @x;
my ($num,$length) = &nn0(50, $total_sequence_length, @x);
    $hash{$file}{$length}=$seq_number;
 # print "$file\t $length\t$seq_number\n";
}

foreach my $key (sort {$b cmp $a} keys %hash){
  foreach my $key2 (sort {$hash{$key}{$a} cmp $hash{$key}{$a}} keys %hash){
#print "$key\n";
    push @arry,$key;
}
}
for my $de(@arry){
   print OUT "$arry[0]\n";
   last;
}


sub nn0  {
    my $nnn = shift @_;
    my $total = shift @_;
    my @nnall = @_;
    my $nnsum = 0;
    my $nnn_num = 0;
    my $nnn_length = 0;
    while(@nnall > 0 && ($nnsum < $nnn/100*$total))  {
        my $x = pop @nnall;
        $nnsum += $x;
        $nnn_num++;
        $nnn_length = $x;
    }
    return $nnn_num,$nnn_length;
}

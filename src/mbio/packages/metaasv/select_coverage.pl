#! /usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;
#use POSIX::ceil;

if(@ARGV!=3) {
    print STDERR "perl $0 dir output\n";
    exit;
}
my ($dir,$output,$coverage)=@ARGV;
open(OUT,">$output") or die;
my @file =glob "$dir/*.fq";
my @arry;
my %hash;
for my $file (@file){
my $seq_number;
my $in = Bio::SeqIO->new(-file => $file, -format => "fastq");
my $total_sequence_length = 0;
my @x;
while(my $seq = $in->next_seq)  {
    $seq_number++;
    push @x, $seq->length;
}
@x = sort {$a <=> $b} @x;
my ($num,$length) = &nn0($coverage, $seq_number, @x);
print OUT "$file\t$length\n";
}
close OUT;

sub nn0  {
    my $nnn = shift @_;
    my $total = shift @_;
    my @nnall = @_;
    my $nnn_num = $nnn/100*$total;
    $nnn_num = printf("%.f", $nnn_num);
#    print("$nnn_num\n");
    my $nnn_length = $nnall[$nnn_num];
    return $nnn_num,$nnn_length;
}

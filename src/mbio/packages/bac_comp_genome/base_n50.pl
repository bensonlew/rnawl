#! /usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;
use Cwd 'abs_path';
##author: qingchen.zhang

if(@ARGV!=8) {
    print STDERR "perl $0 [fasta] [sample] [genome] [type] [cds_gff] [rrna_gff] [trna_gff] [outfile]\n";
    exit;
}
my ($input,$sample,$genome,$type,$cds,$rrna,$trna,$output)=@ARGV;

open (OUT,">$output") || die $!;
print OUT "Sample\tGenome\tN50\tTotal_length\tN_length\tGC\tGene_num\trRNA_num\ttRNA_num\n";
my @arry;
my %hash;
my $seq_number;
my @file;
my $abs_path = abs_path($input);

my $cds_num=0;
open(CDS,"<$cds") or die;
while(<CDS>){
    chomp;
    next if ($_ =~ /^Gene ID/);
    $cds_num++;
}
close CDS;

my $rrna_num=0;
open(RRNA,"<$rrna") or die;
while(<RRNA>){
    chomp;
    next if ($_ =~ /^Gene ID/);
    $rrna_num++;
}
close RRNA;

my $trna_num=0;
open(TRNA,"<$trna") or die;
while(<TRNA>){
    chomp;
    next if ($_ =~ /^Gene ID/);
    $trna_num++;
}
close TRNA;

my $in = Bio::SeqIO->new(-file => $input, -format => "fasta");
my $total_sequence_length = 0;
my $total_special_length = 0;
my $total_n_length = 0;
my $g_length = 0;
my $c_length = 0;
my @x;

while(my $seq = $in->next_seq){
    if($seq->length >=0 )  {
        $seq_number++;
        push @x, $seq->length;
        if ($type eq "draft"){
            $total_sequence_length += $seq->seq =~ tr/atcgnATCGN/atcgnATCGN/;
            $total_special_length += $seq->seq =~ /[^atcgnATCGN]/;
            $total_n_length += $seq->seq =~ tr/nN/nN/;
            $g_length += $seq->seq =~ tr/gG/gG/;
            $c_length += $seq->seq =~ tr/cC/cC/;
        }else{
            $total_sequence_length += $seq->seq =~ tr/atcgnATCGN/atcgnATCGN/;
            $total_special_length += $seq->seq =~ /[^atcgnATCGN]/;
            $total_n_length += $seq->seq =~ tr/nN/nN/;
            $g_length += $seq->seq =~ tr/gG/gG/;
            $c_length += $seq->seq =~ tr/cC/cC/;
        }

    }
}
@x = sort {$b <=> $a} @x;


my ($num,$length) = &nn0(50, $total_sequence_length, @x);

  # print "$file\t $length\n";
my $gc = &nn1($total_sequence_length,$g_length, $c_length);

#¼ÆËãn50
sub nn0 {
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

#¼ÆËãgcº¬Á¿
sub nn1 {
    my $total = shift @_;
    my $g_length = shift @_;
    my $c_length = shift @_;
    my $gc;
    if($total != 0) {
        $gc = ($g_length + $c_length) / $total;
    }
    return $gc;
}

if ($type eq "draft"){
    print OUT "$sample\t$genome\t$length\t$total_sequence_length\t$total_n_length\t$gc\t$cds_num\t$rrna_num\t$trna_num\n";
}else{
    print OUT "$sample\t$genome\t-\t$total_sequence_length\t$total_n_length\t$gc\t$cds_num\t$rrna_num\t$trna_num\n";
}
close OUT;
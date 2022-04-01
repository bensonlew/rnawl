#! /user/bin/perl
use strict;
use Bio::SeqIO;
#author:gaohao
########
##主要是将gff文件转化成phage_info.txt文件
#############

if(@ARGV!=3) {
    print STDERR "perl $0 file scaf output\n";
    exit;
}
my ($input,$scaf,$output)=@ARGV;
my %len=&scaff_len($scaf);
open (IN,$input) || die $!;
open (OUT,">$output") || die $!;
while(<IN>){
  chomp;
    next if(/^Gene id/);
  my @temp=split /\t/;
  my @seq=split /_ORF/,$temp[1];
  if(exists $len{$seq[0]}){
    print OUT "$seq[0]\t$len{$seq[0]}\t$temp[0]\t$temp[2]\t$temp[3]\n";
}
   }
close IN;
close OUT;

sub scaff_len{
my ($input) = @_;
my $in = Bio::SeqIO->new(-file => "$input", -format => "fasta");
my %hash;
while(my $seq = $in->next_seq())  {
      my $id=$seq->id;
    my $len= $seq->length;
$hash{$id}=$len;
}
return %hash;
}

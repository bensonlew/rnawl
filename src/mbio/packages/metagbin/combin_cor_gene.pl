#! /usr/bin/perl -W
use strict;
use Bio::SeqIO;
if(@ARGV != 5){
  print "perl $0 sample_name blastout faa cor_list output \n";
  exit;
}
my ($sample,$input,$fa,$list,$output)=@ARGV;
open (IN,$input) || die $!;
open (OUT,">$output") || die $!;
my %blast;
while(<IN>){
 chomp;
 my @temp =split /\t/;
  $blast{$temp[1]}=$temp[0];
}
close IN;
my %seq;
my $in = Bio::SeqIO->new(-file => "$fa", -format => "fasta");
while(my $seq = $in->next_seq())  {
    my $id = $seq->id;
    if(exists $blast{$id}){
      my $str = $seq->seq;
      my $seq1 = Bio::Seq->new(-id =>$blast{$id}, -seq => $str);
     $seq{$blast{$id}}=$str;
}
}
my @cor_list;
open (IN2,$list) || die $!;
while(<IN2>){
 chomp;
my @temp =split /\t/;
push @cor_list,$temp[0];
}
close IN2;
my $seqs;
for my $de (@cor_list){
  if(exists $seq{$de}){
   $seqs .=$seq{$de};
}
}
print OUT ">$sample\n";
print OUT "$seqs\n";

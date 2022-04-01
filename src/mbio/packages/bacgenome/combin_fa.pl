#! /usr/bin/perl -W
use strict;
use Bio::SeqIO;

if(@ARGV != 2)  {
    print "perl $0 Input blastout\n";
    exit;
}
my ($input,$blast) = @ARGV;
open (IN,$blast) || die $!;
open (OUT,">all.blast_16s.xls") || die $!;
my %ass;
while(<IN>){
  chomp;
  my @temp =split /\t/;
  my @ll =split(/\|/, $temp[1]);
  print($ll[0]);
  my @ll2 = split(/_rRNA/, $ll[0]);
  $ass{$temp[1]}=$ll2[0] . "|" . $ll[1];
  $temp[1] = $ll2[0] . "\t" . $ll[1];
  my $des =join("\t",@temp);
  print OUT "$des\n";
}
close IN;

my %done ;
my $in = Bio::SeqIO->new(-file => "$input", -format => "fasta");
my $out = Bio::SeqIO->new(-file => ">all.fasta", -format => "fasta");
while(my $seq = $in->next_seq())  {
    my $id = $seq->id;
    if(exists $ass{$id}){
        my $str=$seq->seq;
        my $dd = $ass{$id} ;
        my $seq1 = Bio::Seq->new(-id =>$dd, -seq => $str);
        $out->write_seq($seq1);
}
}


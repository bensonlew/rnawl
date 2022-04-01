#! /usr/bin/perl -W
use strict;
use Bio::SeqIO;
if(@ARGV != 5){
  print "perl $0 sample_name blastout faa gff output \n";
  exit;
}
my ($sample,$input,$fa,$gff2,$output)=@ARGV;
open (IN,$input) || die $!;
open (OUT,">$output") || die $!;
print OUT "Bin_id\tLocation\tstart\tEnd\tName\tIndentity\tCoverage\tProtein\n";
my %blast;
while(<IN>){
 chomp;
 my @temp =split /\t/;
 my $des= ($temp[9]-$temp[8])/$temp[3];
 $blast{$temp[1]}=$temp[0] . "\t" .  $temp[2] . "\t" . $des ;
}
close IN;
my %seq;
my $in = Bio::SeqIO->new(-file => "$fa", -format => "fasta");
while(my $seq = $in->next_seq())  {
    my $id = $seq->id;
    if(exists $blast{$id}){
      my $str = $seq->seq;
     $seq{$id}=$blast{$id} . "\t" . $str;
}
}
my %gff;
open (IN2,$gff2) || die $!;
my $name;
while(<IN2>){
 chomp;
 next if(/^Gene/);
 my @temp =split /\t/;
 my @ss =split /_ORF/,$temp[1];
$gff{$temp[0]}=$ss[0] . "\t" . $temp[2] . "\t" . $temp[3];
}
close IN2;
foreach (sort keys %seq){
  if(exists $gff{$_}){
 print OUT "$sample\t$gff{$_}\t$seq{$_}\n";
}
}

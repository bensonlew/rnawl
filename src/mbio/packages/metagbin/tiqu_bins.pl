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
      my ($ids)=$id=~/(.*)_1/;
     $seq{$ids}=$blast{$id} . "\t" . $str;
}
}
my %gff;
open (IN2,$gff2) || die $!;
my $name;
while(<IN2>){
 chomp;
 if(/^#\s+(.*)/){
    $name=$1;
    <IN2>;
    <IN2>;
}else{
  my @temp =split /\t/;
 if ($temp[3] eq "-"){
   $gff{$name . "_" .$temp[0]} =$temp[2] . "\t" . $temp[1];
 }else{
 $gff{$name . "_" .$temp[0]} =$temp[1] . "\t" . $temp[2];
}
}
}
close IN2;
foreach (sort keys %seq){
  if(exists $gff{$_}){
  my @ss =split /\_gene/,$_;
 print OUT "$sample\t$ss[0]\t$gff{$_}\t$seq{$_}\n";
}
}

#! /user/bin/perl
use strict;
#author:gaohao

if(@ARGV!=3) {
    print STDERR "perl $0 type file prefix\n";
    exit;
}
my ($type,$input,$prefix)=@ARGV;
open (IN,$input) || die $!;
my $file="file.txt";
open (OT,">$file") || die $!;
while(<IN>){
  chomp;
if(/^Sequence\s+'(.*)'\s+\(.*\)/){
  print OT "$1\t";
}
if(/^CRISPR\s+([0-9]+)\s+Range:\s+([0-9]+)\s+\-\s+([0-9]+)/){
   my $num=$1;
   my $start=$2;
   my $end=$3;
   print OT "$1\t$2\t$3,";
}
if(/^([0-9]+)\s+([ATGCN]+)\s+([ATGCN]+)\s+\[\s+([0-9]+),\s+([0-9]+)\s+\]/){
  print OT "$1\t$2\t$4\t$3\t$5,"
 }elsif(/^([0-9]+)\s+([ATGCN]+)/){
   print OT "$1\t$2\t-\t-\t-,";
}
if(/^Repeats:\s+([0-9]+)\s+Average Length:\s+([0-9]+)\s+Average Length:\s+([0-9]+)/){
   print OT "$1\t$2\t$3\n";
}
}
close IN;
close OT;
if($type eq 'uncomplete'){
  &uncomplete($file,$prefix);
}elsif($type eq 'complete'){
  &complete($file,$prefix);
}

sub uncomplete{
   my ($file,$prefix)=@_;
open (IN2,$file) || die $!;
open (OUT,">$prefix.summary.xls") || die $!;
open (OUT2,">$prefix.detail.xls") || die $!;
print OUT "Sample_chr\tCrisprs_num\tCrisprs_start\tCrisprs_end\tDR_number\tDR_Average_length\tSPA_Average_length\n";
print OUT2 "Sample_chr\tCrisprs_num\tpostion\tDR Sequence\tDR Length\tSPA Sequence\tSPA Length\n";
my $name;
my %hash;
my $seq;
while(<IN2>){
 chomp;
 if(/^[0-9]+/){
  $hash{$name} .=";" . $_; 
}else{
  my @temp=split /\t/;
  my $des=join("\t",@temp[1..$#temp]);
($seq,$name)=$temp[0]=~/(\D+)(\d+)/;  #Scaffold ¸ÄÎª
   $hash{$name}=$des;
}
}
close IN2;
foreach (sort {$a <=> $b}keys %hash){
  if($hash{$_}=~/;/){
  my @arry=split /;/,$hash{$_};
  for my $de (@arry){
  my @fi=split /,/,$de;
  my @head=split /\t/,$fi[0];
  my $des=$seq . $_;
  my $cre="CRISPR" . $head[0];
  for my $i (1..@fi-2){
 print OUT2 "$des\t$cre\t$fi[$i]\n";
  }
  print OUT "$des\t$cre\t$head[1]\t$head[2]\t$fi[-1]\n";
  }
}else{
  my @fi=split /,/,$hash{$_};
  my @head=split /\t/,$fi[0];
  my $des=$seq . $_;
  my $cre="CRISPR" . $head[0];
  for my $i (1..@fi-2){
   print OUT2 "$des\t$cre\t$fi[$i]\n";
  }
  print OUT "$des\t$cre\t$head[1]\t$head[2]\t$fi[-1]\n";
}
}
close OUT;
close OUT2;
#`rm $file`;
}

sub complete{
   my ($file,$prefix)=@_;
open (IN2,$file) || die $!;
open (OUT,">$prefix.summary.xls") || die $!;
open (OUT2,">$prefix.detail.xls") || die $!;
print OUT "Sample_chr\tCrisprs_num\tCrisprs_start\tCrisprs_end\tDR_number\tDR_Average_length\tSPA_Average_length\n";
print OUT2 "Sample_chr\tCrisprs_num\tpostion\tDR Sequence\tDR Length\tSPA Sequence\tSPA Length\n";
my $name;
my %hash;
my $seq;
while(<IN2>){
 chomp;
 if(/^[0-9]+/){
  $hash{$name} .=";" . $_; 
}else{
  my @temp=split /\t/;
  my $des=join("\t",@temp[1..$#temp]);
     $name=$temp[0];
   $hash{$name}=$des;
}
}
close IN2;
foreach (sort {$a <=> $b}keys %hash){
  if($hash{$_}=~/;/){
  my @arry=split /;/,$hash{$_};
  for my $de (@arry){
  my @fi=split /,/,$de;
  my @head=split /\t/,$fi[0];
  my $des=$seq . $_;
  my $cre="CRISPR" . $head[0];
  for my $i (1..@fi-2){
   print OUT2 "$des\t$cre\t$fi[$i]\n";
  }
  print OUT "$des\t$cre\t$head[1]\t$head[2]\t$fi[-1]\n";
  }
}else{
  my @fi=split /,/,$hash{$_};
  my @head=split /\t/,$fi[0];
  my $des=$seq . $_;
  my $cre="CRISPR" . $head[0];
  for my $i (1..@fi-2){
   print OUT2 "$des\t$cre\t$fi[$i]\n";
  }
  print OUT "$des\t$cre\t$head[1]\t$head[2]\t$fi[-1]\n";
}
}
close OUT;
close OUT2;
#`rm $file`;
}

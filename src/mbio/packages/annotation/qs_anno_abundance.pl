#! /usr/bin/perl -W
use strict;
use warnings;
use Getopt::Long;
my %opts;

GetOptions (\%opts,"i=s","o=s","Identity=i","Bit_score=i","QSclass=s","geneprofile=s");
my $usage = <<"USAGE";
       Program : $0
       Version : 1.0
       Discription: QS class annotation 
       Usage :perl $0 [options]
                   -i* input QS annotation result file
                   -o* output prefix
                   -geneprofile  abundance file of geneprofile
 
      eg: perl $0 -i GP_A.blastM8.xls -geneprofile /mnt/ilustre/users/hao.gao/Gaohao/MJ201508244003_zhouhuiping/out/geneProfile/gene_profile.reads_number.txt -o result
                 
USAGE
die $usage if ( !($opts{i}) || !( $opts{o}) ||  !($opts{geneprofile}));

#define defalts
#$opts{Identity}=$opts{Identity}?$opts{Identity}:50;
#$opts{Bit_score}=$opts{Bit_score}?$opts{Bit_score}:60;


open (IN,$opts{i}) || die $!;
open (OUT,">$opts{o}_class_profile.xls") || die $!;
open (OUT3,">$opts{o}_lowest_profile.xls") || die $!;

my %hash;
my %class;
my %gene;
while(<IN>){
  chomp;
  next if(/^#Query/);
 my @temp=split /\t/;
 if(exists $hash{$temp[1]}){
  $hash{$temp[1]} .=";" . $temp[0] ;
}else{
  $hash{$temp[1]} =$temp[0] ;
}
if($temp[2]=~/;/){
  my @tem=split /;/,$temp[2]; 
  for my $de (@tem){
   if(exists $class{$de}){
  $class{$de} .=";" . $temp[0] ;
}else{
  $class{$de} =$temp[0] ;
}
}
}else{
  if(exists $class{$temp[2]}){
  $class{$temp[2]} .=";" . $temp[0] ;
}else{
  $class{$temp[2]} =$temp[0] ;
}
}
  $gene{$temp[0]} =1;
}
close IN;

open (IN3,$opts{geneprofile}) ||die $!;
my %genepro;
my $name=<IN3>;
chomp($name);
my @name=split /\t/,$name;
shift @name;
my $des=join("\t",@name[0..$#name]);
print OUT "#Class\t$des\n";
print OUT3 "QS_id\t$des\n";
while(<IN3>){
  chomp;
  my @temp=split /\t/;
  my $adu=join("\t",@temp[1..$#temp]);
  $genepro{$temp[0]}=$adu;
}
close IN3;

my %hash5;
foreach (sort keys %gene){
  if(exists $genepro{$_}){
   $hash5{$_}=$genepro{$_};
}
}

foreach (sort keys %hash){  
 print OUT3 "$_";
 my $de=$hash{$_};
  my @arry3;
  if($de=~/;/){
     my @arry=split /;/,$de;
     for my $ds (@arry){
      if(exists $hash5{$ds}){
  my  @arry2=split /\t/,$hash5{$ds};
       push @arry3,[@arry2];
}
}
  }else{
  if(exists $hash5{$de}){
  my  @arry2=split /\t/,$hash5{$de};
       push @arry3,[@arry2];
}
}  
for my $j (0..$#name){
my $sum=0;
 for my $i(0..$#arry3){
  $sum +=$arry3[$i][$j];
}
  #$total +=$sum;
print OUT3 "\t$sum";
}
print OUT3 "\n";
}


foreach (sort keys %class){
   print OUT "$_";
 my $de=$class{$_};
  my @arry3;
 if($de=~/;/){
     my @arry=split /;/,$de;
     for my $ds (@arry){
      if(exists $hash5{$ds}){
  my  @arry2=split /\t/,$hash5{$ds};
       push @arry3,[@arry2];
}
}
  }else{
  if(exists $hash5{$de}){
     print "$hash5{$de}\n";
  my  @arry2=split /\t/,$hash5{$de};
       push @arry3,[@arry2];
}
}  
 for my $j (0..$#name){
 my $sum=0;
 for my $i(0..$#arry3){
  $sum +=$arry3[$i][$j];
}
print OUT "\t$sum";
}
print OUT "\n";
}

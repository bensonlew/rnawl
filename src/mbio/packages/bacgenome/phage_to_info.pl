#! /user/bin/perl
use strict;
use Bio::SeqIO;
#author:gaohao
########
##主要是将gff文件转化成phage_info.txt文件
#############

if(@ARGV!=4) {
    print STDERR "perl $0 type file scaf output\n";
    exit;
}
my ($type,$input,$scaf,$output)=@ARGV;
my %len=&scaff_len($scaf);
my %hash_chr=('Chromosome'=>'Chr','Chromosome1'=>'Chr1',,'Chromosome2'=>'Chr2','Chromosome3'=>'Chr3','Chromosome4'=>'Chr4');
my %hash_pla=('Plasmid'=>'p','PlasmidA'=>'pA',,'PlasmidB'=>'pB','PlasmidC'=>'pC','PlasmidD'=>'pD','PlasmidE'=>'pE',,'PlasmidF'=>'pF','PlasmidG'=>'pG','PlasmidH'=>'pH');
my %hash3;
foreach (sort keys %len){
 if(exists $hash_chr{$_}){
   $hash3{$hash_chr{$_}} =$_ . "\t" . $len{$_};
}elsif( exists $hash_pla{$_}){
   $hash3{$hash_pla{$_}} =$_ . "\t" . $len{$_};
}
}
open (IN,$input) || die $!;
open (OUT,">$output") || die $!;
while(<IN>){
  chomp;
    next if(/^Gene id/);
  my @temp=split /\t/;
  my @seq=split /_/,$temp[1];
if($type == "complete"){
if(exists $hash3{$seq[0]}){
my @temp2 =split /\t/,$hash3{$seq[0]};
    print OUT "$temp2[0]\t$temp2[1]\t$temp[0]\t$temp[2]\t$temp[3]\n";
 }    
}
 if($type == "uncomplete"){
  if(exists $len{$seq[0]}){
    print OUT "$seq[0]\t$len{$seq[0]}\t$temp[0]\t$temp[2]\t$temp[3]\n";
}
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

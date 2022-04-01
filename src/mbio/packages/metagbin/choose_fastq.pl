#! /usr/bin/perl

my ($fa,$inlist,$out)=@ARGV;

open (IN,$inlist) || die $!;
my %hash;
while(<IN>){
   chomp;
   my @temp=split /\t/;
   $hash{$temp[0]}=1;
}
close IN;
open (IN2,$fa)	|| die $!;
open (OUT,">$out")  || die $!;
while(<IN2>){
   chomp;
   if(/^@(.*)/){
    if(exists $hash{$1}){
      my $id =$_;
      my $d =<IN2>;
      my $dc =<IN2>;
      my $ss = <IN2>;
      print OUT "$id\n$d$dc$ss";
  }
  }
}
close IN2;
close OUT;

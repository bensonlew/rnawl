#! /usr/bin/perl -W

my ($input,$output)=@ARGV;
open (IN,$input) || die $!;
open (OUT,">$output") || die $!;
my $total=0;
my $de = <IN>;
print OUT "$de";
while(<IN>){
  chomp;
  my @temp =split /\t/;
  if ($temp[0]<255){
   print OUT "$temp[0]\t$temp[1]\n";
  }else{
   $total +=$temp[1];
  }
}
my $des = sprintf("%.3f",$total);
print OUT "255\t$des\n";

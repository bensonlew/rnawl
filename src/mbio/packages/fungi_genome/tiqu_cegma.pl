#! /usr/bin/perl -W

my ($input,$output) =@ARGV;

open (IN,$input) || die $!;
open (OUT,">$output") || die $!;
print OUT "Complete(%)\tPartial(%)Total\tAverage\n";
my ($com,$par,$total,$average);
while(<IN>){
 chomp;
 if(/\s+Complete\s+([0-9]*)\s+([0-9.]*)\s+-\s+([0-9]*)\s+([0-9.]*)\s+.*/){
  $com=$2;
  $total=$3;
  $average=$4;
}elsif(/\s+Partial\s+([0-9]*)\s+([0-9.]*)\s+-.*/){
  $par=$2;
}
}
print OUT "$com\t$par\t$total\t$average\n";
close IN;
close OUT;

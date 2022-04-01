#! perl -w

my ($input,$output)=@ARGV;
open(IN,"$input") || die $!;
open(OUT,">$output") || die $!;
print OUT "Total Reads No.\tTotal Bases (bp)\tLargest (bp)\tAverage Len (bp)\n";

while(<IN>){
  chomp;
if($_=~/Number of sequences:\s+([0-9]+)/){
print OUT "$1\t";
}elsif($_=~/Total # residues:\s+([0-9]+)/){
print OUT "$1\t";
}elsif($_=~/Largest:\s+([0-9]+)/){
print OUT "$1\t";
}elsif($_=~/Average length:\s+([0-9.]+)/){
print OUT "$1\n";
}else{
next;
}
}
close IN;

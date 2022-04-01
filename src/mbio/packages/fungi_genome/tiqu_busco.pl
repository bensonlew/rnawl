#! /usr/bin/perl -W

my ($input,$output)=@ARGV;

open (IN,$input) || die $!;
open (OUT,">$output") || die $!;
print OUT "Comeplete(%)\tComplete Duplicated(%)\tFragmented(%)\tMissing(%)\n";
my ($com,$dup,$fra,$miss);
while(<IN>){
   chomp;
  if(/\s+C\:([0-9.]*)\%\[S\:([0-9.]*)%,D\:([0-9.]*)\%\],F\:([0-9.]*)\%,M\:-*([0-9.]*)\%,.*/){
     #print "$1\n$3\n$4\n$5\n";
     $com = $1;
     $dup = $3;
     $fra = $4;
     $miss = $5;
}
}
print OUT "$com\t$dup\t$fra\t$miss\n";
close IN;
close OUT;

#!perl -w

my ($input,$domain_type,$output)=@ARGV;
open (IN,$input) || die $!;
open (OUT,">$output")	|| die $!;

while(<IN>){
  chomp;
  if(/^Taxon/){
  print OUT "$_\n"; 
}
  my $de=d__ . $domain_type;
  if($_ =~/$de/){
  print OUT "$_\n";
}
}
close IN;
close OUT;

#! /usr/bin/perl -W
if (@ARGV !=3){
   print "perl $0 align.table anno_file anno\n";
   exit;
}
my ($align,$anno_file,$anno)=@ARGV;

open (IN,$align) || die $!;
open (IN2,$anno_file) || die $!;
open (OUT,">$anno") || die $!;

print OUT "#Query\tQS_id\tClass\tIdentity(%)\tAlign_len\n";
my %hash;
while(<IN>){
  chomp;
  next if(/^Score/);
  my @temp =split /\t/;
  my ($gene_id)=$temp[5]=~/(.*)_1/;  
  my $des = $temp[10] . "\t" . $temp[3] . "\t" . $temp[2];
  $hash{$gene_id} = $des;
}
close IN;
my %anno;
while(<IN2>){
 chomp;
 my @temp =split /\t/;
 $anno{$temp[0]}=$temp[1];
}
close IN2;

foreach (sort keys %hash){
  my @temp =split /\t/,$hash{$_};
 if(exists $anno{$temp[0]}){
  print OUT  "$_\t$temp[0]\t$anno{$temp[0]}\t$temp[1]\t$temp[2]\n";
}
}
close OUT;

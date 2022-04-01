#! /usr/bin/perl

if(@ARGV !=4){
   print "perl $0 check.summary.xls bin_stats_ext.tsv checkm.taxon.xls prefix\n";
   exit;
}

my ($summary,$bin_stat,$taxon_file,$pre) =@ARGV;
open (IN,$summary) || die $!;
open (IN2,$bin_stat) || die $!;
open (IN3,$taxon_file) || die $!;
open (OUT,">$pre.marker.xls") || die $!;
open (OUT2,">$pre.bin.summary.xls") || die $!;
print OUT "Bin Id\tmarkers\tmarker sets\t0\t1\t2\t3\t4\t5+\n";
print OUT2 "Bin Id\tGenome size\tScaffolds Num\tLongest scaffold\tN50 (scaffolds)\tMean scaffold length\tCompleteness\tContamination\tStrain heterogeneity\tDomain\n";
my (%hash,%taxon);
while(<IN3>){
  chomp;
  my @temp =split /\t/;
  $taxon{$temp[0]}=$temp[1];
}
close IN3;
while (<IN2>){
  chomp;
  my @temp =split /\t/;
  my ($scaffolds)=$temp[1]=~/'# scaffolds': ([0-9]*),/;
  my ($gsize)=$temp[1]=~/'Genome size': ([0-9]*),/;
  my ($lscf)=$temp[1]=~/'Longest scaffold': ([0-9]*),/;
  my ($n_scf)=$temp[1]=~/'N50 \(scaffolds\)': ([0-9]*),/;
  my ($mean_scaf)=$temp[1]=~/'\Mean scaffold length': ([0-9.]*)}/; 
  my $des=$gsize . "\t" . $scaffolds . "\t" . $lscf . "\t" . $n_scf . "\t" . $mean_scaf;
  $hash{$temp[0]}=$des;
}
close IN2;
my %hash2;
while(<IN>){
  chomp;
  next if(/^Bin Id/);
  my @temp =split /\t/;
  my @temp2=split /\s/,$temp[1];
  $hash2{$temp[0]}=$temp[-3] . "\t" . $temp[-2] . "\t" . $temp[-1];
  print OUT "$temp[0]\t$temp[3]\t$temp[4]\t$temp[5]\t$temp[6]\t$temp[7]\t$temp[8]\t$temp[9]\t$temp[10]\n";
}
close IN;

foreach (sort {$hash2{$b} <=> $hash2{$a}}keys %hash2){
  if (exists $hash{$_}){
      print OUT2 "$_\t$hash{$_}\t$hash2{$_}\t";
   }
   if(exists $taxon{$_}){
      print OUT2 "$taxon{$_}\n";
   }else{
           print OUT2 "Bacteria\n";
   }
  }


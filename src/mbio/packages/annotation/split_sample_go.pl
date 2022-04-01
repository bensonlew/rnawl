#! /usr/bin/perl -W

if (@ARGV !=4){
  print "perl $0 GO.level2.info.xls  GO.level2.info.xls GO.level2.info.xls output\n";
  exit;
}
my ($file,$file2,$file3,$output)=@ARGV;

&fenli_sample($file,level2);
&fenli_sample($file2,level3);
&fenli_sample($file3,level4);

sub fenli_sample{
  my ($level_file,$prefix)=@_;
 open (IN,$level_file) || die $!;
my $name=<IN>;
my @name=split /\t/,$name;
my @sample=@name[2..$#name-1];
my @abun;
while(<IN>){
  chomp;
 push @abun,$_;
}
close IN;
for my $de (0..$#sample){ 
  my (%hash,$total);
for my $s (@abun){
  my @temp =split /\t/,$s; 
  my $des=$temp[0] . "\t" . $temp[1];
  $hash{$des}=$temp[$de+2];
  $total +=$temp[$de+2];
}
my $file=$sample[$de] . "_" . $prefix;
`mkdir -p $output`;
open (OUT,">$file.info.xls") || die $!;
 if($prefix eq "level2"){
    print OUT "Term type\tGo Term(level2)\tAbundance\tpercent\n";
}elsif($prefix eq "level3"){
   print OUT "Term type\tGo Term(level3)\tAbundance\tpercent\n";
}elsif($prefix eq "level4"){
   print OUT "Term type\tGo Term(level4)\tAbundance\tpercent\n";
}
foreach (sort keys %hash){
  my $de=$hash{$_}/$total;
  my $num=sprintf("%.3f",$de);
  print OUT "$_\t$hash{$_}\t$num\n";
}
close OUT;
}
}

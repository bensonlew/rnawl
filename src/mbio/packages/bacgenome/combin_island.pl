#! perl -W

my ($type,$dimob,$islander)=@ARGV;
open (IN,$dimob) || die $!;
open (IN2,$islander)	|| die $!;
my $out ="all.xls";
open (OUT,">$out")	|| die $!;
my %hash;
while(<IN>){
  chomp;
    next if(/^location/);
  my @temp=split /\t/;
  my $de = $temp[3] . "\t" . $temp[4];
  $hash{$de}=$temp[0] . "\t" . $temp[1];
}
close IN;
my %hash2;
while(<IN2>){ 
  chomp;
    next if(/^#/);
  my @temp=split /\t/;
  my $de = $temp[3] . "\t" . $temp[4];
  $hash2{$de}=$temp[0] . "\t" . $temp[1];
}
close IN2;
foreach (sort keys %hash){
 if(exists $hash2{$_}){
   my @temp =split /\t/,$hash{$_};
   my @temp2 =split /\t/,$hash2{$_};
   if($temp[0] eq $temp2[0]){
    my @arry=split /\t/,$_;
   print OUT "$temp[0]\tGI\tBoth\t$arry[0]\t$arry[1]\n";
}
}else{
  my @temp3 =split /\t/,$hash{$_};
  my @arry2 =split /\t/,$_;
  print OUT "$temp3[0]\tGI\t$temp3[1]\t$arry2[0]\t$arry2[1]\n";
}
}
foreach (sort keys %hash2){
 if(!(exists $hash{$_})){
    my @temp3 =split /\t/,$hash2{$_};
   my @arry2 =split /\t/,$_;
    print OUT "$temp3[0]\tGI\t$temp3[1]\t$arry2[0]\t$arry2[1]\n";
}
}
close OUT;
my $output="all.stat.xls";
&get_result($out,$output);
sub get_result{
    my ($input,$output)=@_;
    open (IN,$input) || die $!;
my %hash;    
while(<IN>){
      chomp;
      my @temp =split /\t/;
       if($type eq 'uncomplete'){
      my ($name,$de)=$temp[0]=~/([a-zA-Z]*)([0-9]*)/;
      $sde =$name . "\t" . $_ . ";";
      $hash{$de} .=$sde;
}elsif($type eq 'complete'){
   $hash{$temp[0]} .=$_ . ";";

}      
}
close IN;
open (OUT3,">all.island.xls") || die $!;
print OUT3 "location\tGI No.\tMethod\tIsland Start\tIsland End\n";
if($type eq 'uncomplete'){
foreach(sort {$a <=> $b} keys %hash){
   my ($ss)=$hash{$_} =~/(.*);$/;
   my @arry =split /;/,$ss;
   for my $dee (@arry){
   my @temp=split /\t/,$dee;
   my $des=join("\t",@temp[2..$#temp]);
   my $de=$temp[0] . $_;
   print OUT3 "$de\t$des\n";
}
}
}elsif($type eq 'complete'){
  foreach(sort {$a cmp $b} keys %hash){
   my ($ss)=$hash{$_} =~/(.*);$/;
   my @arry =split /;/,$ss;
   for my $dee (@arry){
   my @temp=split /\t/,$dee;
   my $des=join("\t",@temp[1..$#temp]);
   my $de=$temp[0];
   print OUT3 "$de\t$des\n";
}
}


}
close OUT3;
}
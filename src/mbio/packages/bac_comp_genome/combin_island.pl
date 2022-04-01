#! perl -W

my ($dimob,$islander)=@ARGV;
open (IN,$dimob) || die $!;
open (IN2,$islander)	|| die $!;
open (OUT,">all.island.xls")	|| die $!;
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

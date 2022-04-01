#! perl -W
use Bio::SeqIO;
my ($type,$input,$anno,$output) =@ARGV;

open (IN2,$anno) || die $!;
open (OUT2,">$output.GI_detail.xls") || die $!;
open (OUT,">$output.GI_summary.xls") || die $!;
print OUT "Location\tGI No.\tIsland Start\tIsland End\tLength(bp)\tMethod\tCDS No.\n";
print OUT2 "Location\tGI No.\tGene Id\tGene Start\tGene End\tStrand\tGene Name\tNR Description\tCOG Type\n";
my %anno;
while(<IN2>){
   chomp;
  next if(/^Gene/); 
  my @temp = split /\t/;
  my $de=$temp[2] . "\t" . $temp[3];
  if($temp[1] eq "-"){
  $anno{$de}=$temp[5] . "\t" . $temp[0] . "\t" . $temp[3] . "\t" . $temp[2] . "\t" . $temp[10] . "\t" . $temp[14] . "\t" . $temp[6] . "\t" . $temp[1];
}elsif($temp[1] eq "+"){
     $anno{$de}=$temp[5] . "\t" . $temp[0] . "\t" . $temp[2] . "\t" . $temp[3] . "\t" . $temp[10] . "\t" . $temp[14] . "\t" . $temp[6] . "\t" . $temp[1];
}    
}
close IN2;
open (IN,$input) || die $!;
my %hash_chr=('Chromosome'=>'Chr','Chromosome1'=>'Chr1',,'Chromosome2'=>'Chr2','Chromosome3'=>'Chr3','Chromosome4'=>'Chr4');
my %hash_pla=('Plasmid'=>'p','PlasmidA'=>'pA',,'PlasmidB'=>'pB','PlasmidC'=>'pC','PlasmidD'=>'pD','PlasmidE'=>'pE',,'PlasmidF'=>'pF','PlasmidG'=>'pG','PlasmidH'=>'pH');
my $num=1;
while(<IN>){
   chomp;
   next if(/^location/);
   my @temp =split /\t/;
   my @gene;
   my $de;
    my $ddd;
   if($num <=9){
          $de="GI0" . $num;
        }else{
          $de="GI" . $num;
             }
  if($type eq 'complete'){
   foreach (sort keys %anno){
     my @temp1=split /\t/,$_;
     my @arry=split /\t/,$anno{$_};
     #if(2>1){
     if($arry[0] eq $hash_chr{$temp[0]} || $arry[0] eq $hash_pla{$temp[0]} || $arry[0] eq $temp[0]){
           $ddd = $arry[0];
     if($temp1[0] >= $temp[3] && $temp1[1] <= $temp[4]){
          push @gene,$arry[1];
         print OUT2 "$arry[0]\t$de\t$arry[1]\t$arry[2]\t$arry[3]\t$arry[7]\t$arry[5]\t$arry[6]\t$arry[4]\n";
  }
}
}
}elsif($type eq 'uncomplete'){
foreach (sort keys %anno){
     my @temp1=split /\t/,$_;
     my @arry=split /\t/,$anno{$_};
 if($arry[0] eq $temp[0]){
      $ddd = $arry[0];
  if($temp1[0] >= $temp[3] && $temp1[1] <= $temp[4]){
          push @gene,$arry[1];
         print OUT2 "$temp[0]\t$de\t$arry[1]\t$arry[2]\t$arry[3]\t$arry[7]\t$arry[5]\t$arry[6]\t$arry[4]\n";
}
}
}
}
   my $count=@gene;
print "$temp[4]\t$temp[3]\n";
my $len =$temp[4] - $temp[3];
   $num++;
    print  OUT "$ddd\t$de\t$temp[3]\t$temp[4]\t$len\t$temp[2]\t$count\n";
}
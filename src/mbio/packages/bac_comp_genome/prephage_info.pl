#! perl -W
use Bio::SeqIO;
my ($type,$input,$gff,$output,$db) =@ARGV;

open (IN2,$gff) || die $!;
open (OUT,">$output.summary.xls") || die $!;
print OUT "Location\tPh No.\tPh Start\tPh End\tStrand\tLength(bp)\tCDS No.\tPossible Phage\tGC(%)\tGene_list\n";

my %id_phage;
open(INF,$db) || die $!;
while(<INF>){
    chomp;
    next if(/^Accession/);
    my @tmp = split /\t/;
    my $accession = $tmp[0];
    my $phage_name = $tmp[1];
    $id_phage{$accession} = $phage_name;
}
close INF;

my %anno;
while(<IN2>){
   chomp;
  next if(/^Gene/);
  my @temp = split /\t/;
  my @temp2 = split /\_ORF/,$temp[1];
  if($temp[4] eq "-"){
   my $de=$temp[3] . "\t" . $temp[2];
  $anno{$de}{$temp2[0]}=$temp2[0] . "\t" . $temp[0];
}elsif($temp[4] eq "+"){
  my $de=$temp[2] . "\t" . $temp[3];
  $anno{$de}{$temp2[0]}=$temp2[0] . "\t" . $temp[0] ;
}
}
close IN2;
open (IN,$input) || die $!;
my $num_pre;
my $total_len;
my $num=1;
while(<IN>){
   chomp;
   next if(/^#/);
   my @temp =split /\t/;
   my @gene;
   my $de;
   my $ddd;
   if($num <=9){
          $de="Ph0" . $num;
        }else{
          $de="Ph" . $num;
        }
   my %hash;
if($type eq 'complete'){
   foreach my $key (sort keys %anno){
       my @temp1=split /\t/,$key;
     foreach my $key2 (sort keys %{$anno{$key}}){
        my @arry=split /\t/,$anno{$_};
       if(2>1){
           $ddd=$arry[0];
           if($temp1[0] >= $temp[3] && $temp1[1] <= $temp[4]){
               if (exists $hash{$de}){
                    $hash{$de} .= ";" . $arry[1];
               }else{
                    $hash{$de} =$arry[1];
               }
           }
       }
     }
   }
}elsif($type eq 'uncomplete'){
foreach my $key (sort keys %anno){
     my @temp1=split /\t/,$key;
    foreach my $key2 (sort keys %{$anno{$key}}){
     my @arry=split /\t/,$anno{$key}{$key2};
 if($arry[0] eq $temp[0]){
         print($arry[0]);
         print($arry[1]);
          $ddd=$arry[0];
          if($temp1[0] >= $temp[3] && $temp1[1] <= $temp[4]){
         if (exists $hash{$de}){
                    $hash{$de} .= ";" . $arry[1];
               }else{
                    $hash{$de} =$arry[1];
               }
}
  }
}
}
}
   my $target_accession = $temp[12];
   my $target_name = $id_phage{$target_accession};
   my $final_possible_phage = $target_name."_".$target_accession;
   print  OUT "$ddd\t$de\t$temp[3]\t$temp[4]\t$temp[21]\t$temp[5]\t$temp[-2]\t".$final_possible_phage."\t$temp[11]\t$hash{$de}\n";
   $num ++;
}
close IN;

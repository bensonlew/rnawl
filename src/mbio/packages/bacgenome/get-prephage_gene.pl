#! perl -W
use Bio::SeqIO;
my ($type,$input,$anno,$output,$db) =@ARGV;

open (IN2,$anno) || die $!;
open (OUT2,">$output.detail.xls") || die $!;
open (OUT,">$output.summary.xls") || die $!;
open (OUT3,">$output.stat.xls") || die $!;
print OUT "Location\tPh No.\tPh Start\tPh End\tLength(bp)\tCDS No.\tPossible Phage\tGC(%)\n";
print OUT2 "Location\tPh No.\tGene Id\tStrand\tGene Start\tGene End\tGene Name\tNR Description\tCOG Type\n";

# modify by ysh in 20190409 for add phage detail name
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
#

my %anno;
while(<IN2>){
   chomp;
  next if(/^Gene/); 
  my @temp = split /\t/;
  my $de=$temp[2] . "\t" . $temp[3];
  if($temp[1] eq "-"){
  $anno{$de}=$temp[5] . "\t" . $temp[0] . "\t" . $temp[3] . "\t" . $temp[2] . "\t" . $temp[10] . "\t" . $temp[14] .  "\t" . $temp[6] . "\t" . $temp[1];
}elsif($temp[1] eq "+"){
$anno{$de}=$temp[5] . "\t" . $temp[0] . "\t" . $temp[2] . "\t" . $temp[3] . "\t" . $temp[10] . "\t" . $temp[14] . "\t" . $temp[6] . "\t" . $temp[1];
}    
}
close IN2;
open (IN,$input) || die $!;
my %hash_chr=('Chromosome'=>'Chr','Chromosome1'=>'Chr1',,'Chromosome2'=>'Chr2','Chromosome3'=>'Chr3','Chromosome4'=>'Chr4');
my %hash_pla=('Plasmid'=>'p','PlasmidA'=>'pA',,'PlasmidB'=>'pB','PlasmidC'=>'pC','PlasmidD'=>'pD','PlasmidE'=>'pE',,'PlasmidF'=>'pF','PlasmidG'=>'pG','PlasmidH'=>'pH');
my $num_pre;
my $total_len;
my $num=1;
my $gene_num=0;
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
if($type eq 'complete'){
   foreach (sort keys %anno){
     my @temp1=split /\t/,$_;
     my @arry=split /\t/,$anno{$_};
     #if($arry[0] eq $hash_chr{$temp[0]} || $arry[0] eq $hash_pla{$temp[0]}){
     if(2>1){
             $ddd=$arry[0];
          if($temp1[0] >= $temp[3] && $temp1[1] <= $temp[4]){
          if (exists $hash_chr{$temp[0]}){
            if ($hash_chr{$temp[0]} eq $arry[0]){
                print OUT2 "$arry[0]\t$de\t$arry[1]\t$arry[7]\t$arry[2]\t$arry[3]\t$arry[5]\t$arry[6]\t$arry[4]\n";
				$gene_num +=1;
            }
          }
          elsif (exists $hash_pla{$temp[0]}){
            if ($hash_pla{$temp[0]} eq $arry[0]){
                print OUT2 "$arry[0]\t$de\t$arry[1]\t$arry[7]\t$arry[2]\t$arry[3]\t$arry[5]\t$arry[6]\t$arry[4]\n";
				$gene_num +=1;
            }
          }
}
  }
}
}elsif($type eq 'uncomplete'){
foreach (sort keys %anno){
     my @temp1=split /\t/,$_;
     my @arry=split /\t/,$anno{$_};
 if($arry[0] eq $temp[0]){
          $ddd=$arry[0];
          if($temp1[0] >= $temp[3] && $temp1[1] <= $temp[4]){
         print OUT2 "$temp[0]\t$de\t$arry[1]\t$arry[7]\t$arry[2]\t$arry[3]\t$arry[5]\t$arry[6]\t$arry[4]\n";
		 $gene_num +=1;
}
  }
}
}
   my $target_accession = $temp[12];
   my $target_name = $id_phage{$target_accession};
   my $final_possible_phage = $target_name."_".$target_accession;
   print  OUT "$ddd\t$de\t$temp[3]\t$temp[4]\t$temp[5]\t$gene_num\t".$final_possible_phage."\t$temp[11]\n";
   $num_pre++;
   $total_num +=$temp[5];
   $num++;
}
close IN;
print OUT3 "Ph No.\t Total Length\n";
print OUT3 "$num_pre\t$total_num\n";
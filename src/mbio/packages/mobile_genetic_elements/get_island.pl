#! perl -W
use Bio::SeqIO;
my ($input,$gff,$output) =@ARGV;

open (IN2,$gff) || die $!;
open (OUT,">$output.GI_summary.xls") || die $!;
print OUT "Location\tGI No.\tIsland Start\tIsland End\tLength(bp)\tMethod\tCDS No.\tGene_list\n";
my %anno;
while(<IN2>){
   chomp;
  next if(/^#/);
  my @temp = split /\t/;
  my @temp1 = split /_ORF/,$temp[1];
  my $de;
 if($temp[6] eq "-"){
      $de=$temp[4] . "\t" . $temp[3];
      $anno{$de}=$temp[1] . "\t" . $temp[0] . "\t" . $temp[3] . "\t" . $temp[4] . "\t" . $temp[6];
 }elsif($temp[6] eq "+"){
      $de=$temp[3] . "\t" . $temp[4];
      $anno{$de}=$temp[1] . "\t" . $temp[0] . "\t" . $temp[3] . "\t" . $temp[4] . "\t" . $temp[6];
 }
}
close IN2;
open (IN,$input) || die $!;
my $num=1;
while(<IN>){
   chomp;
   next if(/^Sample/);
   my @temp =split /\t/;
   my @gene;
   my $de;
   if($num <=9){
       $de="GI0" . $num;
   }else{
       $de="GI" . $num;
   } 
   foreach (sort keys %anno){
       my @temp1=split /\t/,$_;
       my @arry=split /\t/,$anno{$_};
       if($temp[1] eq $arry[1]){
           if($temp1[0] >= $temp[2] && $temp1[1] <= $temp[3]){
               push @gene,$arry[0];
  }
}
}
my $count=@gene;
my $gene_list = join(",",@gene);
my $len;
my $start;
my $end;
if ($temp[3] > $temp[2]){
    $start = $temp[2];
    $end = $temp[3];
    $len =$temp[3] - $temp[2];
  }
elsif ($temp[3] <$temp[2]){
    $len =$temp[2] - $temp[3];
    $start = $temp[3];
    $end = $temp[2];
}
   $num++;
    print  OUT "$temp[1]\t$de\t$start\t$end\t$len\t$temp[4]\t$count\t$gene_list\n";
}

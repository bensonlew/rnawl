#! perl -W
use Bio::SeqIO;
my ($input,$gff,$output) =@ARGV;

open (IN2,$gff) || die $!;
open (OUT,">$output.GI_summary.xls") || die $!;
print OUT "Location\tGI No.\tIsland Start\tIsland End\tLength(bp)\tMethod\tCDS No.\tGene_list\n";
my %anno;
while(<IN2>){
   chomp;
  next if(/^Gene/);
  my @temp = split /\t/;
  my @temp1 = split /_ORF/,$temp[1];
  my $de;
 if($temp[4] eq "-"){
      $de=$temp[3] . "\t" . $temp[2];
      $anno{$de}=$temp1[0] . "\t" . $temp[0] . "\t" . $temp[3] . "\t" . $temp[2] . "\t" . $temp[4];
 }elsif($temp[4] eq "+"){
      $de=$temp[2] . "\t" . $temp[3];
      $anno{$de}=$temp1[0] . "\t" . $temp[0] . "\t" . $temp[2] . "\t" . $temp[3] . "\t" . $temp[4];
 }
}
close IN2;
open (IN,$input) || die $!;
my $num=1;
while(<IN>){
   chomp;
   next if(/^location/);
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
       if($temp[0] eq $arry[0]){
           if($temp1[0] >= $temp[3] && $temp1[1] <= $temp[4]){
               push @gene,$arry[1];
  }
}
}
my $count=@gene;
my $gene_list = join(",",@gene);
my $len =$temp[4] - $temp[3];
   $num++;
    print  OUT "$temp[0]\t$de\t$temp[3]\t$temp[4]\t$len\t$temp[2]\t$count\t$gene_list\n";
}

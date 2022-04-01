#! perl -W
use Bio::SeqIO;
my ($input,$step_clean,$min_len,$type,$prefix)=@ARGV;

#open (IN,$input) || die $!;
open (OUT2,">$prefix.clean.len.xls") || die $!;
print OUT2 "Len\tnum\ttotal_len\n";
if ($type eq "nanopore"){
  open (OUT,">$prefix.Nanopore_statistics.xls") || die $!;
}elsif($type eq "pacbio"){
  open (OUT,">$prefix.PacBio_statistics.xls") || die $!;
}

print  OUT "Total Reads No.\tTotal Bases (bp)\tLargest (bp)\tAverage Len (bp)\n";

my @len_arry2;
my $total_len;
my $num =0;
my $in = Bio::SeqIO->new(-file => "$input", -format => "fastq");
while(my $seq = $in->next_seq())  {
    #print $seq->id, "\t", $seq->length, "\n";
    my $len = $seq->length;
    if($len >=1000){
    #print $seq->id, "\t", $seq->length, "\n"; 
    $total_len += $len;
    push @len_arry2,$len;
      $num ++;
}
}
my @len=sort {$a <=> $b} @len_arry2;
my $largest = $len[-1];
my $average = $total_len/$num;
my $des = sprintf("%.2f",$average);
print OUT "$num\t$total_len\t$largest\t$des\n";

my %hash3=&pacbio_num($step_clean,$min_len,@len_arry2);
my %hash4=&pacbio_len($step_clean,$min_len,@len_arry2);

foreach my $x (sort {$a <=> $b} keys %hash3) {
  if(exists $hash4{$x}){
     my $des=$hash4{$x}/1000000;
        print OUT2 "$x\t$hash3{$x}\t$des\n";
}
}

sub pacbio_len{ 
  my ($step,$min_len,@arry)=@_;
  @arry = sort {$a <=> $b} @arry;
  my %datalist;
  my $range=$step;
  my $r=$min_len;
  my $totallen=0;
  foreach my $data(@arry){
   $totallen +=$data;
   }
  for (my $i=0;$i<int(($arry[-1]-$min_len)/$step);$i++){
  my $num=0;
   foreach my $data (@arry) {
     if($data <=$r+$range*$i) {
	$num +=$data;
     }
   }
   $datalist{$r+$range*$i}=$totallen-$num;
   }
  return %datalist; 
   }

sub pacbio_num{ 
  my ($step,$min_len,@arry)=@_;
  @arry = sort {$a <=> $b} @arry;
  my %datalist;
  my $range=$step;
  my $r=$min_len;
  
for (my $i=0;$i<int(($arry[-1]-$min_len)/$step);$i++){
  my $num=0;
  foreach my $data (@arry) {
     if($data <=$r+$range*$i) {
	$num ++;
     }
  }
  $datalist{$r+$range*$i}=$num;
  } 
   return %datalist; 
   }


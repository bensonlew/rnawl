#perl -W
use strict;
use Bio::SeqIO;

my ($input,$input2,$prefix)=@ARGV;

open (IN,$input) || die $!;
open (IN2,$input2) || die $!;
open (OUT,">$prefix.scaffolds.len.xls") || die $!;
open (OUT2,">$prefix.contigs.len.xls") || die $!;
print OUT "Len\tnum\ttotal_len\n";
print OUT2 "Len\tnum\ttotal_len\n";
my @len_arry = &lenth($input);
my @len_arry2 = &lenth($input2);
@len_arry =sort {$a <=> $b} @len_arry;
@len_arry2 =sort {$a <=> $b} @len_arry2;
my $step_raw = &get_step(@len_arry);
print "$step_raw\n";
my %hash1=&pacbio_num($step_raw,@len_arry);
my %hash2=&pacbio_len($step_raw,@len_arry);

foreach my $x (sort {$b cmp $a} keys %hash1) {
  if(exists $hash2{$x}){
     my $des=$hash2{$x};
        print OUT "$x\t$hash1{$x}\t$des\n";
}
}
my $step_clean = &get_step(@len_arry2);
print "$step_clean\n";
my %hash3=&pacbio_num($step_clean,@len_arry2);
my %hash4=&pacbio_len($step_clean,@len_arry2);

foreach my $x (sort {$b cmp $a} keys %hash3) {
  if(exists $hash4{$x}){
     my $des=$hash4{$x};
        print OUT2 "$x\t$hash3{$x}\t$des\n";
}
}
sub lenth{
    my ($input)=@_;
    my @len;
    my $in = Bio::SeqIO->new(-file => "$input", -format => "fasta");
    while(my $seq = $in->next_seq())  {
       my $len=$seq->length;
       push @len,$len;
    }
    return @len;
}
sub get_step{
  my (@arry,$step)=@_;
   if($arry[-1] =~/^[0-9]{5}$/){
       if (int($arry[-1]/10000) < 3){
           $step =1000;
       }elsif(int($arry[-1]/10000) >= 3){
           $step =int($arry[-1]/30000)*1000;
       }
   }elsif($arry[-1] =~/^[0-9]{6}$/){
       if (int($arry[-1]/100000) < 3){
           $step =10000;
       }elsif(int($arry[-1]/100000) >= 3){
           $step =int($arry[-1]/300000)*10000;
       }
   }
   elsif($arry[-1] =~/^[0-9]{7}$/){
       if (int($arry[-1]/1000000) < 3){
           $step =100000;
       }elsif(int($arry[-1]/1000000) >= 3){
           $step =int($arry[-1]/3000000)*100000;
       }
   }
   return $step
}
sub pacbio_len{ 
  my ($step,@arry)=@_;
  @arry = sort {$a <=> $b} @arry;
  my %datalist;
  my $num=0;
  my $range=$step;
  my $r=$range;
  my $totallen=0;
  foreach my $data(@arry){
   $totallen +=$data;
 }
 my $num2=0;
  foreach my $data (@arry) {
      if($data >= $r*9) {
	      $num2 +=$data;
	      my $de = ">" . $r*9/1000000;
	      my $des = sprintf("%.2f",$num2/$totallen*100);
    	  $datalist{$de}=$des . "\t" . 9;
        }
    }
  for (my $i=0;$i<=8;$i++) {
     my $num=0;
     my $des;
     foreach my $data (@arry) {
        if($data >= $r*(8-$i)) {
	        $num +=$data;
	        if ($i == 8 ){
	            $des =$r*(8-$i) . "-" . $r*(8-$i+1)/1000000;
	        }else{
	            $des =$r*(8-$i)/1000000 . "-" . $r*(8-$i+1)/1000000;
	        }
        }
     }
     my $de = sprintf("%.2f",$num/$totallen*100);
     $datalist{$des}=$de . "\t" . (8-$i);
  }
  return %datalist; 
   }

sub pacbio_num{ 
  my ($step,@arry)=@_;
  @arry = sort {$a <=> $b} @arry;
  my %datalist;
  my $range=$step;
  my $r=$range;
  my $num2=0;
  foreach my $data (@arry) {
      if($data >= $r*9) {
	      $num2 ++;
	      my $de = ">" . $r*9/1000000;
    	  $datalist{$de}=$num2;
        }
    }
  for (my $i=0;$i<=8;$i++) {
     my $num=0;
     my $des;
     foreach my $data (@arry) {
        if($r*(8-$i) < $data && $data<= $r*(8-$i+1)) {
	        $num ++;
	        }
	    if ($i == 8 ){
	        $des =$r*(8-$i) . "-" . $r*(8-$i+1)/1000000;
	    }else{
	        $des =$r*(8-$i)/1000000 . "-" . $r*(8-$i+1)/1000000;
	    }
     }
      $datalist{$des}=$num ;
  }
   return %datalist; 
   }


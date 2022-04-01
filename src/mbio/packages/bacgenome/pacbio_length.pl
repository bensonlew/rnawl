#perl -W

my ($input,$input2,$step_raw,$step_clean,$num_raw,$num_clean,$prefix)=@ARGV;

open (IN,$input) || die $!;
open (IN2,$input2) || die $!;
open (OUT,">$prefix.raw.len.xls") || die $!;
open (OUT2,">$prefix.clean.len.xls") || die $!;
print OUT "Len\tnum\ttotal_len\n";
print OUT2 "Len\tnum\ttotal_len\n";
my @len_arry;
<IN>;
while(<IN>){
  chomp;
  my @temp=split /\,/;
  push @len_arry,$temp[$num_raw];
}
close IN;

my @len_arry2;
<IN2>;
while(<IN2>){
  chomp;
  my @temp=split /\,/;
  push @len_arry2,$temp[$num_clean];
}
close IN2;
my %hash1=&pacbio_num($step_raw,@len_arry);
my %hash2=&pacbio_len($step_raw,@len_arry);

foreach my $x (sort {$a <=> $b} keys %hash1) {
  if(exists $hash2{$x}){
     my $des=$hash2{$x}/1000000;
        print OUT "$x\t$hash1{$x}\t$des\n";
}
}
my %hash3=&pacbio_num($step_clean,@len_arry2);
my %hash4=&pacbio_len($step_clean,@len_arry2);

foreach my $x (sort {$a <=> $b} keys %hash3) {
  if(exists $hash4{$x}){
     my $des=$hash4{$x}/1000000;
        print OUT2 "$x\t$hash3{$x}\t$des\n";
}
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
  foreach my $data (@arry) {
     if($data<=$r) {
	$num +=$data;
    	$datalist{$r}=$totallen-$num;
     }
     else {     		
	        $r+=$range;
	        redo;
    }
   }

   my @keys=sort {$a <=> $b} keys %datalist;
   my $i=0;
   my $s;
   my $total=@keys;
   my $n=0;
   for ($i=1;$i<$total;$i++) {
	if ($datalist{$keys[$i]}<=10) {
		$s +=$datalist{$keys[$i]};
	   if ($s>=$n||$i==$total-1||$datalist{$keys[$i+1]}>10) {
		   $datalist{$keys[$i]}= $s;$s=0;
		} else { delete ($datalist{$keys[$i]});}				
	}  else { $n=$datalist{$keys[$i]};}		
     }
  return %datalist; 
   }

sub pacbio_num{ 
  my ($step,@arry)=@_;
  @arry = sort {$a <=> $b} @arry;
  my %datalist;
  my $num=0;
  my $range=$step;
  my $r=$range;

  foreach my $data (@arry) {
     if($data<=$r) {
	$num ++;
    	$datalist{$r}=$num;
     }
     else {     		    	    	   	        	
    		$num=0;
	        $r+=$range;
	        redo;
    }
  }

  my @keys=sort {$a <=> $b} keys %datalist;
  my $i=0;
  my $s;
  my $total=@keys;
  my $n=0;
  for ($i=1;$i<$total;$i++) {
	if ($datalist{$keys[$i]}<=10) {
		$s +=$datalist{$keys[$i]};
		if ($s>=$n||$i==$total-1||$datalist{$keys[$i+1]}>10) {
		   $datalist{$keys[$i]}= $s;$s=0;
		} else { delete ($datalist{$keys[$i]});}				
	}  else { $n=$datalist{$keys[$i]};}		
     } 
   return %datalist; 
   }


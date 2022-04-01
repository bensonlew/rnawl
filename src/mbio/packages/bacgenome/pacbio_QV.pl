#perl -W

my ($input,$step_raw,$step_clean,$num,$len_num,$prefix)=@ARGV;

open (IN,$input) || die $!;
open (OUT,">$prefix.raw.qv.xls") || die $!;
open (OUT2,">$prefix.clean.qv.xls") || die $!;
print OUT "QV\tnum\ttotal_len\n";
print OUT2 "QV\tnum\ttotal_len\n";
my @len_arry;
my @QV_len;
my @clean_len_arry;
my @clean_QV_len;
<IN>;
while(<IN>){
  chomp;
  my @temp=split /\,/;
  my $de=join("\t",$temp[$num],$temp[$len_num]);
  push @len_arry,$temp[$num];
   push @QV_len,$de;
   if($temp[$num] >= 0.8){
     my $des=join("\t",$temp[$num],$temp[$len_num]);
     push @clean_len_arry,$temp[$num];
     push @clean_QV_len,$des;
}
}
close IN;

my %hash1=&pacbio_num($step_raw,@len_arry);
my %hash2=&pacbio_len($step_raw,@QV_len);

foreach my $x (sort {$a <=> $b} keys %hash1) {
  if(exists $hash2{$x}){
     my $des=$hash2{$x};
     my $ss=sprintf("%.3f",$x);
        print OUT "$ss\t$hash1{$x}\t$des\n";
}

}
my %hash3=&pacbio_num($step_clean,@clean_len_arry);
my %hash4=&pacbio_len($step_clean,@clean_QV_len);

foreach my $x (sort {$a <=> $b} keys %hash3) {
  if(exists $hash4{$x}){
     my $des=$hash4{$x};
      my $ss=sprintf("%.3f",$x);
        print OUT2 "$ss\t$hash3{$x}\t$des\n";
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
  foreach my $data (@arry){
    my @temp=split /\t/,$data;
   $totallen +=$temp[1];
 }
  foreach my $data (@arry) {
     my @temp=split /\t/,$data;
     if($temp[0]<=$r) {
	$num +=$temp[1];
    	$datalist{$r}=$totallen-$num;
     }
     else {     		
	        $r +=$range;
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


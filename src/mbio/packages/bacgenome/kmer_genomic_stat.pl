#! /usr/bin/perl
use strict;
use warnings;
use Math::CDF qw(ppois);
#author:gaohao
my $de=@ARGV;
print "$de\n";
if (@ARGV == 4){
my ($file,$logfile,$bases,$output)=@ARGV;
open (OUT,">$output.summary.xls") || die $!;
open (OUT2,">$output.frequency.xls") || die $!;
print OUT2 "k-mer Depth\tFrequency\n";
print OUT "Method\tK-mer\tK-mer individual number\tk-mer Depth\tGenome size(M)\tRevised genome size(M)\tHeterozygous rate\tRepeat\tProbability having the same kmer\tAccuracy rate\tUsed Bases\tAverage Sequencing Depth\n";
open (IN,$file) || die $!;
my (@kn,@ksp,@kn_out,@ksp_out);
my $more_kn =0;
my $sum =0;
my $sum_ksp =0;
my $kmer=17;
my %hash;
while(<IN>){
  chomp;
  my @arr = split(/\s+/,$_);
  my $numb = $arr[0]*$arr[1]; #caculate kmer individuals number 
  $sum += $numb;
  $sum_ksp += $arr[1];
  if ($arr[0]<1001) {
     $hash{$arr[0]}=$numb;
    push (@kn,$numb);
    push (@ksp,$arr[1]);
  }else{
    $more_kn += $numb;
  }
}
close IN;
foreach (sort{$a<=>$b} keys %hash){
  my $de=$hash{$_}/$sum;
  my $des=sprintf("%.3f",$de*100);
  print OUT2 "$_\t$des\n";
}
push (@kn,$more_kn);
my $acc_rate=1-$kn[0]/$sum; #accuracy rate (consider kmer_depth=1 is error number)
my $base=$sum*17;
@kn_out  = &est_freq_Hr(1,@kn); #kn_out[5] is max_num
my $rep = &est_repeat ($kn_out[2],@kn);
&output (@kn_out);
@ksp_out = &est_freq_Hr(0,@ksp); #kmer species results
&output (@ksp_out);

sub est_freq_Hr{
	my ($flag,@in) = @_;
	my ($max_freq,$half_freq,$max_mean,$half_mean,$Hr);
	my $max_num = 0;
    ## compare to find the max freq
	for my $i(1..@in-2){
		if ($in[$i-1] < $in[$i] && $in[$i] > $in[$i+1]){
			if ($max_num < $in[$i]){
				$max_num = $in[$i];
				$max_freq = $i+1;
			}
		}
	}
	
	my $Gsize = $sum/$max_freq;

    ## estimate the average of heterozygous and homozygous peak
	if ($max_freq %2 ==0) {
		$half_freq = $max_freq/2;
		$half_mean = &mean($half_freq,@in);
	}else {
		$half_freq = ($max_freq-1)/2;
		$half_mean = &mean($half_freq,@in);
	}	
	$max_mean = &mean($max_freq,@in); 

    ## Calculate heterozygous rate
	if ($flag == 1){
		$Hr = $half_mean/$kmer/($half_mean + $max_mean); #kn Hr
		return ("kn",$sum,$max_freq,$Gsize,$Hr,$max_num);
	}else {
		$Hr = $half_mean/$kmer/($half_mean + 2*$max_mean); #ksp Hr
		return ("ksp",$sum,$max_freq,$Gsize,$Hr,$max_num);
  }	
}

## calculate mean of three numbers
sub mean{
	my ($mid,@mean_in) = @_;
	my $avg = ($mean_in[$mid-2]+$mean_in[$mid-1]+$mean_in[$mid])/3;
	return $avg;
}

## calculate repeat using area
sub est_repeat{
	my ($peak,@repeat_in) = @_;
	my $eq = $peak;
	my $cumul_pois = 0;
	for my $i($peak..@repeat_in){
		my $std0 = &ppois($i-1,$peak);
		my $std1 = &ppois($i,$peak);
		my $dif  = $std1 - $std0;
		my $pois = $repeat_in[$i-1]/$sum;
		if ($pois < $dif){
			$eq = $i;
		}else {$cumul_pois += $pois;}
	}
	my $cumul_std = 1- &ppois ($eq,$peak);
	return $cumul_pois - $cumul_std;
}
## Genomeye
my @Ge_out; #Genomeye results
    $Ge_out[0]="Genomeye";
    open IN_Ge, "<$logfile" or die $!;
    while (<IN_Ge>){
	$Ge_out[1]=$1 if ($_=~/n_kmer:\s+(\S+)/); #sum kmer number 
	$Ge_out[2]=$1 if ($_=~/Coverage\s+:\s+(\S+)/); #coverage 
        $Ge_out[3]=$1 if ($_=~/Genome size1\s+:\s+(\S+)/); #Gsize
        $Ge_out[4]=$1 if ($_=~/heterozygous\s+rate\s+:\s+(\S+)/); #heterozygous rate
    }
    close IN_Ge;
    &output(@Ge_out);

sub output {
	my @out = @_;
	$out[5] = $out[3]*$acc_rate; #revised Gsize 
	$out[6] = $out[5]/(4**17); #probability having the same kmer
          my $x = $bases/$out[5];
	my $out_str="$out[0]\tk$kmer\t$out[1]\t$out[2]\t$out[3]\t$out[5]\t$out[4]\t$rep\t$out[6]\t$acc_rate\t$bases\t$x";
	my @out_RS = split (/\s+/,$out_str);
	for my $i(0..@out_RS-1){
		next if ($out_RS[$i]=~/^\D+/); #filter non numeric element
		if ($i==3 && $out_RS[$i]=~/[\.]/){$out_RS[$i]=$out_RS[$i]} #[2]peak
		if ($i==4||$i==5) {$out_RS[$i]=$out_RS[$i]/1000000}#Gsize,Revised Gsize
        #[6]Hr,[7]rep,[8]probability,[9]acc_rate
		if ($i==6||$i==7||$i==8||$i==9) {
			$out_RS[$i]=$out_RS[$i]*100;
		} 
		if ($i==10){$out_RS[$i]=$out_RS[$i]/1000000000} #[10]Used bases
		if ($i==11){$out_RS[$i]=$out_RS[$i]} #[11]X
	}
	my $out_RS_str = join "\t",@out_RS;
	#print OUT "$out_str\n";
	print OUT "$out_RS_str\n";
}

close OUT;
#close OUT_RS;
}elsif(@ARGV ==3){
my ($file,$bases,$output)=@ARGV;
open (OUT,">$output.summary.xls") || die $!;
open (OUT2,">$output.frequency.xls") || die $!;
print OUT2 "k-mer Depth\tFrequency\n";
print OUT "Method\tK-mer\tK-mer individual number\tk-mer Depth\tGenome size(M)\tRevised genome size(M)\tHeterozygous rate\tRepeat\tProbability having the same kmer\tAccuracy rate\tUsed Bases\tAverage Sequencing Depth\n";
open (IN,$file) || die $!;
my (@kn,@ksp,@kn_out,@ksp_out);
my $more_kn =0;
my $sum =0;
my $sum_ksp =0;
my $kmer=17;
my %hash;
while(<IN>){
  chomp;
  my @arr = split(/\s+/,$_);
  my $numb = $arr[0]*$arr[1]; #caculate kmer individuals number 
  $sum += $numb;
  $sum_ksp += $arr[1];
  if ($arr[0]<1001) {
     $hash{$arr[0]}=$numb;
    push (@kn,$numb);
    push (@ksp,$arr[1]);
  }else{
    $more_kn += $numb;
  }
}
close IN;
foreach (sort{$a<=>$b} keys %hash){
  my $de=$hash{$_}/$sum;
  my $des=sprintf("%.3f",$de*100);
  print OUT2 "$_\t$des\n";
}
push (@kn,$more_kn);
my $acc_rate=1-$kn[0]/$sum; #accuracy rate (consider kmer_depth=1 is error number)
my $base=$sum*17;
@kn_out  = &est_freq_Hr1(1,@kn); #kn_out[5] is max_num
my $rep = &est_repeat1 ($kn_out[2],@kn);
&output1 (@kn_out);
@ksp_out = &est_freq_Hr1(0,@ksp); #kmer species results
&output1 (@ksp_out);

sub est_freq_Hr1{
	my ($flag,@in) = @_;
	my ($max_freq,$half_freq,$max_mean,$half_mean,$Hr);
	my $max_num = 0;
    ## compare to find the max freq
	for my $i(1..@in-2){
		if ($in[$i-1] < $in[$i] && $in[$i] > $in[$i+1]){
			if ($max_num < $in[$i]){
				$max_num = $in[$i];
				$max_freq = $i+1;
			}
		}
	}
	
	my $Gsize = $sum/$max_freq;

    ## estimate the average of heterozygous and homozygous peak
	if ($max_freq %2 ==0) {
		$half_freq = $max_freq/2;
		$half_mean = &mean1($half_freq,@in);
	}else {
		$half_freq = ($max_freq-1)/2;
		$half_mean = &mean1($half_freq,@in);
	}	
	$max_mean = &mean1($max_freq,@in); 

    ## Calculate heterozygous rate
	if ($flag == 1){
		$Hr = $half_mean/$kmer/($half_mean + $max_mean); #kn Hr
		return ("kn",$sum,$max_freq,$Gsize,$Hr,$max_num);
	}else {
		$Hr = $half_mean/$kmer/($half_mean + 2*$max_mean); #ksp Hr
		return ("ksp",$sum,$max_freq,$Gsize,$Hr,$max_num);
  }	
}

## calculate mean of three numbers
sub mean1{
	my ($mid,@mean_in) = @_;
	my $avg = ($mean_in[$mid-2]+$mean_in[$mid-1]+$mean_in[$mid])/3;
	return $avg;
}

## calculate repeat using area
sub est_repeat1{
	my ($peak,@repeat_in) = @_;
	my $eq = $peak;
	my $cumul_pois = 0;
	for my $i($peak..@repeat_in){
		my $std0 = &ppois($i-1,$peak);
		my $std1 = &ppois($i,$peak);
		my $dif  = $std1 - $std0;
		my $pois = $repeat_in[$i-1]/$sum;
		if ($pois < $dif){
			$eq = $i;
		}else {$cumul_pois += $pois;}
	}
	my $cumul_std = 1- &ppois ($eq,$peak);
	return $cumul_pois - $cumul_std;
}

sub output1 {
	my @out = @_;
	$out[5] = $out[3]*$acc_rate; #revised Gsize 
	$out[6] = $out[5]/(4**17); #probability having the same kmer
          my $x = $bases/$out[5];
	my $out_str="$out[0]\tk$kmer\t$out[1]\t$out[2]\t$out[3]\t$out[5]\t$out[4]\t$rep\t$out[6]\t$acc_rate\t$bases\t$x";
	my @out_RS = split (/\s+/,$out_str);
	for my $i(0..@out_RS-1){
		next if ($out_RS[$i]=~/^\D+/); #filter non numeric element
		if ($i==3 && $out_RS[$i]=~/[\.]/){$out_RS[$i]=$out_RS[$i]} #[2]peak
		if ($i==4||$i==5) {$out_RS[$i]=$out_RS[$i]/1000000}#Gsize,Revised Gsize
        #[6]Hr,[7]rep,[8]probability,[9]acc_rate
		if ($i==6||$i==7||$i==8||$i==9) {
			$out_RS[$i]=$out_RS[$i]*100;
		} 
		if ($i==10){$out_RS[$i]=$out_RS[$i]/1000000000} #[10]Used bases
		if ($i==11){$out_RS[$i]=$out_RS[$i]} #[11]X
	}
	my $out_RS_str = join "\t",@out_RS;
	#print OUT "$out_str\n";
	print OUT "$out_RS_str\n";
}
close OUT;
}


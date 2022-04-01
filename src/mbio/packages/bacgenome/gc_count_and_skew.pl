#!/usr/bin/perl -w
use strict;
use warnings;
use POSIX;

die "Usage: perl $0 [scaffold.fna] [windows] [step]\n" unless(@ARGV==3);
my $in = shift;
my $wind_num = shift;
my $step=shift;
open IN,"<$in"||die $!;
my %hash;my %sort;
$/=">";<IN>;
my $n= 1;
while(<IN>){
	chomp;
	my @arr = split /\n/,$_;
	my $sequence = join("",@arr[1..$#arr]);
	$hash{$arr[0]}=$sequence;
	$sort{$n} = $arr[0];
	$n++;
}
close IN;
$/="\n";

my $gc_sum=0;my $len_sum=0;my %gc;my %skew;my @sort;my %wind_sort;
for my $k (sort {$a <=> $b} keys %sort){
	my $seq = $hash{$sort{$k}};
	push (@sort,$sort{$k});
	my $id = $sort{$k};
	my $len = length($seq);
	$len_sum = $len_sum+$len;
	my $gc =($seq =~tr/GCgc/GCgc/);
	$gc_sum=$gc_sum+$gc;
	my $max = POSIX::ceil(1 + ($len-$wind_num)/$step);
#	my $max = POSIX::ceil($len/$wind_num);
	for(my $i=0;$i<$max;$i++){
		if ($i*$step+$wind_num <= $len){
			my $wind = substr($seq,$i*$step,$wind_num);
			my $g_wind = ($wind =~ tr/Gg/Gg/);
			my $c_wind = ($wind =~ tr/Cc/Cc/);
			my $gc_skew;
			if ($g_wind == 0 && $c_wind == 0){
				$gc_skew = 0;
			}else{
				$gc_skew = ($g_wind-$c_wind)/($g_wind+$c_wind);
			}
			my $gcwind_count = ($g_wind+$c_wind)/$wind_num;
			my $key = ($i*$step+1)." ".($i*$step+$wind_num);
			$gc{$id}{$key}=$gcwind_count;
			$skew{$id}{$key}=$gc_skew;
			push (@{$wind_sort{$id}},$key);
		}elsif($i*$step+$wind_num > $len){
			my $last = $len-$i*$step;
			my $wind_last = substr($seq,$i*$step,$last);
			my $g_wind_last = ($wind_last =~ tr/Gg/Gg/);
			my $c_wind_last = ($wind_last =~ tr/Cc/Cc/);
			my $gc_skew_last;
			if($g_wind_last == 0 && $c_wind_last == 0){
				$gc_skew_last = 0;
			}else{
				$gc_skew_last = ($g_wind_last-$c_wind_last)/($g_wind_last+$c_wind_last);
			}
			my $gcwind_count_last = ($g_wind_last+$c_wind_last)/$last;
			my $key_last = ($i*$step+1)." ".($i*$step+$last);
			$gc{$id}{$key_last}=$gcwind_count_last;
			$skew{$id}{$key_last}=$gc_skew_last;
			push (@{$wind_sort{$id}},$key_last);
		}
	}
}

open OUT1,">positive_gc_skew.txt" || die $!;
open OUT2,">negative_gc_skew.txt" || die $!;
open OUT3,">positive_gc_count.txt" || die $!;
open OUT4,">negative_gc_count.txt" || die $!;
my $gc_averge = $gc_sum/$len_sum;
for my $k1(@sort){
	for my $k2(@{$wind_sort{$k1}}){
		if ($skew{$k1}{$k2}>=0){
			print OUT1 "$k1 $k2 $skew{$k1}{$k2}\n";
		}else{
			print OUT2 "$k1 $k2 $skew{$k1}{$k2}\n";
		}
		my $temp = $gc{$k1}{$k2}-$gc_averge;
		if($temp >= 0){
			print OUT3 "$k1 $k2 $temp\n";
		}else{
			print OUT4 "$k1 $k2 $temp\n";
		}
	}
}
close OUT1;
close OUT2;
close OUT3;
close OUT4;

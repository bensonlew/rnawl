#!usr/bin/perl
use strict;
use warnings;
die "usage: perl $0 binRout\n" unless(@ARGV == 1);
my $file = $ARGV[0];#annotation rout for each bin
my %allScore;
my %levelScore;
my $t = 0.7;      #choose annotation result which score is higher than t value
my %levelName;    #每一水平的哈希，存放各级中的名字为数组
my %name2name;    #以上一级为key，存放下一集的名称为数组
my @class = ("none","none","d_","p_","c_","o_","f_","g_","s_");   #taxon level:domain is index 2
my ($level,$string,$subStr);
#store Score and taxon 
open IN,"$file/result.phylotype" or die "Can't find input file $!\n";
while(<IN>){
	next if (/Query/);
	chomp;
	my @tmp = split /\t/,$_;
	for my $flag(2..$#tmp){
		$tmp[$flag] =~ s/\)//;
		my @tmp2 = split(/\(/,$tmp[$flag]);
		if($tmp2[1] >= 0.7){
			$levelScore{$flag}{$tmp2[0]} += $tmp2[1];
			$allScore{$flag} += $tmp2[1];
			unless (grep(/^$tmp2[0]$/,@{$levelName{$flag}})){
				push @{$levelName{$flag}},$tmp2[0];
			}
		}
		if($tmp[$flag+1]){
			my $childName = (split(/\(/,$tmp[$flag+1]))[0];
			my $childScore = (split(/\(/,$tmp[$flag+1]))[1];
			$childScore =~ s/\)//;
			push @{$name2name{$class[$flag].$tmp2[0]}},$childName unless (grep(/^$childName$/,@{$name2name{$class[$flag].$tmp2[0]}}) || $childScore < 0.7);
		} 
	}
}
close IN;
#write result
open OT,">$file/bin.anno.xls" or die "$!\n";
print OT "Superkingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\n";
for my $nameFlag(0..$#{$levelName{2}}){   #ergodic domain name
	$level = 2;
	func($levelName{2}[$nameFlag],$level);
}
close OT;

#function:calculate score percent,then write them into anno.xls
sub func{
	my @in = @_;
	$level = pop @in;
	for my $i(@in){
		my $score = ($levelScore{$level}{$i})/$allScore{$level};
		$score = sprintf "%.2f",$score;
		$string .= $i."($score)\t";
		if(ref($name2name{$class[$level].$i}) eq "ARRAY" && @{$name2name{$class[$level].$i}} > 0){
			func(@{$name2name{$class[$level].$i}},$level+1);
		}else{
			my $test = $class[$level].$i;
			print OT "$string\n";
		}
        	$string = (split(/$i/,$string))[0];
	}
	$level -= 1;
}

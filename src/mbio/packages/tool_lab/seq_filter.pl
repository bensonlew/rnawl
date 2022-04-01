#! /usr/bin/env perl
use warnings;
use strict;

if(@ARGV < 4) {
    print STDERR "fasta-rename.pl  fastaFile sort(Y/N)\tminLen\tmaxLen\toutfile\n";
    exit;
}
my $file = shift;
my $sort = shift;
my $minLen = shift;
my $maxLen = shift;
my $outfile = shift;
open(FILE,"<$file") or die;
my $num =1;
my $fas;
my %seq;
my %len;
my %index;
$fas = <FILE>;
chomp $fas;
while(<FILE>){
	chomp $_;
	if($_ =~/\>/){
		$index{$_} = $num;
		$num++;
		$fas = $_;
	}else{
		$seq{$fas} .= $_;
		$len{$fas} = length($seq{$fas});
	}
}
close(FILE);

open(OUT,">$outfile") or die;
my $index = 0;
if($sort eq "Y"){
	foreach my $fas (sort {$len{$b}<=>$len{$a}} keys %len ){
		if($len{$fas}>=$minLen && $len{$fas}<=$maxLen){
			$index++;
			print OUT "$fas\n",uc($seq{$fas}),"\n";
		}
	}
}else{
	foreach my $fas (sort {$index{$a}<=>$index{$b}} keys %len){
		if($len{$fas}>$minLen){
			$index++;
			print OUT "$fas\n",uc($seq{$fas}),"\n";
		}
	}
}
close(OUT);

#!/usr/bin/perl -w
use strict;

die "usage: perl $0 fastq.file min.length output\n" unless(@ARGV==3);

my ($fq,$min,$out)=@ARGV;

if($fq=~/\.gz$/){
	open IN,"gzip -dc $fq |" or die "read $fq: $!\n";
}else{
	open IN,$fq or die "read $fq: $!\n";
}
open OT,">$out" or die "write $out: $!\n";

while(<IN>){
	chomp;
	my $id=(split /\s+/)[0];
	chomp(my $l2=<IN>);
	my $l3=<IN>;
	my $l4=<IN>;
	my $len=length($l2);
	if($len >= $min){
		print OT "$id\n$l2\n$l3$l4";
	}
}
close IN;
close OT;

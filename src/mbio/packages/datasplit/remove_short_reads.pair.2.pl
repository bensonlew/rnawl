#!/usr/bin/perl -w
use strict;

die "usage: perl $0 fastq1 fastq2 min.length out.prefix\n" unless(@ARGV==4);

my ($fq1,$fq2,$min,$out)=@ARGV;

if($fq1=~/\.gz$/){
	open I1,"gzip -dc $fq1 |" or die "read $fq1: $!\n";
}else{
	open I1,$fq1 or die "read $fq1: $!\n";
}
if($fq2=~/\.gz$/){
	open I2,"gzip -dc $fq2 |" or die "read $fq2: $!\n";
}else{
	open I2,$fq2 or die "read $fq2: $!\n";
}

open O1,">$out.1.fq" or die "write $out.1.fq: $!\n";
open O2,">$out.2.fq" or die "write $out.2.fq: $!\n";

while(<I1>){
	chomp;
	my $lid=(split /\s+/)[0];
	chomp(my $l2=<I1>);
	my $l3=<I1>;
	my $l4=<I1>;
	my $llen=length($l2);
	chomp(my $y1=<I2>);
	my $yid=(split /\s+/,$y1)[0];
	chomp(my $y2=<I2>);
	my $y3=<I2>;
	my $y4=<I2>;
	my $ylen=length($y2);
	if($llen >= $min && $ylen >= $min){
		print O1 "$lid\n$l2\n$l3$l4";
		print O2 "$yid\n$y2\n$y3$y4";
	}
}
close I1;
close I2;
close O1;
close O2;

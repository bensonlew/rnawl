#! /usr/bin/perl -w
use strict;
use warnings;

if(@ARGV < 2) {
    print STDERR "contigONscaffold.pl  scaffoldFile prefix\n";
    exit;
}
my $file = shift;
my $prefix = shift;
open(FILE,"<$file") or die;
my $query;
my %seq;
my @array;
while(<FILE>){
	chomp $_;
	if($_ =~ /^>(\S+)/){
		$query = $1;
		$array[scalar(@array)] = $1;
	}else{
		$seq{$query} .= $_;
	}
}
close(FILE);
my $order_out=0;
open(AGP,">$prefix".".agp") or die;
open(CONTIG,">$prefix".".scaf2contig") or die;
foreach my $fas (@array){
	my $row_in=0;
	my $pos=1;
	while($seq{$fas} =~ /[Nn]{1,}/g){
		$row_in++;
		$order_out++;
		print CONTIG ">contig$order_out\n".substr($seq{$fas},$pos-1,pos($seq{$fas})-length($&)-$pos+1)."\n";
		print AGP $fas."\t".$pos."\t".(pos($seq{$fas})-length($&))."\t".$row_in."\tW\tcontig".$order_out."\t1\t".(pos($seq{$fas})-length($&)-$pos+1)."\t+\n";
		$row_in++;
		print AGP $fas."\t".(pos($seq{$fas})-length($&)+1)."\t".pos($seq{$fas})."\t".$row_in."\tN\t".length($&)."\tscaffold\tyes\tpaired-ends\n";
		$pos = pos($seq{$fas})+1;
	}
	$order_out++;
	$row_in++;
	print CONTIG ">contig$order_out\n".substr($seq{$fas},$pos-1,length($seq{$fas})-$pos+1)."\n";
	print AGP $fas."\t".$pos."\t".length($seq{$fas})."\t".$row_in."\tW\tcontig".$order_out."\t1\t".(length($seq{$fas})-$pos+1)."\t+\n";
}

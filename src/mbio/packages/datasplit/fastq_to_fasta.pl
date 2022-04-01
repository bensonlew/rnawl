#!/usr/bin/perl -w

use strict;
use Getopt::Long;
my %opts;
GetOptions (\%opts,"i=s","o=s");
my $usage = <<"USAGE";
	Discription:change fq file to fa.
	Usage:perl $0 [options]
	-i    input.fq
	-o    output.fa
USAGE
die $usage if ( !($opts{i} && $opts{o}));

open IN,"<$opts{i}" or die "File $opts{i} not found error here\n";
open OUT,">$opts{o}" or die "File $opts{o} not found error here\n";
my $c = 0;
my $id ='';
my $seq = '';

while(<IN>){
	chomp;
	$c++;
	if($c == 1){
		$id = (split /\s+/,$_)[0];
		$id =~ s/\@//;
	}elsif($c == 2){
		$seq = $_;
	}elsif($c == 4){
		$c = 0;
		print OUT ">$id\n$seq\n";
	}
}
close IN;
close OUT;

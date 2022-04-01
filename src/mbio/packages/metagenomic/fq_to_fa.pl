#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
my %opts;
GetOptions (\%opts,"i=s","o=s");

my $usage = <<"USAGE";
	Usage:perl $0 [options] 
                -i	    input a fastq file;
                -o	    output a fasta file;

USAGE
die $usage if ( !$opts{i} && !$opts{o} );

open IN, "<$opts{i}"or die;
open OUT, ">$opts{o}"or die;

while(<IN>){
	if(/^@(.*)/){
     		chomp;
	     my	$id=">" . $1;
	     my	$seq=<IN>;
		print OUT "$id\n$seq";
}
}
close(IN);
close(OUT);

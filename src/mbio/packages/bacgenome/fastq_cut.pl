#! /usr/bin/perl
use strict;
use warnings;

if(@ARGV!=3) {
    print STDERR "fastq_cut.pl fastq front_num output\n";
    exit;
}
my ($file,$front,$output)=@ARGV;
open(FILE,"<$file") or die;
open(OUT,">$output") or die;

while(<FILE>){
	chomp $_;
	print OUT "$_\n";
	my $seq = <FILE>;
	chomp $seq;
	my $subSeq = substr($seq,0,$front);
	print OUT "$subSeq\n";
	my $qname= <FILE>;
	print OUT $qname;
	my $qual = <FILE>;
	chomp $qual;
	my $subQual = substr($qual,0,$front);
	print OUT "$subQual\n";	
	}	
	close(FILE);
close OUT;



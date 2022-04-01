#!/usr/bin/perl -w
#
#Author: Wangzhaoyue
use warnings;
use strict;
unless (@ARGV ==2) {
	&usage;
	exit;
}
open IN,$ARGV[0] or die "$ARGV[0] $!\n";
open OT,">$ARGV[1].scaftig" or die "$ARGV[1].scaftig $!\n";

my $min_length = 0;
my $name = '';
my $seq = '';

while(<IN>){
	if(/^>(\S+)/){
		&print_scafftig($name, $seq) if($seq);
		$name = $1;
		$seq  = '';
	} else {
		chomp;
		$seq .= $_;
	}
}
&print_scafftig($name, $seq) if($seq);

sub print_scafftig {
	my $name = shift;
	my $seq  = shift;
	my $id = 1;
	while($seq=~/([ATGCatgc]+)/g){
		my $s = $1;
		next if(length($s) < $min_length);
		print OT ">$name\_$id length=".length($s)."\n";
		while($s=~/(.{1,60})/g){
			print OT "$1\n";
		}
		$id ++;
	}
}
$/ = "\n";
close IN;
close OT;
sub usage{
	print STDERR<<USAGE;
	
	perl $0 [scafseq.file] [scaftig.file]
	
USAGE
}
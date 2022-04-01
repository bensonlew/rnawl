#!/usr/bin/perl -w
use strict;
use warnings;

if(@ARGV != 2){
	print STDERR "faa2m.pl <faa> <faa_m>\n";
	exit;
}

my $fas = shift;
my $faa = shift;

(-s $fas) || die "$fas is not exists!\n";
open (FAS, "<$fas") or die $!;
open (FAA, ">$faa") or die $!;
print FAA ">";
$/=">",<FAS>;
while(<FAS>) {
	my @all_lines = split /\n/,$_;
	$all_lines[0] =~ s/_1 / /;
	my $str=substr($all_lines[1],1,length($all_lines[1])-1);
	$all_lines[1] ="M$str";
	print FAA (join "\n",@all_lines);
}
print FAA "\n";
close FAS;
close FAA;

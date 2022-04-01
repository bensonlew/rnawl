#! /usr/bin/perl -w
use warnings;
use strict;

if(@ARGV < 4) {
    print STDERR "compare_reads_name.pl  file1.list file2.list same.list left2.list \n";
    exit;
}
my $file_1 = shift;
my $file_2 = shift;
my $same = shift;
#my $left1 = shift;
my $left2 = shift;
open(FILE1,"<$file_1") or die;
open(FILE2,"<$file_2") or die;
open(SAME,">$same ") or die;
#open(LEFT1,">$left1") or die;
open(LEFT2,">$left2") or die;
my %file1;
while(<FILE1>){
    chomp;
    my @q = split;
    $file1{$q[0]}=$q[1];
}
while(<FILE2>){
    chomp;
    my @a = split;
    if (exists $file1{$a[0]}){
        print SAME "$a[0]\n";
    }else{
        #print LEFT1 "$file1{$q[0]}\n";
        print LEFT2 "$a[0]\n";
    }
}
close FILE1;
close FILE2;
close SAME;
close LEFT2;
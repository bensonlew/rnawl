#!/usr/bin/perl -w
##author:qingchen.zhang
use strict;

die "usage: perl $0 input output\n" unless(@ARGV==2);

my ($in,$out)=@ARGV;

open IN,$in or die "read $in: $!\n";
open OUT,">$out" or die "write $out: $!\n";
while(<IN>){
    if(/^>/){
        chomp;
        my $id=(split /\s+/)[0];
        print OUT $id,"\n";
    }else{
        print OUT $_;
    }
}
close IN;
close OUT;
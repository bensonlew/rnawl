#!/usr/bin/perl

use strict;
use warnings;

die "perl $0 <input1> <output>" if ($#ARGV < 1);
open IN1,$ARGV[0] or die "$!";
#my $file = (split /.gbk/,$ARGV[0])[0];
open OUT,">$ARGV[1]" or die "$!";

$/= "  gene";
<IN1>;

my $n = 0;
my @ids ;
while (<IN1>){
        chomp;
        $_ =~ s/\n\s+/ /g;
        if ($_ =~ /.*\/product="(.*)"\s+\/protein_id="(\S+)".*\/translation="(.*)"/){
            my $des = $1;
            my $protein_id = $2;
            unless(grep {$_ eq $protein_id} @ids){
            push @ids, $protein_id;
            my $seq = $3;
            $seq =~ s/\s//g;
            print OUT ">$protein_id $des\n$seq\n";
            }else{
            print "$protein_id\n";
            }

        }else{
            #print "$_\n";
            $n++;
            #print "$n\n";
        }
}
close IN1;
close OUT;

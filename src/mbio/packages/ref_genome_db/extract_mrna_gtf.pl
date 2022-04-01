#!/usr/bin/perl
use strict;
use warnings;

die "Usage: perl $0 [total.gtf] [gene_biotype] [trans_biotype] [mrna.gtf]\n" unless @ARGV==4;

open FILEA, $ARGV[0];
open FILEB, $ARGV[1];
open FILED, $ARGV[2];
open FILEC, ">$ARGV[3]";

my %hash;
my %hash_trans;

while (<FILEB>){
    chomp;
    my @f1 = split /\t/;
    if ($f1[1] eq 'mRNA') {
        $hash{$f1[0]} = 'mRNA';
    }
}

while (<FILED>){
    chomp;
    my @f1 = split /\t/;
    if ($f1[1] eq 'mRNA') {
        $hash_trans{$f1[0]} = 'mRNA';
    }
}

while (<FILEA>){
    if ($_ =~ /^#/){
        print FILEC $_;
    }else{
        my @f1 = split /\t/;
#        if ($f1[8] =~ /gene_id "(.*?)";/){
#            #print $1, "\n";
#            if (exists $hash{$1}){
#                print FILEC $_;
#            }
#        }
        if ($f1[8] =~ /transcript_id "(.*?)";/){
            #print $1, "\n";
            if (exists $hash_trans{$1}){
                print FILEC $_;
            }
        }
    }
}

close FILEA;
close FILEB;
close FILEC;
close FILED;

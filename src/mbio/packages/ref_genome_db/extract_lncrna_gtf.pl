#!/usr/bin/perl
use strict;
use warnings;

die "Usage: perl $0 [total.gtf] [gene_biotype] [trans_biotype] [lncRNA.gtf] [lncRNA.list]\n" unless @ARGV==5;

open FILEA, $ARGV[0];
open FILEB, $ARGV[1];
open FILED, $ARGV[2];
open FILEC, ">$ARGV[3]";
open FILEE, ">$ARGV[4]";

my %hash;
my %hash_trans;

while (<FILEB>){
    chomp;
    my @f1 = split /\t/;
    if ($f1[1] eq 'lncRNA') {
        $hash{$f1[0]} = 'lncRNA';
    }
}

while (<FILED>){
    chomp;
    my @f1 = split /\t/;
    if ($f1[1] eq 'lncRNA') {
        $hash_trans{$f1[0]} = 'lncRNA';
        print FILEE $f1[0], "\n";
    }
}

while (<FILEA>){
    if ($_ =~ /^#/){
        print FILEC $_;
    }else{
        my @f1 = split /\t/;
        if ($f1[8] =~ /gene_id "(.*?)";/){
            #print $1, "\n";
            if (exists $hash{$1}){
                print FILEC $_;
            }
        }elsif ($f1[8] =~ /transcript_id "(.*?)";/){
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
close FILEE;

#!/usr/bin/perl -w
use strict;
use warnings;

die "usage: perl $0 bwa.sam sam_type fq_filter prefix\n" unless(@ARGV == 4);
my ($sam,$sam_type,$fq_filter,$prefix) = @ARGV;

my (%unalignSeq);
my $unalignNum = 0;
my $unalignBas=0;
my $alignNum = 0;
my @arr=split /\//,$prefix;
my $path=join "/",@arr[0..($#arr-1)];
my $sample=$arr[-1];
open INS,$sam or die "$!\n";
open LIS,">>$path/rehost_stat.list.txt"or die "$!\n";
if ($sam_type =~/pe/){
    open OUT1,"> $prefix.1.fq" or die "$!\n";
    open OUT2,"> $prefix.2.fq" or die "$!\n";
    while(<INS>){
        chomp;
        next if(/^@/);
        my @R1 = split /\t/;
        my $fq2_line;
        chomp($fq2_line = <INS>);
        my @R2 = split /\t/, $fq2_line;
        if ($fq_filter eq "map" and  ($R1[2] ne "*" or $R2[2] ne "*")){
            $alignNum =$alignNum + 2;
            print OUT1 "@",$R1[0],"\n",$R1[9],"\n","+","\n",$R1[10],"\n";
            print OUT2 "@",$R2[0],"\n",$R2[9],"\n","+","\n",$R2[10],"\n";
        }elsif( $R1[2] eq "*" and $R2[2] eq "*" and $fq_filter ne "map"){
            $unalignNum= $unalignNum +2;
            $unalignBas= $unalignBas+length($R1[9])+length($R2[9]);
            print OUT1 "@",$R1[0],"\n",$R1[9],"\n","+","\n",$R1[10],"\n";
            print OUT2 "@",$R2[0],"\n",$R2[9],"\n","+","\n",$R2[10],"\n";
        }
    }
    print LIS $sample,"\t",$unalignNum,"\t",$unalignBas,"\tpe\n";
    close OUT1;
    close OUT2;
}

if ($sam_type =~/se/){
    open OUTS1,"> $prefix.s.fq" or die "$!\n";
    while(<INS>){
        chomp;
        next if(/^@/);
        my @Rs = split /\t/;
        if ($Rs[2] ne "*" and $fq_filter eq "map"){
            $alignNum++;
            print OUTS1 "@",$Rs[0],"\n",$Rs[9],"\n","+","\n",$Rs[10],"\n";
        }elsif($Rs[2] eq "*" and $fq_filter ne "map"){
            $unalignNum++;
            $unalignBas= $unalignBas+length($Rs[9]);
            print OUTS1 "@",$Rs[0],"\n",$Rs[9],"\n","+","\n",$Rs[10],"\n";
        }
    }
    print LIS $sample,"\t",$unalignNum,"\t",$unalignBas,"\tse\n";
    close OUTS1;
}
close LIS;
close INS; 

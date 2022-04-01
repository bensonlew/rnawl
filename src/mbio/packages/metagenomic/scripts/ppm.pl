#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my %opts;
GetOptions(
        \%opts,"i=s","stat=s","o=s","help!");

my $usage = <<"USAGE";
              Version : 20180926
              Writer : shaohua.yuan\@majorbio.com

Usage: perl $0 [options]
                -i     *gene abundance file
                -stat  *clean reads stat
                -o     *output
                -help  Display this usage information
         eg:perl $0  -i reads_number.xls -stat reads.cleanData.stat -o out

USAGE

die $usage if ( !( $opts{i} ) || !( $opts{stat} ) || !( $opts{o} ));

#`mkdir -p $opts{o}` unless -e ($opts{o});

my (%sample_stat_dict);

open (INF,$opts{stat}) || die "can not open $opts{stat} !\n";
while(<INF>){
    chomp;
    next if (/^#/);
    my @tmp = split /\t/;
    my $sample = $tmp[0];
    my $readnum = $tmp[1];
    $sample_stat_dict{$sample} = $readnum;
}
close INF;

my (@sams);
open (ING,$opts{i}) || die "can not open $opts{i} !\n";
open OUT,"> $opts{o}  " or die "can not open $opts{o}!\n";
while(<ING>){
    chomp;
    if (/^GeneID/){
        my @head = split /\t/;
        @sams = split /\t/;
        @sams = grep { $_ ne "GeneID" } @sams;
        @sams = grep { $_ ne "Total" } @sams;
        print join(",",@sams),"\n";
        my $new_head = join("\t",@head);
        print OUT "$new_head\n";
    }else{
        my @tmp = split /\t/;
        my $gene = $tmp[0];
        my $gene_total;
        $gene_total = 0;
        print OUT "$gene";
        for(my $i = 0;$i < @sams ;$i++){
            my $sam = $sams[$i];
            my $sam_abu = $tmp[$i+1];
            #print "$sam_abu\t$sam\n";
            my $ppm_abu = ($sam_abu * (10**6))/ $sample_stat_dict{$sam} ;
            print OUT "\t$ppm_abu";
            $gene_total += $ppm_abu;
        }
        print OUT "\t$gene_total\n";
    }
}
close ING;
close OUT;

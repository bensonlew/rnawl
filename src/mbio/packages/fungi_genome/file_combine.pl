#!/usr/bin/perl
use strict;
use Getopt::Long;
my ($selrank,$row,$symbol,$output,$help);
GetOptions(
        "l:s"=>\$selrank,
        "r"=>\$row,
        "s:s"=>\$symbol,
        "h"=>\$help,
        "o:s"=>\$output
);
if(!@ARGV || $help || !$selrank){
        die"Name: file_combine.pl
Describe: to combine select ranks or rows in seval file to one file
Author: liuwenbin, liuwenbin\@genomic.org.cn
Version: 1.0, Date: 2011-03-19
Usage: perl file_com.pl <file1 file2 ... filen> [-option]
        files...    input files to combine
        -l <str>    selrank: f1:r1,r2..;f2:r1,r2,..;fn:ri,rj.., no default
        -r          sel row form infile, default sel rank
        -s <str>    symbol name add at eatch rank head
        -o          result of file
        -h          output help information to screen\n";
}
open (OUT,">$output") || die $!;
$symbol && (($symbol =~ s/,/\t/g),print OUT  "$symbol\n");
my @ranks = read_data($selrank,$row);#sub1
my $rankn = 0;
foreach(@ranks){
  my $temn = @{$_};
  ($temn > $rankn) && ($rankn = $temn);
}
foreach my $i(0..$rankn-1){
  my $out = ${$ranks[0]}[$i];
  foreach my $j(1..$#ranks){
    $out .= (${$ranks[$j]}[$i] ? "\t${$ranks[$j]}[$i]" : "\t0");
  }
  print OUT  "$out\n";
}

#sub1
#============#
sub read_data
#============#
{
    #usage: read_data(rank,row)
    my ($rank, $row) = @_;
    my @out_array;
    foreach my $i (split /;/, $rank){
        my @a = split /:|,/, $i;
        my $j = shift @a;
        foreach my $k (@a){
            my @out;
            if ($row){
                my $l = `sed -n '$k,${k}p' $ARGV[$j]`;
                $l =~ s/^\s+//;
                @out = split /\s+/, $l;
            }else{
                chomp(@out = `awk '{print \$$k}' $ARGV[$j]`);
            }
            push @out_array,\@out;
        }
    }
    @out_array;
}

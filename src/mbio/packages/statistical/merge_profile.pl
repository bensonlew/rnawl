#!/usr/bin/perl -w
use List::Util qw/sum/;
use strict;

my $script_name = $0;

my $arg;
my $in1;
my $sams;
my $out;
my $fa;

while ($arg=shift) {
        if    ($arg eq "-h" ) { print_usage(); }
	elsif ($arg eq "-i" ) { $in1 = shift; }
	elsif ($arg eq "-f" ) { $fa = shift; }
	elsif ($arg eq "-s" ) { $sams = shift; }
	elsif ($arg eq "-o" ) { $out = shift; }
}

($in1 and $sams and $out and $fa) || print_usage();
my (@genename,%lengths);
`mkdir $out` unless(-e "$out");

if($fa=~/\.gz$/){
    open I1,"gzip -dc $fa |" or die "can't read $fa: $!\n";
}else{
    open I1,$fa or die "can't read $fa: $!\n";
}

$/=">";
<I1>;
my $i=0;
while(<I1>){
    $i++;
    chomp;
    my @list=split /\n/;
    my $name=(split /\s+/,(shift @list))[0];
    my $seq=join "",@list;
    my $len=length($seq);
    $lengths{$name}=$len;
    push @genename,$name;
}
close I1;
$/="\n";

open O1,">$in1/head.txt" or die "can't write $in1/head.txt: $!\n";
print O1 "GeneID\n";
foreach my $gen(@genename){
   print O1 "$gen\n";
}
close O1;

open F3,">$out/gene.uniGeneset.fa.length.txt" or die "can't write $out/gene.uniGeneset.fa.length.txt: $!\n";
print F3 "GeneID\tgene_length\n";
my @files = ("reads_number.xls","reads_number_relative.xls","reads_length_ratio_relative.xls","reads_length_ratio.xls","RPKM.xls","TPM.xls","PPM.xls");
my $n =1;
my @samps = split/,/,$sams;
my $parts = 100;
my $split = int(@samps / $parts);
foreach my $file(@files){
    my $paste = "$in1/head.txt";
    for ($i = 0; $i <= $split; $i += 1){
        my @s_samps = @samps[$i * $parts .. (($i + 1) * $parts -1)];
        foreach my $samp(@s_samps){
            $paste = $paste." $in1/$samp\_$file" if $samp;
        }
        print "$paste\n";
        `paste $paste > $out/tmp_$file`;
        `mv $out/tmp_$file $out/bak_$file`;
        $paste = "$out/bak_$file";
    }
    open F1,"$out/bak_$file" or die "can't read $out/bak_$file: $!\n";
    open F2,">$out/$file" or die "can't write$out/$file: $!\n";
    print F2 "GeneID\t".join("\t", @samps)."\tTotal\n";
    readline F1; # skip the first line
    while (my $line=<F1>){
    chomp($line);
    my @list = split /\t/,$line;
    my $sum = sum @list[1..$#list];
    if ($sum <= 0){
         next;
    }else{
         print F2 "$line\t$sum\n";
         if($n == 1){
         print F3 "$list[0]\t$lengths{$list[0]}\n";
         }
    }
    }
    $n += 1;
    close F1, close F2,close F3;
}

my $num = @samps + 2;
`head -n 1 $out/reads_number.xls > $out/top100_reads_number.xls`;
`head -n 1 $out/reads_number_relative.xls > $out/top100_reads_number_relative.xls`;
`awk '{if(NR!=1) print }' $out/reads_number.xls | sort -gr -k $num | head -n 100 >> $out/top100_reads_number.xls`;
`awk '{if(NR!=1) print }' $out/reads_number_relative.xls |sort -gr -k $num  | head -n 100 >> $out/top100_reads_number_relative.xls`;

sub print_usage {
        print <<EOD;
Usage: perl $script_name options
    -i input samples profile dir
    -s input sample names ,splited by ,
    -f fasta sequence
     -o output directory, required
     -h print this help

EOD
exit;
}
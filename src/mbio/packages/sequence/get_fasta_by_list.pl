#! /usr/bin/perl -w
use strict;
die "perl $0 <lst><fa><out> [faa]\n" unless  @ARGV>=3;
my ($lst,$fa,$out)=@ARGV;
my $faa = $ARGV[3];
open IN,$lst||die;
readline IN;
my %ha;
map{chomp;$ha{(split('\t'))[0]}=1}<IN>;
close IN;

$fa=~/gz$/?(open IN,"gzip -cd $fa|"||die):(open IN,$fa||die);
$/=">";<IN>;$/="\n";
my %out;
open O2,">$out" or die "$out: $!\n";
while(<IN>){
    my ($info, $posf);
    $posf = '';
    $info=$1 if(/^(\S+)/);
    ($info, $posf) = ($1, $2) if(/^(\S+)(\_\d)$/ and $faa);
    $/=">";
    my $seq=<IN>;
    $/="\n";
    $seq=~s/>|\r|\*//g;
    print O2">$info$posf\n$seq" if(exists $ha{$info} && ! exists $out{$info});
    $out{$info}=1;
}
close O2;
close IN;

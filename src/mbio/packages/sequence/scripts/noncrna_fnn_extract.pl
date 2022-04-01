#!/usr/bin/perl -w 

use strict;
use warnings;

(@ARGV==4) || die "Usage: perl $0 [input.fasta] [gff] [type] [outdir]\n";
my($fna,$gff,$type,$outdir)=@ARGV;
my ($line,@inf,%seq,%scaffold,%num,$start,$end,%gff,@scaff,%accumu);
my $accumulate = 0;

open FNA, "$fna" or die "can not open file: $fna\n";
$/=">",<FNA>;
while($line=<FNA>){
    chomp $line;
    @inf=split /\n/,$line;
    my @ele=split /\s/,$inf[0];
    my $temp="";
    for(my $i=1;$i<=$#inf;$i++){
        $temp.=$inf[$i];
    }
    my @scaffold_name = split / /,$inf[0];
    my $scaffold = $scaffold_name[0];
    push (@scaff,$scaffold);
    $accumu{$scaffold} = $accumulate;
    $accumulate = $accumu{$scaffold} + length($temp);
    $seq{$ele[0]}=$temp;
}
close FNA;

$/="\n";
open IN, "$gff" or die "can not open file: $gff\n";
my @arr = split /\//,$gff;
my @file = split /\.gff/,$arr[-1];
my ($fnn_name,$gbkff,$gff_new);
if ($type =~ /rrna/){
    $fnn_name = $file[0].".rRNA.fnn";
    $gbkff = $file[0].".rRNA.gbk.gff";
    $gff_new = $file[0].".rRNA.gff";
}elsif ($type =~ /trna/){
    $fnn_name = $file[0].".tRNA.fnn";
    $gbkff = $file[0].".tRNA.gbk.gff";
    $gff_new = $file[0].".tRNA.gff";
    <IN>;<IN>;<IN>;
}
open FNNOUT,">$outdir/$fnn_name" or die;
open GBKFF,">$outdir/$gbkff" or die;
open GFF,">$outdir/$gff_new" or die;
if ($type =~ /rrna/){
    print GFF "Sequence Name\tSequence id\tType\tStart\tEnd\tE-values\tStrand\tPhase\tAttributes\n";
     print GBKFF "Sequence Name\tSequence id\tType\tStart\tEnd\tE-values\tStrand\tPhase\tAttributes\tBegin(gbk)\tEnd(gbk)\n";
}elsif ($type =~ /trna/){
    print GFF "Sequence Name\tSequence id\tStart\tEnd\ttRNA Type\tAnti Codon\tIntron Begin\tBounds End\tScore\n";
    print GBKFF "Sequence Name\tSequence id\tStart\tEnd\ttRNA Type\tAnti Codon\tIntron Begin\tBounds End\tScore\tBegin(gbk)\tEnd(gbk)\n";
}

#my $sum = 0;
while(<IN>){
        chomp;
        next if ($_ =~ /^##gff/);
        if ($type =~ /trna/){
            @inf = split /\s+/,$_;
        }elsif ($type =~ /rrna/){
            @inf= split /\t/,$_;
        }
        my @scaffold_name = split /_/,$inf[0];
        my $scaffold = pop(@scaffold_name);
        print($scaffold);
        my $sample = join "_", @scaffold_name;
        print($sample);
        if (exists $num{$scaffold}){
            $num{$scaffold}++;
        }else{
            $num{$scaffold} = 1;
        }
        if ($type =~ /rrna/){
            my $temp=substr($seq{$inf[0]},$inf[3]-1,$inf[4]-$inf[3]+1);
            if($inf[6] eq "-"){
                $temp=reverse($temp);
                $temp=~ tr/atcgATCG/tagcTAGC/;
            }
            #$sum ++;
            #$inf[0] = $sample."_rRNA".("0"x(2-length($sum))).$sum;
            $inf[0] = $sample."_rRNA".("0"x(2-length($num{$scaffold}))).$num{$scaffold};
            $inf[1] = $scaffold."_rRNA".("0"x(2-length($num{$scaffold}))).$num{$scaffold};
            my @arr = split /[=;_]/,$inf[8];
            my $rrna_type = $arr[1];
		    $start = $accumu{$scaffold} + $inf[3];
		    $end = $accumu{$scaffold} + $inf[4];
            my $num = $num{$scaffold};
            $gff{$num}{$scaffold}{"gff"} = (join "\t",@inf)."\n";
            $gff{$num}{$scaffold}{"gbk"} = (join "\t",@inf)."\t$start\t$end\n";
            $gff{$num}{$scaffold}{"fnn"} = ">$inf[0] $start $end $inf[1] $inf[3] $inf[4] $rrna_type\n$temp\n";
        }elsif ($type =~ /trna/){
            my $temp = "";
            if ($inf[2] > $inf[3]){
                $temp=substr($seq{$inf[0]},$inf[3]-1,$inf[2]-$inf[3]+1);
                $temp=reverse $temp;
                $temp=~ tr/atcgATCG/tagcTAGC/;
            }else{
                $temp=substr($seq{$inf[0]},$inf[2]-1,$inf[3]-$inf[2]+1);
            }
            $start = $accumu{$scaffold} + int($inf[2]);
            $end = $accumu{$scaffold} + int($inf[3]);
            #$sum ++;
            #$inf[0] = $sample."_tRNA".("0"x(3-length($sum))).$sum;
            $inf[0] = $sample."_tRNA".("0"x(3-length($inf[1]))).$inf[1];
            my $num = $inf[1];
            $inf[1] = $scaffold."_tRNA".("0"x(3-length($inf[1]))).$inf[1];
            if ($inf[$#inf] =~ /^[0-9\.]{1,}/){
                my $ok = 0;
            }else{
                splice (@inf, -1);
            }
            $gff{$num}{$scaffold}{"gff"} = (join "\t",@inf)."\n";
            $gff{$num}{$scaffold}{"gbk"} = (join "\t",@inf)."\t$start\t$end\n";
            $gff{$num}{$scaffold}{"fnn"} = ">$inf[0] $start $end $inf[1] $inf[2] $inf[3]\n$temp\n";
        }
}

my $sum = 0;
for my $scaffold (@scaff){
    for my $k1 (sort { $gff{$a} <=> $gff{$b} } keys %gff){
        if (exists $gff{$k1}{$scaffold}){
            $sum ++;
            my @gff = split /\t/,$gff{$k1}{$scaffold}{"gff"};
            $gff[0] = (split /RNA/,$gff[0])[0]."RNA".("0"x(3-length($sum))).$sum;
            my @fnn = split / /,$gff{$k1}{$scaffold}{"fnn"};
            $fnn[0] = ">".(split /RNA/,$gff[0])[0]."RNA".("0"x(3-length($sum))).$sum;
            my @gbk = split /\t/,$gff{$k1}{$scaffold}{"gbk"};
            $gbk[0] = (split /RNA/,$gff[0])[0]."RNA".("0"x(3-length($sum))).$sum;
            print FNNOUT (join " ",@fnn);
            print GFF (join "\t",@gff);
            print GBKFF (join "\t",@gbk);
            }
    }
}
close IN;
close FNNOUT;
close GFF;
close GBKFF;

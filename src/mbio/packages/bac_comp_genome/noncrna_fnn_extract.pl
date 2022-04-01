#!/usr/bin/perl -w

##此脚本主要是用于生成标准格式的gff文件，修改基因前缀和基因组位置

use strict;
use warnings;

(@ARGV==6) || die "Usage: perl $0 [input.fasta] [gff] [type] [outdir] [genome_name] [gene_tag]\n";
my($fna,$gff,$type,$outdir, $genome_name, $gene_tag)=@ARGV;
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
    #my @scaffold_name = split /_/,$inf[0];
    my @scaffold_name = $ele[0];
    my $scaffold = pop(@scaffold_name);
    #print "$scaffold\n";
    push (@scaff,$scaffold); ###将原始序列的名称提取出来，并将原始序列分条进行统计
    $accumu{$scaffold} = $accumulate;
    $accumulate = $accumu{$scaffold} + length($temp);
    $seq{$ele[0]}=$temp;
}
close FNA;

$/="\n";
open IN, "$gff" or die "can not open file: $gff\n";
my @arr = split /\//,$gff;
my @file = split /\.gff/,$arr[-1];
my ($fnn_name,$gbkff,$gff_new, $s16);
if ($type =~ /rrna/){
    $fnn_name = $genome_name ."_rRNA.fna";
    $gbkff = $genome_name ."_rRNA.gbk.gff";
    $gff_new = $genome_name ."_rRNA.gff";
    $s16 = $genome_name ."_16S.ffn";
}elsif ($type =~ /trna/){
    $fnn_name = $genome_name ."_tRNA.fna";
    $gbkff = $genome_name ."_tRNA.gbk.gff";
    $gff_new = $genome_name ."_tRNA.gff";
    <IN>;<IN>;<IN>;
}
if ($type =~ /rrna/){
    open RNA,">$outdir/$s16" or die;
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
my (@rrna_s, $s16_num, %rrna_seq, %rrna_id);
while(<IN>){
        chomp;
        next if ($_ =~ /^##gff/);
        if ($type =~ /trna/){
            @inf = split /\s+/,$_; #因为tRNA预测结果中格式不固定可能会出现多个空格的情况并不只是\t的情况；
        }elsif ($type =~ /rrna/){
            @inf= split /\t/,$_;
        }
        #my @scaffold_name = split /_/,$inf[0]; #如果序列为draft类型，那么将以draft的每条序列编号为location信息
        #my $scaffold = pop(@scaffold_name);
        my $scaffold = $inf[0];#如果序列为draft类型，那么将以draft的每条序列编号为location信息
        #my $scaffold = $genome_name;
        #my $sample = join "_", @scaffold_name;
        my $sample = $gene_tag;
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
            $inf[0] = $sample.("0"x(2-length($num{$scaffold}))).$num{$scaffold};
            $inf[1] = $scaffold."_rRNA".("0"x(2-length($num{$scaffold}))).$num{$scaffold};
            my @arr = split /[=;_]/,$inf[8];
            my $rrna_type = $arr[1];
            $start = $accumu{$scaffold} + $inf[3];
            $end = $accumu{$scaffold} + $inf[4];
            my $num = $num{$scaffold};
            $gff{$num}{$scaffold}{"gff"} = (join "\t",@inf)."\n";
            $gff{$num}{$scaffold}{"gbk"} = (join "\t",@inf)."\t$start\t$end\n";
            $gff{$num}{$scaffold}{"fnn"} = ">$inf[0] $start $end $inf[1] $inf[3] $inf[4] $rrna_type\n$temp\n";
            ###下面用于提取16S用，按16S长度进行排序
            if ($rrna_type eq "16S"){
		        $s16_num=length($temp);
                push @rrna_s, $s16_num;
                $rrna_id{$inf[0]} = "$inf[0] $start $end $inf[1] $inf[3] $inf[4]";
                $rrna_seq{$inf[0]} = $temp;
		    }
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
            $inf[0] = $sample.("0"x(3-length($inf[1]))).$inf[1];
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
###这个函数作用是将多条scaffold或者多个基因组序列循环输出结果
for my $scaffold (@scaff){
    for my $k1 (sort { $gff{$a} <=> $gff{$b} } keys %gff){
        if (exists $gff{$k1}{$scaffold}){
            $sum ++;
            my @gff = split /\t/,$gff{$k1}{$scaffold}{"gff"};
            my @fnn = split / /,$gff{$k1}{$scaffold}{"fnn"};
            my @gbk = split /\t/,$gff{$k1}{$scaffold}{"gbk"};
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

#先生成一个大小经过排序的数组，然后遍历哈希，从哈希中取大小排序过后的序列及序列id
@rrna_s = sort {$b <=> $a} @rrna_s;
my @seq_id;
if ($type =~ /rrna/){
    foreach my $seq_size(@rrna_s){
        foreach my $key (sort keys %rrna_seq){
            if (length($rrna_seq{$key})==$seq_size){
                if (not grep {$key eq $_} @seq_id){
                    push @seq_id, $key;
                    print RNA "\>$rrna_id{$key}\n$rrna_seq{$key}\n";
                }
            }
        }
    }
    close RNA;
}
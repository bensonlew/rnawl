#!/usr/bin/perl -w

##�˽ű���Ҫ���������ɱ�׼��ʽ��gff�ļ����޸Ļ���ǰ׺�ͻ�����λ��

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
    push (@scaff,$scaffold); ###��ԭʼ���е�������ȡ����������ԭʼ���з�������ͳ��
    $accumu{$scaffold} = $accumulate;
    $accumulate = $accumu{$scaffold} + length($temp);
    $seq{$ele[0]}=$temp;
}
close FNA;

$/="\n";
open IN, "$gff" or die "can not open file: $gff\n";
my @arr = split /\//,$gff;
my @file = split /\.gff/,$arr[-1];
my ($s16);
if ($type =~ /rrna/){
    $s16 = $genome_name ."_16S.ffn";
}
if ($type =~ /rrna/){
    open RNA,">$outdir/$s16" or die;
}

#my $sum = 0;
my (@rrna_s, $s16_num, %rrna_seq, %rrna_id);
while(<IN>){
        chomp;
        next if ($_ =~ /^##gff/);
        if ($type =~ /trna/){
            @inf = split /\s+/,$_; #��ΪtRNAԤ�����и�ʽ���̶����ܻ���ֶ���ո���������ֻ��\t�������
        }elsif ($type =~ /rrna/){
            @inf= split /\t/,$_;
        }
        my $scaffold = $inf[0];#�������Ϊdraft���ͣ���ô����draft��ÿ�����б��Ϊlocation��Ϣ
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
            $gff{$num}{$scaffold}{"fnn"} = ">$inf[0] $start $end $inf[1] $inf[3] $inf[4] $rrna_type\n$temp\n";
            ###����������ȡ16S�ã���16S���Ƚ�������
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
            $gff{$num}{$scaffold}{"fnn"} = ">$inf[0] $start $end $inf[1] $inf[2] $inf[3]\n$temp\n";
        }
}


#������һ����С������������飬Ȼ�������ϣ���ӹ�ϣ��ȡ��С�����������м�����id
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
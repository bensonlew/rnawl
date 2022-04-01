#! /usr/bin/perl
use strict;
use warnings;

if(@ARGV != 6){
	print STDERR "get_fa_bygff.pl <genome_fna> <faa> <Predict> <out_dir> <type>\n";
	print STDERR "This script is designed to extract fnnFile for predicting_gff\n";
	print STDERR "The third parameter N_num is the number of N setting between contigs\n";
	exit;
}

my $fas = shift;#核酸序列--基因组水平上的
my $faa = shift; #蛋白序列
my $gene_fna = shift; #基因序列文件
my $predict = shift; #预测的gff文件
#my $N_num = shift; #N的个数
my $prefix = (split /\_/,(split /\//,$predict)[-1])[-1];
my $pre = (split /\./, $prefix)[0];
my $out_dir = shift; #输出文件夹
my $type = shift; #是complete、chromosome还是draft
my $name = (split /\//,$fas)[-1];
$name =~ s/\.fna$//;
#print $name, "\n";

(-s $faa) || die "$faa is not exists!\n";
open (FAA, "<$faa") or die $!;
my %faa_seq;
my $faa_name;
while(<FAA>){
    chomp;
    if(/^>(.*)$/){
    my @faa_tmp = split (/\s+/, $1);
        if($1 =~ /(.+)locus_tag=(.+)$/){
            my $faa_n = $2;
            my @faa_sp = split(/\]/, $faa_n);
            $faa_name = $faa_sp[0];
        }else{
            $faa_name = $faa_tmp[0];
        }
    }else{
        $faa_seq{$faa_name}.=$_;
    }
}
close FAA;

open (GENE, "<$gene_fna") or die $!;
my %gene_fna_seq;
my $gene_fna_name;
while(<GENE>){
    chomp;
    if(/^>(.*)$/){
        my @fna_tmp = split (/\s+/, $1);
        if($1 =~ /(.+)locus_tag=(.+)$/){
            my $fna_n = $2;
            #print "$fna_n\n";
            my @fna_sp = split(/\]/, $fna_n);
            $gene_fna_name = $fna_sp[0];

        }else{
            $gene_fna_name = $fna_tmp[0];
        }
    }else{
        $gene_fna_seq{$gene_fna_name}.=$_;
    }
}
close GENE;


(-s $fas) || die "$fas is not exists!\n";
open (FAS, "<$fas") or die $!;
my ($line,@inf,%seq,%rev_complement,%scaffold,%gff,@scaff,%accumu);
my $accumulate = 0;
$/=">",<FAS>;
while($line=<FAS>){
    chomp $line;
    @inf=split /\n/,$line;
    my @ele=split /\s/,$inf[0];
    my $temp="";
    for(my $i=1;$i<=$#inf;$i++){
        $temp.=$inf[$i];
    }
    my $scaffold = $ele[0];
    push (@scaff,$scaffold);
    $accumu{$scaffold} = $accumulate;
    $accumulate = $accumu{$scaffold} + length($temp);
    $seq{$scaffold}=$temp;
    $rev_complement{$scaffold} = reverse $temp;
    $rev_complement{$scaffold} =~ tr/agctAGCT/TCGATCGA/;
}
close FAS;


my ($astart,$aend,$seq,@rrna_s, $s16, %rrna_seq, %rrna_id,@aa, @arr);
open (PRE, "<$predict") or die $!; #基因预测的gff文件
$/="\n";
chomp(my $li_name = <PRE>);
open (FFN,">$out_dir"."/$name"."_$pre".".fna") or die $!;
open (FFA,">$out_dir"."/$name"."_$pre".".faa") or die $!;
open (RNA,">$out_dir"."/$name"."_16S".".fna") or die $!;
open (LIST,">$out_dir"."/$name"."_$pre".".gff") or die $!;
if ($pre eq "CDS"){
	print LIST "$li_name\tA.start\tA.end\tInitiator Codon\tTerminator Codon\n";
}else{
	print LIST "$li_name\tA.start\tA.end\n";
	}
while(<PRE>){
	chomp $_;
	@arr = split /\t/, $_;
	my @sequence_id = split /_/,$arr[1];
	my $scaffold = $sequence_id[0] . "_" . $sequence_id[1];
	if (exists $seq{$scaffold}){
        if ($type eq "complete"){
            #$arr[1] =~ s/plasmid/p/;
            $astart = $arr[2];
            $aend = $arr[3];
        }elsif($type eq "chromosome"){
            #$arr[1] =~ s/chromosome/chr/;
            $astart = $arr[2];
            $aend = $arr[3];
        }else{
            $astart = $accumu{$scaffold} + $arr[2];
		    $aend = $accumu{$scaffold} + $arr[3];
        }
        if($arr[2]<$arr[3] and $arr[5] eq "-"){
            $seq = substr($rev_complement{$scaffold},(length($seq{$scaffold})-$arr[3]),$arr[3]-$arr[2]+1);
        }elsif($arr[2]<$arr[3] and $arr[4] eq "-"){
            $seq = substr($rev_complement{$scaffold},(length($seq{$scaffold})-$arr[3]),$arr[3]-$arr[2]+1);
        }elsif($arr[2]<$arr[3] and $arr[5] eq "+"){
            $seq = substr($seq{$scaffold},$arr[2]-1,$arr[3]-$arr[2]+1);
        }elsif($arr[2]<$arr[3] and $arr[4] eq "+"){
            $seq = substr($seq{$scaffold},$arr[2]-1,$arr[3]-$arr[2]+1);
        }elsif($arr[2]<$arr[3]){
            $seq = substr($seq{$scaffold},$arr[2]-1,$arr[3]-$arr[2]+1);
        }else{
            $seq = substr($rev_complement{$scaffold},(length($seq{$scaffold})-$arr[2]),$arr[2]-$arr[3]+1);
        }
        $seq =~ tr/agct/AGCT/;
		my $initiator = substr($seq, 0 ,3);
		my $terminator = substr($seq, -3);
		if(($arr[2]<$arr[3] and $arr[5] =~ /^-$/) or ($arr[2]<$arr[3] and $arr[4] =~ /^-$/ and $pre != "rRNA")){
		if ($pre eq "CDS"){
		    if (exists $faa_seq{$arr[0]}){
		        print FFA "\>$arr[0] $astart $aend $arr[1] $arr[3] $arr[2]\n$faa_seq{$arr[0]}\n";
		    }
		    if (exists $gene_fna_seq{$arr[0]}){
		        print FFN "\>$arr[0] $astart $aend $arr[1] $arr[3] $arr[2]\n$gene_fna_seq{$arr[0]}\n";
		    }
		}else{
		    print FFN "\>$arr[0] $astart $aend $arr[1] $arr[3] $arr[2]\n$seq\n";
		}
		if ($pre eq "rRNA"){
		    if ($arr[-1] =~ /16S/){
		        $s16=length($seq);
                push @rrna_s, $s16;
                $rrna_id{$arr[0]} = "$arr[0] $astart $aend $arr[1] $arr[3] $arr[2]";
                $rrna_seq{$arr[0]} = $seq ;
		    }
		}
		my $ss_sart=$arr[2];
		my $ss_end=$arr[3];
		@aa=@arr;
		$aa[3]=$ss_sart;
		$aa[2]=$ss_end;
		if ($pre eq "CDS"){
		    print LIST (join "\t", @aa),"\t", $astart,"\t", $aend,"\t", $initiator,"\t", $terminator, "\n";
		}else{
		    print LIST (join "\t", @aa),"\t", $astart,"\t", $aend,"\t", "\n";
		    }
		}else{
		if ($pre eq "CDS"){
		    if (exists $faa_seq{$arr[0]}){
		        print FFA "\>$arr[0] $astart $aend $arr[1] $arr[2] $arr[3]\n$faa_seq{$arr[0]}\n";
		    }
		    if (exists $gene_fna_seq{$arr[0]}){
		        print FFN "\>$arr[0] $astart $aend $arr[1] $arr[2] $arr[3]\n$gene_fna_seq{$arr[0]}\n";
		    }

		}else{
		    print FFN "\>$arr[0] $astart $aend $arr[1] $arr[2] $arr[3]\n$seq\n";
		}
		if ($pre eq "rRNA"){
		    if ($arr[-1] =~ /16S/){
		        $s16=length($seq);
                push @rrna_s, $s16;
                $rrna_id{$arr[0]} = "$arr[0] $astart $aend $arr[1] $arr[2] $arr[3]";
                $rrna_seq{$arr[0]} = $seq ;
		    }
		}

		if ($pre eq "CDS"){
		    print LIST (join "\t", @arr),"\t", $astart,"\t", $aend,"\t", $initiator,"\t", $terminator, "\n";
		}else{
		    print LIST (join "\t", @arr),"\t", $astart,"\t", $aend,"\t", "\n";
		    }
		}
 	}
}
close PRE;
close FFN;
close FFA;

#先生成一个大小经过排序的数组，然后遍历哈希，从哈希中取大小排序过后的序列及序列id
@rrna_s = sort {$b <=> $a} @rrna_s;
my @seq_id;
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
#!/usr/bin/perl -w
use strict;
use warnings;

if(@ARGV != 7){
	print STDERR "get_fa_bygff.pl <fasta> <Predict>  <N_num> <out_dir> <genome_name> <type> <genome_type>\n";
	print STDERR "This script is designed to extract fnnFile for predicting_gff\n";
	print STDERR "The third parameter N_num is the number of N setting between contigs\n";
	exit;
}

my $fas = shift;
my $predict = shift;
my $N_num = shift;
my $out_dir = shift;
my $name = shift;
my $prefix = shift;
my $type = shift;

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

my ($astart,$aend,$seq, @rrna_s, $s16, %rrna_seq, %rrna_id,@sequence_id);
open (PRE, "<$predict") or die $!;
$/="\n";
chomp(my $li_name = <PRE>);
if ($prefix eq 'rRNA'){
    open (RNA,">$out_dir"."/$name"."_16S".".fna") or die $!;
}
open (FFN,">$out_dir"."/$name"."_$prefix".".fna") or die $!;
open (LIST,">$out_dir"."/$name"."_$prefix".".gff") or die $!;
print LIST "$li_name\tA.start\tA.end\tInitiator Codon\tTerminator Codon\n";
while(<PRE>){
	chomp $_;
	my @arr = split /\t/, $_;
	if ($prefix eq "CDS"){
    	@sequence_id = split /_ORF/,$arr[1];

    }elsif($prefix eq "rRNA"){
        @sequence_id = split /_rRNA/,$arr[1];
    }elsif($prefix eq "tRNA"){
        @sequence_id = split /_tRNA/,$arr[1];
    }
	my $scaffold = $sequence_id[0];
	if (exists $seq{$scaffold}){
        if ($type eq 'chromosome'){
            $astart = $arr[2];
            $aend = $arr[3];
        }elsif($type eq 'complete'){
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
		if(($arr[2]<$arr[3] and $arr[5] =~ /^-$/) or ($arr[2]<$arr[3] and $arr[4] =~ /^-$/ and $prefix != "rRNA")){
		print FFN "\>$arr[0] $astart $aend $arr[1] $arr[3] $arr[2]\n$seq\n";
            if ($prefix eq "rRNA"){
                if ($arr[-1] =~ /16S/){
                    $s16=length($seq);
                    push @rrna_s, $s16;
                    $rrna_id{$arr[0]} = "$arr[0] $astart $aend $arr[1] $arr[3] $arr[2]";
                    $rrna_seq{$arr[0]} = $seq ;
                }
            }
		my $ss_sart=$arr[2];
		my $ss_end=$arr[3];
		my @aa=@arr;
		$aa[3]=$ss_sart;
		$aa[2]=$ss_end;
		print LIST (join "\t",@aa),"\t", $astart,"\t", $aend,"\t", $initiator,"\t", $terminator, "\n";
		}else{
            if ($prefix eq "rRNA"){
                if ($arr[-1] =~ /16S/){
                    $s16=length($seq);
                    push @rrna_s, $s16;
                    $rrna_id{$arr[0]} = "$arr[0] $astart $aend $arr[1] $arr[3] $arr[2]";
                    $rrna_seq{$arr[0]} = $seq ;
                }
            }
		print FFN "\>$arr[0] $astart $aend $arr[1] $arr[2] $arr[3]\n$seq\n";
		print LIST (join "\t",@arr),"\t", $astart,"\t", $aend,"\t", $initiator,"\t", $terminator, "\n";
		}
 	}
}
close PRE;
close FFN;


#������һ����С������������飬Ȼ�������ϣ���ӹ�ϣ��ȡ��С�����������м�����id
@rrna_s = sort {$b <=> $a} @rrna_s;
my @seq_id;
if ($prefix eq "rRNA"){
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

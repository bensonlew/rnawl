#!/usr/bin/perl -w
use strict;
use warnings;

if(@ARGV != 4){
	print STDERR "get_fa_bygff.pl <fasta> <Predict>  <N_num> <out_dir>\n";
	print STDERR "This script is designed to extract fnnFile for predicting_gff\n";
	print STDERR "The third parameter N_num is the number of N setting between contigs\n";
	exit;
}

my $fas = shift;
my $predict = shift;
my $N_num = shift;
my $prefix = (split /\./,(split /\//,$predict)[-1])[0];
my $out_dir = shift;
my $name = (split /\//,$fas)[-1];
$name =~ s/\.(\S+)$//;

(-s $fas) || die "$fas is not exists!\n";
open (FAS, "<$fas") or die $!;
my ($line,@inf,%seq,%rev_complement,%scaffold,%gff,@scaff,%accumu);
my $accumulate = 0;
$/=">",<FAS>;
while($line=<FAS>){
    chomp $line;
    @inf=split /\n/,$line;
    my @ele=split /\s+/,$inf[0];
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

my ($astart,$aend,$seq);
open (PRE, "<$predict") or die $!;
$/="\n";
chomp(my $li_name = <PRE>);
open (FFN,">$out_dir"."/$name"."_$prefix".".fnn") or die $!;
open (LIST,">$out_dir"."/$name"."_$prefix".".gff") or die $!;
print LIST "$li_name\tA.start\tA.end\tInitiator Codon\tTerminator Codon\n";
while(<PRE>){
	chomp $_;
	my @arr = split /\t/, $_;
	my @sequence_id = split /_ORF/,$arr[1];
	my $scaffold = $sequence_id[0];
	print $scaffold;
	if (exists $seq{$scaffold}){
        if ($scaffold =~ /^Plasmid{0,}/){
            $arr[1] =~ s/Plasmid/p/;
            $astart = $arr[2];
            $aend = $arr[3];
        }elsif($scaffold =~ /^Chromosome{0,}/){
            $arr[1] =~ s/Chromosome/Chr/;
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
		if($arr[2]<$arr[3] and $arr[5] =~ /^-$/ or $arr[2]<$arr[3] and $arr[4] =~ /^-$/){
		print FFN "\>$arr[0] $astart $aend $arr[1] $arr[3] $arr[2]\n$seq\n";
		my $ss_sart=$arr[2];
		my $ss_end=$arr[3];
		my @aa=@arr;
		$aa[3]=$ss_sart;
		$aa[2]=$ss_end;
		print LIST (join "\t",@aa),"\t", $astart,"\t", $aend,"\t", $initiator,"\t", $terminator, "\n";
		}else{
		print FFN "\>$arr[0] $astart $aend $arr[1] $arr[2] $arr[3]\n$seq\n";
		print LIST (join "\t",@arr),"\t", $astart,"\t", $aend,"\t", $initiator,"\t", $terminator, "\n";
		}
 	}
}
close PRE;

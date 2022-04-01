#!/usr/bin/perl -w
use strict;
use warnings;

if(@ARGV != 5){
	print STDERR "glimmer2ffn.pl <fasta> <glimmerPredict> <orf_prefix> <N_num> > <prefix>\n";
	print STDERR "This script is designed to extract ffnFile  for draft genome used glimmer software for it\'s orf predicting\n";
	print STDERR "The third parameter N_num is the number of N setting between contigs\n";
	exit;
}

my $fas = shift;
my $glimmerPredict = shift;
my $orf_prefix = shift;
my $N_num = shift;
my $prefix = shift;

(-s $fas) || die "$fas is not exists!\n";
open (FAS, "<$fas") or die $!;

my %id;
my %fas;
my $name = "";
while(<FAS>) {
	chomp($_);
	if($_ =~ /^>/) {
		$name = $_;
		$_ =~ /^>(\S+)/;
		$id{$1} = $name;
		$fas{$name} = "";
	}else{
		$_ =~ s/[^a-zA-z]//g;
		if($name ne "") {
			$fas{$name} .= $_;
		}
	}
}
close(FAS);
my %rev_complement;
foreach my $name (keys %fas) {
	$rev_complement{$name} = reverse $fas{$name};
	$rev_complement{$name} =~ tr/agctAGCT/TCGATCGA/;
}

my ($orf_id, $sequence_id,$orf) ;
my $contig = "";
my $accumulate = 0;
my $index = 0;
my $index2 = 0;
open (PREDICT, "<$glimmerPredict") or die;
open (FFN,">$prefix".".fna") or die;
open(LIST,">$prefix".".predict.gff") or die;
print LIST "Gene id\tSequence id\tStart\tEnd\tStrand\tGene Length(bp)\tProtein Length\n";
while(<PREDICT>) {
	chomp $_;
	if($_ =~ /\>(\S+)/){
		if(exists $id{$contig}){
			$accumulate = $accumulate + length($fas{$id{$contig}}) + $N_num;
		}else{
			$accumulate = $accumulate + 0 + $N_num;
		}
		$contig = $1;
		$index = 0;
	}else{
		$index ++;
		$index2 ++;
		my @line = split /\s+/, $_;
        my @str = split(/[0-9]/,$line[3],2);
        my $strand = $str[0];
		my $orf_num = ("0" x (4-length($index))).$index;
		my $orf_num2 = ("0" x (4-length($index2))).$index2;
        my $scaffold = $contig;
        my @scaffold_name = split /_/,$contig;
        my $sample = join "_", @scaffold_name;
        if (length($orf_prefix) != 0) {
        	$orf_id = $orf_prefix.$orf_num2;
        	$sequence_id = $scaffold."_ORF".$orf_num;
        }else{
        	$orf_id = $sample."_ORF".$orf_num2;
        	$sequence_id =  $scaffold."_ORF".$orf_num;
        }
        my $gene_length = abs($line[2]-$line[1]) + 1;
        my $protein_length = $gene_length/3 - 1;
		print LIST "$orf_id\t$sequence_id\t$line[1]\t$line[2]\t$strand\t$gene_length\t$protein_length\n";
		my $start = $accumulate + $line[1];
		my $end = $accumulate + $line[2];
		my $orf = "ORF".$orf_num;
		#if (length($orf_prefix) != 0) {
			#$orf = $orf_id;
		#}else{
			#$orf = "ORF".$orf_num;
		#}
		if($line[1]<$line[2]){
			my $seq = substr($fas{$id{$contig}},$line[1]-1,$line[2]-$line[1]+1);
			print FFN "\>$orf $start $end $contig $line[1] $line[2]\n$seq\n";
		}else{
			my $seq = substr($rev_complement{$id{$contig}},(length($fas{$id{$contig}})-$line[1]),$line[1]-$line[2]+1);
			print FFN "\>$orf $start $end $contig $line[1] $line[2]\n$seq\n";	
		}
 	}
}
close(PREDICT);

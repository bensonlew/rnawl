#!/usr/bin/perl -w
# __author__: qingchen.zhang @20200303

use warnings;
use strict;
use Bio::SeqIO;

die "usage: perl transgenescan_stat.pl ffn summary\n" unless(@ARGV == 2);

my ($input,$prefix)=@ARGV;
my $summary=$prefix.".summary";
fa_stat($input, $summary);

sub fa_stat{
	my $min=0;
	my $max=0;
	my $gene_num=0;
	my $GC_content=0;
	my $GC_number=0;
	my $Total_gene_length=0;
	my $avage_gene_length=0;
	my @x;
	my ($input,$output) =@_;
	my $in = Bio::SeqIO->new(-file => $input, -format => "fasta");
	while(my $seq =$in->next_seq){
    	push @x,$seq->length;
    	$gene_num++;
    	$GC_number += $seq->seq =~ tr/cgCG/cgCG/;
    	$Total_gene_length +=$seq->seq =~ tr/atcgATCG/atcgATCG/;
    	}
	@x =sort {$a<=>$b} @x;
	$min=$x[0];
	$max=$x[-1];
	$gene_num=@x;
	$GC_content =$GC_number/$Total_gene_length;
	$avage_gene_length=$Total_gene_length/$gene_num;
	open (OUT,">$output") || die $!;
	print OUT "gene_count\ttotal_length\taverage_length\tmax_length\tmin_length\tGC_content\n";
	print OUT "$gene_num\t$Total_gene_length\t$avage_gene_length\t$max\t$min\t$GC_content\n";
}

#system("rm -f $prefix.gene.as $prefix.gene.sn $prefix.gene.out $prefix.gene.gff");
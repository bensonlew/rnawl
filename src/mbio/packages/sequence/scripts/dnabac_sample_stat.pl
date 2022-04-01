#!/usr/bin/perl -w
use strict;
use warnings;

(@ARGV==4) || die "Usage: perl $0 [fasta] [fnn] [sample] [outdir]\nDiscription: Statistical gene predict information\n";

my($fna,$fnn,$sample,$outdir)=@ARGV;

my $contig_num=0;my $leng_string=0;my ($num_g,$num_c,$num_a,$num_t,$num_n)=(0)x5;
my ($location,%num_g,%num_c,%num_a,%num_t,%num_n);
(-s $fna) || die "$fna is not exists!\n";
(-s $fnn) || die "$fnn is not exists!\n";
open(FNA,"<$fna") or die $!;
while(<FNA>){
	chomp $_;
	if($_ =~ /^\>(\S+)/){
		$contig_num++;
		$location = $1;
		if ($location !~ /[sS]caffold[0-9].*/ and not exists $num_a{$location}){
			$num_a{$location} += (0)x5;
			$num_t{$location} += (0)x5;
			$num_g{$location} += (0)x5;
			$num_c{$location} += (0)x5;
			$num_n{$location} += (0)x5;
		}
	}else{
		$leng_string += length($_);
		$num_g += tr/gG/gG/;
		$num_c += tr/Cc/Cc/;
		$num_a += tr/Aa/Aa/;
		$num_t += tr/Tt/Tt/;
		$num_n += tr/Nn/Nn/;
		if ($location !~ /[sS]caffold[0-9].*/){
			$num_a{$location} += tr/Aa/Aa/;
			$num_t{$location} += tr/Tt/Tt/;
			$num_g{$location} += tr/gG/gG/;
			$num_c{$location} += tr/Cc/Cc/;
			$num_n{$location} += tr/Nn/Nn/;
		}
	}
}
close(FNA);
my $total_len = $num_a+$num_t+$num_g+$num_c+$num_n;
if($leng_string != $total_len){
       die "The sequence file contains orther chars expect AGCTN\n";
}

my %orf;my $name; my $gene_num=0;my($gene_a,$gene_t,$gene_g,$gene_c,$gene_n)=(0)x5;
my(%gene_a,%gene_t,%gene_g,%gene_c,%gene_n,%gene_num);
open(FFN,"<$fnn") or die $!;
while(<FFN>){
	chomp;
	if(/^>/){
		$_ =~ /^>(\S+).*\s(\S+)_.*/;
		$name = $1;
		$location = $2;
		if ($location !~ /[sS]caffold[0-9].*/ and not exists $gene_a{$location}){
			$gene_a{$location} += (0)x5;
			$gene_t{$location} += (0)x5;
			$gene_g{$location} += (0)x5;
			$gene_c{$location} += (0)x5;
			$gene_n{$location} += (0)x5;
			$gene_num{$location} = 0;
		}
		$orf{$name} = "";
		$gene_num++;
	}else{
		$_ =~ s/[^a-zA-z]//g;
		if($name ne ""){
			$orf{$name} .= $_;
		}
		$gene_a += tr/Aa/Aa/;
		$gene_t += tr/Tt/Tt/;
		$gene_g += tr/gG/gG/;
		$gene_c += tr/Cc/Cc/;
		$gene_n += tr/Nn/Nn/;
		if ($location !~ /[sS]caffold[0-9].*/){
			$gene_a{$location} += tr/Aa/Aa/;
			$gene_t{$location} += tr/Tt/Tt/;
			$gene_g{$location} += tr/gG/gG/;
			$gene_c{$location} += tr/Cc/Cc/;
			$gene_n{$location} += tr/Nn/Nn/;
			$gene_num{$location}++;
		}
	}
}

if (%gene_a){
	for my $k1 (keys %gene_a){
		my $loc = $k1;
		open OUT,">$outdir"."/$sample"."_$loc"."_CDS_statistics.xls" or die $!;
		print OUT "Sample\tGene num\tGene total length(bp)\tGene average length(bp)\tGene density(number/kb)\tGC content in gene region(%)\tGene/Genome(%)\tIntergenetic region length(bp)\tGC content in intergenetic region(%)\tIntergenetic length/Genome(%)\n";
		my $gene_total = $gene_g{$k1} +$gene_c{$k1}+$gene_a{$k1}+$gene_t{$k1}+$gene_n{$k1};
		my $gene_percent = sprintf("%.2f",($gene_total/$total_len)*100);
		my $inter_percent = sprintf("%.2f",($total_len-$gene_total)*100/$total_len);
		my $gene_average = sprintf("%.2f",$gene_total/$gene_num{$k1});
		my $inter_length = $total_len - $gene_total;
		my $gc_gene = ($gene_g{$k1}+$gene_c{$k1})/($gene_a{$k1}+$gene_t{$k1}+$gene_g{$k1}+$gene_c{$k1})*100;
		my $gc_gene_percent= sprintf("%.2f",$gc_gene);
		my $gc_inter = ($num_g+$num_c-$gene_g{$k1}-$gene_c{$k1})/($num_a+$num_t+$num_g+$num_c-$gene_g{$k1} -$gene_c{$k1}-$gene_a{$k1}-$gene_t{$k1})*100;
		my $gc_inter_percent= sprintf("%.2f",$gc_inter);
		my $CDS_density = sprintf("%.2f",($gene_num{$k1}*1000/$total_len));
		print OUT "$loc\t$gene_num{$k1}\t$gene_total\t$gene_average\t$CDS_density\t$gc_gene_percent\t$gene_percent\t$inter_length\t$gc_inter_percent\t$inter_percent\n";
		close OUT;
	}
}
my $gene_total = $gene_g +$gene_c+$gene_a+$gene_t+$gene_n;
my $gene_percent = sprintf("%.2f",($gene_total/$total_len)*100);
my $inter_percent = sprintf("%.2f",($total_len-$gene_total)*100/$total_len);
my $gene_average = sprintf("%.2f",$gene_total/$gene_num);
my $inter_length = $total_len - $gene_total;
my $gc_gene = ($gene_g+$gene_c)/($gene_a+$gene_t+$gene_g+$gene_c)*100;
my $gc_gene_percent= sprintf("%.2f",$gc_gene);
my $gc_inter = ($num_g+$num_c-$gene_g-$gene_c)/($num_a+$num_t+$num_g+$num_c-$gene_g -$gene_c-$gene_a-$gene_t)*100;
my $gc_inter_percent= sprintf("%.2f",$gc_inter);
my $CDS_density = sprintf("%.2f",($gene_num*1000/$total_len));

if (%gene_a){
	open OUT,">$outdir"."/$sample"."_whole_genome_CDS_statistics.xls"  || die $!;
	}else{
	open OUT,">$outdir"."/$sample"."_CDS_statistics.xls"  || die $!;
	}
print OUT "Sample\tGene num\tGene total length(bp)\tGene average length(bp)\tGene density(number/kb)\tGC content in gene region(%)\tGene/Genome(%)\tIntergenetic region length(bp)\tGC content in intergenetic region(%)\tIntergenetic length/Genome(%)\n";
my @path;
@path = split /\//,$fnn;
my $sample_name= $sample;
print OUT "$sample_name\t$gene_num\t$gene_total\t$gene_average\t$CDS_density\t$gc_gene_percent\t$gene_percent\t$inter_length\t$gc_inter_percent\t$inter_percent\n";
close FNA;
close OUT;


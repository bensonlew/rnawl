#!/usr/bin/perl -w
#�˽ű��ô�
#��������ĺ���͵����ļ��������ƻ���һ�£���Ҫ�������ļ��õ��ĵ��׽���ͺ��������Ʋ���׼��Ҫ����У��
use strict;
use warnings;

(@ARGV==5) || die "Usage: perl $0 [fasta] [fnn] [faa] [prefix] [outdir]\nDiscription: Statistical gene predict information\n";

my($fna,$fnn,$faa,$prefix,$outdir)=@ARGV;

my $contig_num=0;my $leng_string=0;my ($num_g,$num_c,$num_a,$num_t,$num_n)=(0)x5;
(-s $fna) || die "$fna is not exists!\n";
(-s $fnn) || die "$fnn is not exists!\n";
(-s $faa) || die "$faa is not exists!\n";
open(FNA,"<$fna") or die $!;
while(<FNA>){
	chomp $_;
	if($_ =~ /^\>(\S+)/){
		$contig_num++;
	}else{
		$leng_string += length($_);
		$num_g += tr/gG/gG/;
		$num_c += tr/Cc/Cc/;
		$num_a += tr/Aa/Aa/;
		$num_t += tr/Tt/Tt/;
		$num_n += tr/Nn/Nn/;
	}
}
close(FNA);
my $total_len = $num_a+$num_t+$num_g+$num_c+$num_n;
#if($leng_string != $total_len){
#       die "The sequence file contains orther chars expect AGCTN\n";
#}

my %orf;my $name;my $gene_num=0;my($gene_a,$gene_t,$gene_g,$gene_c,$gene_n)=(0)x5;
open(FFN,"<$fnn") or die $!;
while(<FFN>){
	chomp;
	if(/^>/){
		$_ =~ /^>(\S+)/;
		$name = $1;
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

system ("mkdir $outdir ") unless (-d $outdir);
open OUT,">>$outdir/gene_statistics.xls" || die $!;
print OUT "Sample\tGene num\tGene total length(bp)\tGene average length(bp)\tGene density(number/kb)\tGC content in gene region(%)\tGene/Genome(%)\tIntergenetic region length(bp)\tGC content in intergenetic region(%)\tIntergenetic length/Genome(%)\n";
my @path;
@path = split /\//,$fnn;
my $sample_name=$path[-1];
$sample_name=~s/\.fnn//g if($sample_name=~/\.fnn/);
$sample_name=~s/\.fnn//g if($sample_name=~/\.fnn/);
print OUT "$sample_name\t$gene_num\t$gene_total\t$gene_average\t$CDS_density\t$gc_gene_percent\t$gene_percent\t$inter_length\t$gc_inter_percent\t$inter_percent\n";
close FNA;
close OUT;

open(IN,"<$faa") or die $!;
my @file;
@file= split /\//,$faa;
open(FAAOUT,">$outdir/$file[-1]")||die;
$/=">";
<IN>;
while(my $block=<IN>){
	chomp $block;
    if ($block =~ /^gene_([0-9]*)\S+\|([0-9]*)\|([0-9]*)\s/){
        my $sec = <IN>;
        chomp $sec;
        my @all_lines = split /\n/,$sec;
        my @scaffold_name = split /_/,$all_lines[0];
        my $scaffold = pop(@scaffold_name);
        #my $scaffold = $all_lines[0];
        my $sample = join "_", @scaffold_name;
        my $orf_num = ("0"x(4-length($1))).$1;
        if (length($prefix) != 0) {
            $all_lines[0] = $prefix.$orf_num." ".$2." ".$3." ".$scaffold."_".$prefix.$orf_num." ".$2." ".$3;
        }else{
            $all_lines[0] = $sample."_ORF".$orf_num." ".$2." ".$3." ".$scaffold."_ORF".$orf_num." ".$2." ".$3;
        }
        my $str=substr($all_lines[1],1,length($all_lines[1])-1);
        $all_lines[1] ="M$str";
        print FAAOUT ">",(join "\n",@all_lines),"\n";
    }elsif ($block =~ /^ORF([0-9]*)_*1*\s+([0-9]*)\s+([0-9]*)\s+(\S+)\s+([0-9]*)\s+([0-9]*)/){
        my @all_lines = split /\n/,$block;
        my $orf_num = $1;
        my @scaffold_name = split /_/,$4;
        #my $scaffold = pop(@scaffold_name);
        my $scaffold = $4;
        my $sample = join "_", @scaffold_name;
        if (length($prefix) != 0) {
            $all_lines[0] = $prefix.$orf_num." ".$2." ".$3." ".$scaffold."_".$prefix.$orf_num." ".$5." ".$6;
        }else{
            $all_lines[0] = $sample."_ORF".$orf_num." ".$3." ".$4." ".$scaffold."_ORF".$orf_num." ".$6." ".$7;
        }
        my $str=substr($all_lines[1],1,length($all_lines[1])-1);
        $all_lines[1] ="M$str";
        print FAAOUT ">",(join "\n",@all_lines),"\n";
    }
}
close IN;
close FAAOUT;

open(FNN,"<$fnn") or die $!;
my @file2;
@file2 = split /\//,$fnn;
open(FNNOUT,">$outdir/$file2[-1]")||die;
$/=">";
<FNN>;
while(my $block=<FNN>){
	chomp $block;
    if ($block =~ /^gene_([0-9]*)\S+\|([0-9]*)\|([0-9]*)\s/){
        my $sec = <FNN>;
        chomp $sec;
        my @all_lines=split /\n/,$sec;
        my @scaffold_name = split /_/,$all_lines[0];
        my $scaffold = pop(@scaffold_name);
        my $sample = join "_", @scaffold_name;
        my $orf_num = ("0"x(4-length($1))).$1;
        if (length($prefix) != 0) {
            $all_lines[0] = $prefix.$orf_num." ".$2." ".$3." ".$scaffold."_".$prefix.$orf_num." ".$2." ".$3;
        }else{
            $all_lines[0] = $sample."_ORF".$orf_num." ".$2." ".$3." ".$scaffold."_ORF".$orf_num." ".$2." ".$3;
        }
        print FNNOUT ">",(join "\n",@all_lines),"\n";
    }elsif ($block =~ /^ORF([0-9]*)\s+([0-9]*)\s+([0-9]*)\s+(\S+)\s+([0-9]*)\s+([0-9]*)/){
        my @all_lines = split /\n/,$block;
        my $orf_num = $1;
        my @scaffold_name = split /_/,$4;
        #my $scaffold = pop(@scaffold_name);
        my $scaffold = $4;
        my $sample = join "_", @scaffold_name;
        if (length($prefix) != 0) {
            $all_lines[0] = $prefix.$orf_num." ".$2." ".$3." ".$scaffold."_".$prefix.$orf_num." ".$6." ".$7;
        }else{
            $all_lines[0] = $sample."_ORF".$orf_num." ".$2." ".$3." ".$scaffold."_ORF".$orf_num." ".$6." ".$7;
        }
    print FNNOUT ">",(join "\n",@all_lines),"\n";
    }
}
 close FNN;
 close FNNOUT;

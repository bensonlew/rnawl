#!/usr/bin/perl
use strict;
use warnings;
 
if(@ARGV!=3){
	print STDERR "perl orthomcl_stat.pl mclFile faa.list outDir\n";
	exit;
}

my $mclFile = shift;
my $faalist= shift;
my $outDir= shift;

if (!-d $outDir){
    system ("mkdir -p $outDir");
}

(-s $faalist) || die "$faalist is not exists!\n";
open LI,"<$faalist" || die $!;
my (@spec_names,%choose, %total_seq);
while(<LI>){ ###获取所有样本的gene_id
	chomp;
	my @l=split /\t/;
	push(@spec_names,$l[0]);
	(-s $l[1]) || die "$l[1] is not exists!\n";
	open FA,"<$l[1]" || die $!;
	my $id;
	while(<FA>){
		chomp;
		if(/^>/){
			$id = (split /\s/,$_)[0];
			$id =~ s/>//g;
			$choose{$l[0]}{$id}=0;
		}
	}
	close FA;
}
close LI;

my $sample_len = @spec_names;
#print "$sample_len\n";

my (%family_number,%gene_number);
open OUT1,"> $outDir/homologues_cluster.xls" || die $!;
print OUT1 "#ClusterID\tSample_number\tGene_number\t",join("\t",@spec_names),"\n";
open(ORTHOMCL,"<$mclFile") or die;
my $cluster_num=0;
while(<ORTHOMCL>){
	my %species_orf;
	for(my $j=0;$j<@spec_names;$j++){
		$species_orf{$spec_names[$j]} = "";
	}
	chomp $_;
	my @line = split /\t/, $_;
	my @clsuter = split(/\(/,$line[0]);
	$cluster_num ++;
	my $newcluster = "CLUSTER".$cluster_num ;
	print OUT1 "$newcluster";
	my @line2 = split /\s+/,$line[1];
	for(my $i=0;$i<@line2;$i++){	
		if($line2[$i] =~ /(\S+)\((\S+)\)/){
			#my @gene_name = split(/\|/,$1);
			#$choose{$2}{$gene_name[1]} = 1;
			$species_orf{$2} .= $1.",";
		}
	}
	my (@tmp_gene_num);
	my $tmp_len =0;
	my $tmp_sample_len=0;
	my (@total_sample,$sample_num );
	for (my $m=0;$m<@spec_names;$m++){
		if ($species_orf{$spec_names[$m]} ne ""){
			my $sample_genes = $species_orf{$spec_names[$m]};
			$sample_genes =~ s/;$//;
			$sample_genes =~ s/,$//;
			push (@total_sample, $sample_genes);
			@tmp_gene_num = split(/,/,$sample_genes);
			$sample_num = @tmp_gene_num;
			if ($sample_num != 0){
				$tmp_len += @tmp_gene_num;
				$tmp_sample_len += 1;
			}
			#push (@sample_exists,$spec_names[$m]);
		}else{
			push (@total_sample, "-");
		}
	}
	my $sample_exists = join("\t",@total_sample);
	print OUT1 "\t$tmp_sample_len\t$tmp_len\t$sample_exists\n";
}

=pod
{
	my $sample_exists = join(";",@sample_exists);
	print OUT1 "\t$sample_exists\n";
	if (!exists $family_number{$sample_exists}){
		$gene_number{$sample_exists} = $tmp_len;
		$family_number{$sample_exists} = 1;
	}else{
		$gene_number{$sample_exists} += $tmp_len;
		$family_number{$sample_exists} += 1;
	}
}

for(my $i=0;$i<@spec_names;$i++){
	my $single_num = 1;
	for my $j(keys %{$choose{$spec_names[$i]}}){
		if($choose{$spec_names[$i]}{$j} == 0){
			#print OUT1 "SINGLE"."-"x($i+1)."\t$j\n";
			print OUT1 "SINGLE".$single_num."\t-"x$i."\t"."$spec_names[$i]|".$j."\t-"x($sample_len - $i-1)."\t$spec_names[$i]\n";
			$single_num += 1;
			if (!exists $family_number{$spec_names[$i]}){
				$family_number{$spec_names[$i]} = 1;
				$gene_number{$spec_names[$i]} = 1;
			}
			else{
				$family_number{$spec_names[$i]} += 1;
				$gene_number{$spec_names[$i]} += 1;
			}
		}
	}
}

=cut
close(ORTHOMCL);
close OUT1;



#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
my $in;
my $snp_list="/mnt/ilustre/users/changliang.qiu/danni.li/source/SOGG_snp.xls";
#my $tree="/mnt/clustre/users/danni.li/tree/programe/20170728.txt";
#my $tree="/mnt/ilustre/users/changliang.qiu/danni.li/source/20181027_sourse_tree.txt";//0310
#my $tree="/mnt/ilustre/users/changliang.qiu/danni.li/source/20200127tree.txt";//0505
# my $tree="/mnt/ilustre/users/changliang.qiu/danni.li/source/20200505_tree.txt";//0606
# my $tree="/mnt/ilustre/users/changliang.qiu/danni.li/source/20200606_tree.txt";
# my $pos_snp = "/mnt/ilustre/users/changliang.qiu/danni.li/source/Y_20200606_list.txt";//0606fix
# my $tree="/mnt/ilustre/users/changliang.qiu/danni.li/source/20200801_tree.txt";//0801
# my $pos_snp = "/mnt/ilustre/users/changliang.qiu/danni.li/source/Y_20200801_list.txt";//0801
my $tree="/mnt/ilustre/users/changliang.qiu/danni.li/source/20200801_tree_fix.txt";
my $pos_snp = "/mnt/ilustre/users/changliang.qiu/danni.li/source/Y_20200801_list_fix.txt";
GetOptions(
'vcf=s'=>\$in,
'SNP=s'=>\$snp_list,
'tree=s'=>\$tree,
'pos=s' =>\$pos_snp,
);
my $usage=<<USAGE;
-vcf	<str> input Yfull vcf
-SNP	<str> eg:/mnt/clustre/users/danni.li/sourse/Y_full/test_result_full/programe/SOGG_snp.xls(SOGG SNP)
-tree   <str> tree file 
-pos	<str> pos_list
USAGE
if(!$in ||!$snp_list||!$tree){
	print "$usage\n";
	exit;
	}
my %pos_alt_SNP;
# get_SNP();###get all known SNP type###  #for get all pos to _YFULL_SNP.txt from vcf   del at 20200701
get_pos_snp();   ###get pos_snp from tree  add at 20200701
my %SNP_Y;
my %SNP_N;
my %SNP_ND;
my %sample_SNPpos;
my %pos_snp;
my %pos_ref;
open(IN,"$in");
print "input $in\n";
while(<IN>){
	chomp;
	next if($_=~/^#/);
	my @array=split"\t",$_;
	$array[7]=~/DP=(\d+);/;
#	next if($1<3);
	my $ALT;
	if($array[4] eq "."){
		$ALT=$array[3];
		}
	else{
		my @tt=split",",$array[4];
		$ALT=$tt[0];
		}
	$pos_alt_SNP{$array[1]}->{$ALT}=$pos_snp{$array[1]};
	$sample_SNPpos{$array[1]}=$ALT;
	if(defined $pos_ref{$array[1]}){
		my $pos = $array[1];
		if($pos_ref{$array[1]} eq lc($ALT)){
			$SNP_N{$array[1]}->{$ALT}++;
			}
		else{
			$SNP_Y{$array[1]}->{$ALT}++;
			}
		# if($pos_ref{$array[1]} =~ /$ALT/){
		# 	$SNP_Y{$array[1]}->{$ALT}++;
		# 	}
		# else{
		# 	$SNP_N{$array[1]}->{$ALT}++;
		# 	}
	}
	else{
		$SNP_ND{$array[1]}->{$ALT}++;
	}
}
close IN;

open(OUT,">YFULL_SNP_reslt.xls");
print OUT "pos\tSNP\tresult\talt\n";
foreach my $key (sort {$a<=>$b} keys %pos_alt_SNP){
		foreach my $key1 (keys %{$pos_alt_SNP{$key}}){
			if(defined $SNP_Y{$key}->{$key1}){
				print  OUT "$key\t$pos_alt_SNP{$key}->{$key1}\t+\t$sample_SNPpos{$key}\n";
				}
			elsif(defined $SNP_N{$key}){
				 print OUT "$key\t$pos_alt_SNP{$key}->{$key1}\t-\t$sample_SNPpos{$key}\n";
				}
			elsif(defined $SNP_ND{$key}){
				# print OUT "$key\t$pos_alt_SNP{$key}->{$key1}\tnot detected\t$sample_SNPpos{$key}\n";
				print OUT "$key\t$pos_alt_SNP{$key}->{$key1}\tnot detected\n";
				}
			else{
				print  OUT "$key\t$pos_alt_SNP{$key}->{$key1}\tnot detected\n";
				}
		}
	}

sub get_SNP{
	open(IN,"$snp_list");
	while(<IN>){
		chomp;
		my @array=split"\t",$_;
		$pos_alt_SNP{$array[0]}->{$array[2]}=$array[1];
		}
	close IN;
	
	}

sub get_pos_snp{
	open(IN,"$pos_snp");
	while(<IN>){
		chomp;
		my @array=split"\t",$_;
		my @ref_list;
		if(scalar(@array) > 5){
			$pos_ref{$array[0]}=lc($array[2]);
			# if($array[4] eq "*"){
			# 	my @ref = ("A", "G", "C", "T");
			# 	foreach my$rr (@ref){
			# 		if($rr ne $array[2]){
			# 			push(@ref_list, $rr);
			# 		}
			# 	}
			# 	$pos_ref{$array[0]}=join",",@ref_list;
			# }
			# else{
			# 	push(@ref_list, $array[2]);
			# 	$pos_ref{$array[0]}=join",",@ref_list;
			# }
		}
		else{
			$pos_ref{$array[0]}=lc($array[2]);
		}
		$pos_snp{$array[0]}=$array[1];
	}
	close IN;
	
}

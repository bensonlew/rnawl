#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
my $in;
my $snp_list;
#my $snp_list="/mnt/clustre/users/danni.li/sourse/Y_full/test_result_full/programe/8Msnp.xls";
# my $snp_list="/mnt/ilustre/users/changliang.qiu/danni.li/source/SOGG_snp.xls";
GetOptions(
'vcf=s'=>\$in,
'SNP=s'=>\$snp_list,
);
my $usage=<<USAGE;
-vcf	<str> input 8M vcf 
-SNP	<str> 8M SNP eg:/mnt/clustre/users/danni.li/sourse/Y_full/test_result_full/programe/8Msnp.xls(8M bed SNP)
USAGE
if(!$in||!$snp_list){
	print "$usage\n";
	exit;
	}

my %pos_alt_SNP;
get_SNP();###get 8M known SNP type###
open(IN,"$in");
print "input $in\n";
open(OUT1,">sample_SNP.txt");
open(OUT2,">sample_list.txt");
print OUT1 "pos\tref\talt\n";
while(<IN>){
	chomp;
	next if($_=~/^#/);
	my @array=split"\t",$_;
	$array[7]=~/DP=(\d+);/;
	next if($1<5);
	my $ALT;
	if($array[4] eq "."){
		$ALT=$array[3];
		}
	else{
		$ALT=$array[4];
		}
	if(defined $pos_alt_SNP{$array[1]}){
		}
	elsif($array[9]=~/1\/1/){
		print OUT1 "$array[1]\t$array[3]\t$array[4]\n";
		print OUT2 "$array[1]$array[3]>$array[4];";
		}
	}
close OUT1;
close OUT2;

close IN;
=pod
open(OUT,">YFULL_SNP_reslt.xls");
foreach my $key (sort {$a<=>$b} keys %pos_alt_SNP){
		foreach my $key1 (keys %{$pos_alt_SNP{$key}}){
			if(defined $SNP_Y{$key}->{$key1}){
				print  OUT "$key\t$pos_alt_SNP{$key}->{$key1}\t+\n";
				}
			elsif(defined $SNP_N{$key}){
				 print OUT "$key\t$pos_alt_SNP{$key}->{$key1}\t-\n";
				}
			else{
				print  OUT "$key\t$pos_alt_SNP{$key}->{$key1}\tnot been detected\n";
				}
			}
		}
=cut
sub get_SNP{
	open(IN,"$snp_list");
	while(<IN>){
		chomp;
		my @array=split"\t",$_;
		$pos_alt_SNP{$array[0]}->{$array[2]}=$array[1];
		}
	close IN;
	}

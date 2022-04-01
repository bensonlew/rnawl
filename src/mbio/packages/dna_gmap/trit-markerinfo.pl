#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($trit,$vcf,$out,$fid,$mid,$pop);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"trit:s"=>\$trit,
	"vcf:s"=>\$vcf,
	"fid:s"=>\$fid,
	"mid:s"=>\$mid,
	"out:s"=>\$out,
	"pop:s"=>\$pop,
			) or &USAGE;
&USAGE unless ($vcf and $trit and $out and $fid and $mid);

$pop||="F2";
my %stat;
my $trait;
open IN,$trit;
my(@sample,@value);
if($pop ne "CP"){
	while(<IN>){
		chomp;
		next if($_ eq ""|| /^$/);
		my($ni,@indi)=split(/\,/,$_);
		if(/^Genotype/){
			push @sample,@indi;
			next;
		}else{
			$trait=$ni;
			push @value,@indi;
		}
		print join("\t",@sample),"\n";
	}
	for(my$i=0;$i<scalar@sample;$i++){
		my $sample=$sample[$i];
		my $value=$value[$i];
		$stat{$sample}=$value;
	}
}else{
	while(<IN>){
		chomp;
		next if($_ eq ""|| /^ntrt/|| /^nind/|| /^miss/);
		if(/^sampleID/){
			$trait=(split(/\s+/,$_))[1];
			next;
		}else{
			my($sample,$value)=split(/\s+/,$_);
			$stat{$sample}=$value;
		}
	}
}
close IN;
open OUT,">$out/$trait.assocation.xls";
# print OUT "#Genotype  Frequency Distribution\tTrait\tChr\tPos\tRef\tAlt\tType\tF1 Genotype\tF1 Depth\tM1 Genotype\tM1 Depth\n";
print OUT "Trait\tChr\tPos\tRef\tAlt\tType\tF1 Genotype\tF1 Depth\tM1 Genotype\tM1 Depth\n";
my @Indi;
my %region;
open IN,$vcf;
while(<IN>){
	chomp;
	next if($_ eq ""|| /^$/|| /^#@/|| /^@/);
	my($chr,$pos,$id,$ref,$alt,$qual,$filt,$info,$format,@undi)=split(/\t/,$_);
	if(/^#CHROM/){
		push @Indi,@undi;
	}else{
		my $defp=join(",",$ref,$alt);
		my @alt=split(/,/,$defp);
		my($type,$geno,$d,$fgeno,$fd,$mgeno,$md,$eff);
		if($alt[0] eq $alt[1]){
			$type="SNP";
		}else{
			$type="InDel";
		}
		my @format=split(/:/,$format);
		if($info=~/ANN=([^\;]*)/){
			$eff=(split(/\|/,$1))[2];
		}
		my($gt,$ad);
		for (my $i=0;$i<@format;$i++) {
			$gt=$i if ($format[$i] eq "GT");
			$ad=$i if ($format[$i] eq "AD");
		}
		open Mark,">$out/$trait-$id-info.xls";
		print Mark "#SampleID\tGenotype\tTrait\tTrait Value\n";
		for(my$i=0;$i<@Indi;$i++){
			my $sam=$Indi[$i];
			my ($g1,$g2)=split(/\//,(split(/\:/,$undi[$i]))[$gt]);
			my ($d1,$d2)=split(/\,/,(split(/\:/,$undi[$i]))[$ad]);
			if($g1 eq $g2){
				if($g1 eq "\."){
					$geno="--";
					$d="--";
				}else{
					$geno=$alt[$g1];
					$d=$d1 + $d2;
				}
			}else{
				$geno=join(",",$alt[$g1],$alt[$g2]);
				$d=join(",",$d1,$d2);
			}
			if($sam eq $fid){
				$fgeno=$geno;
				$fd=$d;
			}elsif($sam eq $mid){
				$mgeno=$geno;
				$md=$d;
			}
			foreach my$sample(sort keys %stat){
				if($sample eq $sam){
					print Mark "$sample\t$geno\t$trait\t$stat{$sample}\n";
				}
			}
		}
		close Mark;
		# print OUT "$out/$trait.$id.info.xls\t$trait\t$chr\t$pos\t$ref\t$alt\t$type\t$eff\t$fgeno\t$fd\t$mgeno\t$md\n";
		# print OUT "$trait\t$chr\t$pos\t$ref\t$alt\t$type\t$eff\t$fgeno\t$fd\t$mgeno\t$md\n";
		print OUT "$trait\t$chr\t$pos\t$ref\t$alt\t$type\t$fgeno\t$fd\t$mgeno\t$md\n";
	}
}
close IN;



#######################################################################################
#print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub USAGE {#
        my $usage=<<"USAGE";
Contact:        caixia.tian\@majorbio.com;
Script:			$Script
Description:
	fq thanslate to fa format
	eg:
	perl $Script -i -o -k -c

Usage:
  Options:
  -vcf	<file>	trait*.vcf.total
  -trit	<file>	Trait file
  -fid	<str>	fid
  -mid	<str>	mid
  -out	<dir>	output result dir
  -pop	<str>	pop type|F2 or CP..
  -h         Help

USAGE
        print $usage;
        exit;
}

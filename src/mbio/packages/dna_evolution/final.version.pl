#!/mnt/bin/env perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($ref,$gff,$vcf,$list,$gene,$out);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
        "help|?" =>\&USAGE,
        "ref:s"=>\$ref,
        "gff:s"=>\$gff,
        "vcf:s"=>\$vcf,
		"out:s"=>\$out,
		"list:s"=>\$list,
		"gene:s"=>\$gene,
                        ) or &USAGE;
&USAGE unless ($ref and $vcf and $list);
my %comb=(
	"10"=>"GGCTAC",
	"A\/G"=>"R",
	"G\/A"=>"R",
	"A\/C"=>"M",
	"C\/A"=>"M",
	"A\/T"=>"W",
	"T\/A"=>"W",
	"C\/T"=>"Y",
	"T\/C"=>"Y",
	"G\/T"=>"K",
	"T\/G"=>"K",
	"G\/C"=>"S",
	"C\/G"=>"S",
);
my ($genename,$chr,$pos1,$pos2);
if(! $gff){
	open Gene,$gene;
	while (<Gene>){
		chomp;
		next if ($_ eq "" || /^$/);
		($genename,$chr,$pos1,$pos2)=split(/\t/,$_);
	}
	close Gene;
}else{
	open Gff,$gff;
	if($gff=~/gz$/){
		close Gff;
		open Gff,"gunzip -c $gff|";
	}
	$/="\n";
	while (<Gff>){
		chomp;
		next if ($_ eq "" || /^$/|| /^#/);#chr1	RefSeq	gene	3631	5899	.	+	.	ID=gene0;Dbxref=Araport:AT1G01010,TAIR:AT1G01010,GeneID:839580;Name=NAC001;gbkey=Gene;gene=NAC001;gene_biotype=protein_coding;gene_synonym=ANAC001,NAC	domain	containing	protein	1,T25K16.1,T25K16_1;locus_tag=AT1G01010
		my ($chrid,$seqty,$type,$start,$end,undef,undef,undef,$info)=split(/\t/,$_);
		my $geneid;
		# next if($type ne "gene");  # 传进来的gene不一定是gene，可能是rna，因此注释
		# if($info=~/ID=(\w+)\;(.*)/){  # 若是有中划线的时候-匹配不到，因此进行修改
    if($info=~/ID=(.*?)\;(.*)/){
			$geneid=$1;
			if($gene eq $geneid){
				$genename=$gene;
				$chr=$chrid;
				$pos1=$start;
				$pos2=$end;
			}else{
				next;
			}
		}
	}
	close Gff;
}
my $geneseq;
$/=">";
open Ref,$ref;
if($ref=~/gz$/){
	close Ref;
	open Ref,"gunzip -c $ref|";
}
while(<Ref>){
	chomp;
	next if ($_ eq ""|| /^$/);
	my($id,@seq)=split(/\s+/,$_);
	next if($id ne $chr);
	my$seq=join("",@seq);
	my $length=$pos2 - $pos1 + 1;
	my$start=$pos1 - 1;
	$geneseq=substr($seq,$pos1,$length);
	}
close Ref;
$/="\n";
my @samples;
open List,$list;
while(<List>){
	chomp;
	next if($_ eq "" || /^$/);
	push @samples,$_;
}
close List;

$/="\n";
open In,$vcf;
if($vcf=~/gz$/){
	close In;
	open In,"gunzip -c $vcf|";
}
my @sample;
my %stat;
while (<In>){
	chomp;
	next if ($_ eq  ""|| /^$/|| /^##/);
	#my($vchr,$pos,$id,$ref,$alt,$qual,$filter,$info,$format,@undi)=split(/\t/,$_);#chr1	1072	chr1_1072	A	C	154.84	PASS	AC=2;AF=1.00;AN=2;ANN=C|upstream_gene_variant|MODIFIER|NAC001|gene0|transcript|rna0|protein_coding||c.-2688A>C|||||2559|,C|intergenic_region|MODIFIER|CHR_START-NAC001|CHR_START-gene0|intergenic_region|CHR_START-gene0|||n.1072A>C||||||;DP=6;ExcessHet=3.0103;FS=0;MLEAC=2;MLEAF=1;MQ=60;QD=25.81;SOR=2.303;set=variant	GT:AD:DP:GQ:PL	1/1:0,6:6:18:183,18,0
	if($_ =~ /^#/){
		my($vchr,$pos,$id,$ref,$alt,$qual,$filter,$info,$format,@undi)=split(/\t/,$_);
		@sample=@undi;
		#print join("\t",@sample),"\n";
		#next;
	}else{
		my($vchr,$pos,$id,$ref,$alt,$qual,$filter,$info,$format,@undi)=split(/\t/,$_);
		next if($vchr ne $chr);	#!=
		next if($pos < $pos1);
		next if($pos > $pos2);
		#if($pos>=$pos1 and $pos<= $geneend)
		my @alt=split(/,/,join(",",$ref,$alt));
		#print $vchr,"\t",$pos,join("\t",@undi),"\n";
		#die;

		my $newpos=$pos - $pos1;
		my @format=split(/:/,$format);

		for (my $i=0;$i<@sample;$i++) {
			my $sample=$sample[$i];
			foreach my$samples(@samples){
				next if($samples ne $sample);
				my @info=split(/:/,$undi[$i]);
				for (my $j=0;$j<@info;$j++) {
					if ($format[$j] eq "GT") {
						next if ($info[$j] eq "./." || $info[$j] eq "0/0");
						my ($g1,$g2)=split(/\//,$info[$j]);
						if ($g1 eq $g2) {
							my $newinfo=join("\t",$alt[$g1],$newpos);
							push @{$stat{$samples}},$newinfo;
						}else{
							my$var="$alt[$g1]\/$alt[$g2]";
							if(length$var eq "3"){
								my $newinfo=join("\t",$comb{$var},$newpos);
								push @{$stat{$samples}},$newinfo;
							}else{
								my $newinfo=join("\t",$alt[$g1],$newpos);
								push @{$stat{$samples}},$newinfo;
							}
						}
					}
				}
			}
		}
	}
}
close In;
open Out,">$out";
print Out ">$genename\t$chr\t$pos1\t$pos2\n$geneseq\n";
foreach my $samples (keys %stat){
	print Out ">$samples\n";
	my@vars=();
	my@segment=();
	my$last=0;
	foreach my$info(@{$stat{$samples}}){
		my($var,$pos)=split(/\t/,$info);
		my $len=$pos - $last;
		push @vars,$var;
		push @segment,substr($geneseq,$last,$len);
		$last=$pos;
	}
	for (my$i=0;$i<scalar@vars;$i++){
		print Out $segment[$i],$vars[$i];
	}
	print Out "\n";
	#print Out "$samples\n",join("\n",@{$stat{$samples}}),"\n";
}
close Out;

#######################################################################################
#print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
########################################################################################
sub USAGE {#
        my $usage=<<"USAGE";
Contact:        caixia.tian\@majorbio.com;
Script:                 $Script
Description:
        get gene's seq
        eg:
        perl $Script -ref -gff -vcf -out -gene -list

Usage:
  Options:
	-ref	<file>  input ref.fa
	-gff	<file>  input ref.gff
	-vcf	<file>	pop.final.vcf
	-out	<file>	output file name
	-gene	<str>	geneid
	-list	<file>	input sample.list
	-h			Help

USAGE
        print $usage;
        exit;
}

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
&USAGE unless ($ref and $gff and $vcf and $list);
my %ref;
open Ref,$ref;
if($ref=~/gz$/){
        close Ref;
        open Ref,"gunzip -c $ref|";
}
$/=">";
while(<Ref>){
        chomp;
        next if ($_ eq ""|| /^$/);
        my($chr,@seq)=split(/\s+/,$_);
        my $seq=join("",@seq);
        $ref{$chr}{seq}=$seq;
}
close Ref;
my $geneinfo;
my $refs;
open Gff,$gff;
if($gff=~/gz$/){
	close Gff;
	open Gff,"gunzip -c $gff|";
}
$/="\n";
while (<Gff>){
	chomp;
	next if ($_ eq "" || /^$/|| /^#/);#chr1	RefSeq	gene	3631	5899	.	+	.	ID=gene0;Dbxref=Araport:AT1G01010,TAIR:AT1G01010,GeneID:839580;Name=NAC001;gbkey=Gene;gene=NAC001;gene_biotype=protein_coding;gene_synonym=ANAC001,NAC	domain	containing	protein	1,T25K16.1,T25K16_1;locus_tag=AT1G01010
	my ($chr,$seqty,$type,$start,$end,undef,undef,undef,$info)=split(/\t/,$_);
	next if ($type eq "region");
	if ($type eq "gene" |$type eq "mRNA") {
		if ($info =~ /ID=([^;]*)/) {
			next if ($1 ne $gene);
			$geneinfo=join("_",$chr,$start,$end);
			$refs=substr($ref{$chr}{seq},$start,$end-$start+1);
		}
		last;
	}
}
close Gff;

my @samples;
open List,$list;
while(<List>){
	chomp;
	next if($_ eq "" || /^$/);
	push @samples,$_;
}
close List;
die "没有找到gene_id：$gene" if(!defined $geneinfo);
my %sample;
foreach my $sample (@samples) {
	@{$sample{$sample}}=split(//,$refs);
}
my %stat;
open In,$vcf;
if($vcf=~/gz$/){
	close In;
	open In,"gunzip -c $vcf|";
}
my @sample;
while (<In>){
	chomp;
	next if ($_ eq "" || /^$/|| /##/);
	my($chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$format,@undi)=split(/\t/,$_);#chr1	1072	chr1_1072	A	C	154.84	PASS	AC=2;AF=1.00;AN=2;ANN=C|upstream_gene_variant|MODIFIER|NAC001|gene0|transcript|rna0|protein_coding||c.-2688A>C|||||2559|,C|intergenic_region|MODIFIER|CHR_START-NAC001|CHR_START-gene0|intergenic_region|CHR_START-gene0|||n.1072A>C||||||;DP=6;ExcessHet=3.0103;FS=0;MLEAC=2;MLEAF=1;MQ=60;QD=25.81;SOR=2.303;set=variant	GT:AD:DP:GQ:PL	1/1:0,6:6:18:183,18,0
	if($chr eq "#CHROM"){
		@sample=@undi;
		#print join("\t",@sample),"\n";
		#die;
	}else{
    my ($genechr,$genestart,$geneend)=split(/\_/,$geneinfo);
		next if($chr ne $genechr);	#!=
		if($pos>=$genestart and $pos<= $geneend){
			my @alt=split(/,/,join(",",$ref,$alt));
			my (%ale,%len);
			for (my $i=0;$i<@alt;$i++) {
				$ale{$alt[$i]}=$i;
				$len{length($alt[$i])}=1;
			}
			my $snptype="SNP";
			$snptype="INDEL"  if (scalar keys %len > 1);

			my $snppos=$pos - $genestart;
			my @format=split(/:/,$format);
			for (my $i=0;$i<@sample;$i++) {
				my $sample=$sample[$i];
				my @info=split(/:/,$undi[$i]);
				for (my $j=0;$j<@info;$j++) {
					if ($format[$j] eq "GT") {
						#$stat{$sample}{$snppos}{type}=$ref if ($info[$j] eq "./." || $info[$j] eq "0/0");
						my ($g1,$g2)=split(/\//,$info[$j]);
						#print "$id\t$alt[$g1]\/$alt[$g2]\t$snppos\n";
						if ($g1 eq $g2) {
							$sample{$sample}[$snppos]=$alt[$g1];
						}else{
							my $type=length$ref ;
							if($snptype eq "SNP"){
								my $basetype=join("\/",$alt[$g1],$alt[$g2]);
								$sample{$sample}[$snppos]="R" if($basetype eq "A\/G");
								$sample{$sample}[$snppos]="M" if($basetype eq "A\/C");
								$sample{$sample}[$snppos]="W" if($basetype eq "A\/T");
								$sample{$sample}[$snppos]="Y" if($basetype eq "C\/T");
								$sample{$sample}[$snppos]="K" if($basetype eq "G\/T");
								$sample{$sample}[$snppos]="S" if($basetype eq "G\/C");
							}else{
								my $lang;
								if(length$alt[$g1] >= length$alt[$g2]){
									$lang = $alt[$g1];
								}else{
									$lang = $alt[$g2];
								}
								$sample{$sample}[$snppos]=$lang;
							}
						}
					}
				}
			}
		}
	}
}
close In;
open Out,">$out/gene.fa";
foreach my $sample (sort keys %sample) {
	print Out ">$sample\n";
	print Out join("",@{$sample{$sample}}),"\n";
}
close Out;
# my $job="clustalo -i $out/$gene.fa -o $gene.diff.phy --outfmt=phy ";
# `$job`;

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

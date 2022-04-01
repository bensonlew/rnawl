#! /usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
my %opts;
GetOptions (\%opts,"rfam=s","mrd=s","arf=s","max_map=i","min_count=i","om=s","oa=s");
my $usage = <<"USAGE";

        Program :    $0
        Discription: prepare map and fa files for mireap in plant sRNA project
        Usage:       perl $0 [options]
	Basic:
                -rfam    	file    ex: sRNA seq fa after rfam trimming, rfam_trimed.fa
                -mrd  	  	file    ex: mature miRNA align file from mirdeep, miRBase.mrd 
		-arf    	file    ex: genome alignment form mirdeep preprocess, reads_vs_genome.arf_parsed.arf
        Opintional:
		-max_map 	numic	ex: max mapped loc num, defult:20
		-min_count	numic	ex: min novel miRNA counts, defult:10
                -om      	file    ex: output aln file, defult: map.txt
        	-oa	 	file	ex: output fa file, defult: filtered.fa
        Usage: perl $0 -rfam <rfam_trimed.fa> -mrd <miRBase.mrd> -arf <reads_vs_genome.arf_parsed.arf> -max_map 20 -min_count 10 -om <map.txt> -oa <filtered.fa>
        Contect: yuntao.guo\@majorbio.com
        Data: 2016-04-28
USAGE
die $usage if (!($opts{rfam} && $opts{mrd} && $opts{arf}));
$opts{max_map}=$opts{max_map}? $opts{max_map}:"20";
$opts{min_count}=$opts{min_count}? $opts{min_count}:"10";
$opts{om}=$opts{om}? $opts{om}:"map.txt";
$opts{oa}=$opts{oa}? $opts{oa}:"filtered.fa";

my (%fa,%mrd,%arf,%exp,%list,%align,$fa_seq,$fa_id);
#open fa and store into %fa
open FA,"<$opts{rfam}" || die "can't open $opts{rfam}$!\n";
while (<FA>){
        chomp;
        if(/^>/){
                if ($fa_seq){
                        $fa{$fa_id}=$fa_seq;
                        $fa_seq='';
                }
                $fa_id=(split /\s+/,$_)[0];
        }else{
                $fa_seq.=$_;
        }
}
$fa{$fa_id}=$fa_seq;
close FA;
#open mrd file and store known miR id into %mrd
open MRD,"<$opts{mrd}" || die "$! in read $opts{mrd}\n";
while(<MRD>){
	chomp;
	if (/^(\S+_\S+_x\d+)\s+/){
		$mrd{$1}='';
	}
}
close MRD;
#open arf file and filter by mapped times and counts num
open ARF,"<$opts{arf}" || die "$! in read $opts{arf}\n";
while (<ARF>){
	chomp;
	my @line=split /\t/;
	my $ex=($line[0]=~/_x(\d+)$/)? $1:die "$. bad format!";
	$exp{$line[0]}=$ex;
	$align{$line[0]}+=1;
	$arf{$line[0]}=$line[5]."\t".$line[7]."\t".$line[8]."\t".$line[10];
}
close ARF;
#out map.txt
open MAP,">$opts{om}"|| die "can't open $opts{om} $!\n";
for my $nov (keys %arf){
	next if (($align{$nov} > $opts{max_map} ) or ($exp{$nov} < $opts{min_count}) or (exists $mrd{$nov}));
	$list{$nov}='';
	my $id_nov=$nov;
	$id_nov =~ s/_x\d+$//g;
	print MAP "$id_nov\t$arf{$nov}\n";
}
close MAP;
#out filted.fa
open FILTER,">$opts{oa}" || die "can't open $opts{oa} $!\n";
for my $nov (keys %list){
	my $id_nov=">".$nov;
	die ("not cantain $nov in $opts{rfam}\n") if (!exists $fa{$id_nov});
	my $ids=$id_nov;
	$ids =~ s/_x/ /g;
	print FILTER "$ids\n$fa{$id_nov}\n";
}
close FILTER;

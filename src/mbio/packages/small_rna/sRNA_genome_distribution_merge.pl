#!/usr/bin/perl -w
use strict;
use warnings;
use Carp;
use List::Util qw ( sum);
use Math::Round qw(:all);
use Getopt::Long;
use Config::IniFiles;

my %opts;
my $VERSION="2.0";
GetOptions( \%opts,"fa=s", "arf=s","o=s","ref_fai=s","config=s","bar!", "h!");

my $usage = <<"USAGE";
       Program : $0
       Version : $VERSION
       Contact : yuntao.guo\@majorbio.com
       Lastest modify:2018-04-27 liubinxu
       Discription: calculate seq distribution for each sample after blast to rfam and prepare files for mirdeep2 and mireap
       Usage :perl $0 [options]
                -fa*		rfam_trimed.fa		fasta file input for  mapper.pl
                -config* 	Uniq.config.ini		config file for uniqueFastaqs.pl         
                -arf*		reads_vs_genome.arf	read mapped to genome arf file create by mapper.pl from miRdeep2
                -o		map_stat.xls		out stat prefix,default: map_stat.xls
                -bar					bar plot
                -ref_fai
                -h					Display this usage information
                * 					must be given Argument      
USAGE

die $usage if ( !( $opts{fa} && $opts{config} && $opts{arf}) && $opts{ref_fai} || $opts{h} );
$opts{o} = $opts{o}? $opts{o} : "map_stat.xls";

my $cfg = Config::IniFiles->new(-file => $opts{config});
my @names=$cfg->Parameters("FASTA");
my %titles;
foreach my $f (@names){
	my $name;
	if($cfg->SectionExists("NAME")&& $cfg->val('NAME',$f)){
		$name=$cfg->val('NAME',$f);
		$titles{$f}=$name;
	}else{
		$name=$f;
		$titles{$f}=$f;
	}
}

	my %seqs_uniq;
	my %seqs_num;
open UNIQ, "< $opts{fa}" or die "Error:Cannot open file  $opts{fa} : $! \n";
while(<UNIQ>){
	chomp;
	if(/^>(\S+)/){
		my $name =$1;
		$name=~/^(\w+)\_\d+\_x(\d+)$/;
		$seqs_uniq{$titles{$1}}++;
		$seqs_num{$titles{$1}}+=$2;
	}
}
close UNIQ;

my @chromosome;

open REF, "< $opts{ref_fai}" or die "Error:Cannot open file  $opts{ref_fai} : $! \n";
while(<REF>){
    $_=~s/[\n\r]//g;
	next if ($_ eq ""||/^$/||/^#/);
    my ($chr,$length)=split;
    push(@chromosome, $chr);
}
close REF;


my (%matchs,%matchsall,%sample_match,%sample_matchall);
open ARF, "< $opts{arf}" or die "Error:Cannot open file  $opts{arf} : $! \n";
while(<ARF>){
	chomp;
	my @a=split("\t",$_);
	my $c;
	if($a[5]=~/chromosome/i){
		$c=$a[5];
	}elsif($a[5]=~/^\w{3,}/){
		$c=$a[5];
	}
	else{
            $c=$a[5];
            # $c="chromosome_".$a[5];
	}
	$a[0]=~/^(\w+)\_\w+_x(\d+)$/;
	my $n=$2;
	my $name=$titles{$1};
	$matchs{$c}{$a[10]}{$a[0]}=$n;
	$matchsall{$a[10]}{$a[0]}=$n;
	$matchs{$c}{non}{$a[0]}=$n;
	$matchsall{non}{$a[0]}=$n;
	
			$sample_match{$name}{$c}{$a[10]}{$a[0]}=$n;
			$sample_matchall{$name}{$a[10]}{$a[0]}=$n;
			$sample_match{$name}{$c}{non}{$a[0]}=$n;
			$sample_matchall{$name}{non}{$a[0]}=$n;
	
}
close ARF;


my %hash;
my @chromosome_in=sort(keys(%matchs));

#out all summary
open ALL, ">All_$opts{o}";
print ALL "Chromosome\tTotal_num\tTotal_percent\tForword\tReverse\n";
my $total_all_samples=sum(values(%seqs_num));
print ALL "All_reads\t$total_all_samples\t100%\t-\t-\n";

my $match_total_non=sum(values(%{$matchsall{non}}));
print ALL "All_mapped\t".$match_total_non."\t".nearest(.01, $match_total_non/$total_all_samples * 100)."%\t".sum(values(%{$matchsall{'+'}}))."\t".sum(values(%{$matchsall{'-'}}))."\n";

foreach my $c (@chromosome){
        if(! exists $matchs{$c}){
            next;
        }
	my $match_total_all=sum(values(%{$matchs{$c}{non}}));
	$match_total_all=$match_total_all?$match_total_all:0;
	my $f=sum(values(%{$matchs{$c}{'+'}}));
	$f=$f?$f:0;
	my $r=sum(values(%{$matchs{$c}{'-'}}));
	$r=$r?$r:0;
	print ALL "$c\t".$match_total_all."\t".nearest(.01, $match_total_all/$total_all_samples * 100)."%\t".$f."\t".$r."\n";
	if($c=~/chromosome_(\w+)/){
		$hash{'All'}{$1}{"Forword"}="$f";
		$hash{'All'}{$1}{"Reverse"}="$r";
	}
}
close ALL;

#out sample summary
foreach my $s (sort values(%titles)){
	open OUT,">$s\_$opts{o}";
	print OUT "Chromosome\tUniq_num\tUniq_percent\tForword\tReverse\tTotal_num\tTotal_percent\tForword\tReverse\n";
	my $total_uniq=0;
	$total_uniq=$seqs_uniq{$s};
	$total_uniq=$total_uniq?$total_uniq:0;
	my $total_all=0;
	$total_all=$seqs_num{$s};
	$total_all=$total_all?$total_all:0;
	print OUT "$s\_reads\t$total_uniq\t100%\t-\t-\t$total_all\t100%\t-\t-\n";
	
	my $match_total_uniq_non=scalar(keys(%{$sample_matchall{$s}{non}}));
	my $match_total_all_non=sum(values(%{$sample_matchall{$s}{non}}));
	my $total_forword_all_pos;
		if(exists($sample_matchall{$s}{'+'})){
			$total_forword_all_pos=sum(values(%{$sample_matchall{$s}{'+'}}));
		}else{
			$total_forword_all_pos=0;
		}
	my $total_reverse_all_nag;
		if(exists($sample_matchall{$s}{'-'})){
			$total_reverse_all_nag=sum(values(%{$sample_matchall{$s}{'-'}}));
		}else{
			$total_reverse_all_nag=0;
		}
        my $match_total_pct = 0; 
        if($total_all != 0){
            $match_total_pct = $match_total_all_non/$total_all * 100;
        }
        my $match_total_uniq_pct = 0; 
        if($total_uniq != 0){
            $match_total_uniq_pct = $match_total_uniq_non/$total_uniq * 100;
        }
	print OUT "All_mapped\t".$match_total_uniq_non."\t".nearest(.01, $match_total_uniq_pct)."%\t".scalar(keys(%{$sample_matchall{$s}{'+'}}))."\t".scalar(keys(%{$sample_matchall{$s}{'-'}}))."\t".$match_total_all_non."\t".nearest(.01, $match_total_pct)."%\t".$total_forword_all_pos."\t".$total_reverse_all_nag."\n";
	
	foreach my $c (@chromosome){
		my $match_total_uniq=0;			
		$match_total_uniq=scalar(keys(%{$sample_match{$s}{$c}{non}}));
		next if $match_total_uniq==0;
		
		my $match_total_all=sum(values(%{$sample_match{$s}{$c}{non}}));
		
		my $total_forword;
		if(exists($sample_match{$s}{$c}{'+'})){
			$total_forword=sum(values(%{$sample_match{$s}{$c}{'+'}}));
		}else{
			$total_forword=0;
		}
		my $total_reverse;
		if(exists($sample_match{$s}{$c}{'-'})){
			$total_reverse=sum(values(%{$sample_match{$s}{$c}{'-'}}));
		}else{
			$total_reverse=0;
		}
		
		print OUT "$c\t".$match_total_uniq."\t".nearest(.01, $match_total_uniq/$total_uniq * 100)."%\t".scalar(keys(%{$sample_match{$s}{$c}{'+'}}))."\t".scalar(keys(%{$sample_match{$s}{$c}{'-'}}))."\t".$match_total_all."\t".nearest(.01, $match_total_all/$total_all * 100)."%\t".$total_forword."\t".$total_reverse."\n";
		if($c=~/chromosome_(\w+)/){
			my $c=$1;
			my $ss=$s;
			$ss=~s/\s+/_/g;
			$ss=~s/[^\w]//g;
			$hash{$ss}{$c}{"Forword"}="$total_forword";
			$hash{$ss}{$c}{"Reverse"}="$total_reverse";
		}
	}
	close OUT; 
}

#bar out

if($opts{bar} && keys(%hash)>=2){
	open RS,">cmd.r";
	print RS "library(\"ggplot2\")\n";
	for my $s (keys %hash){
		open DATA,">$s.stat.tmp";
		print DATA "Chr\tReads\tStrand\n";
		for my $c ( sort_by_num (keys %{$hash{$s}}) ){
			print DATA "$c\t$hash{$s}{$c}{\"Forword\"}\tForword\n$c\t-$hash{$s}{$c}{\"Reverse\"}\tReverse\n";
		}
		close DATA;
		print RS"
dat<-read.table(\"$s.stat.tmp\",header=T)
dat\$Chr<-factor(dat\$Chr,levels<-unique(as.character(dat\$Chr)))
pdf (\"$s\_map.pdf\",w=10,h=6)
ggplot(data=dat,aes(x=Chr,y=Reads,fill=Strand))+geom_bar(stat=\"identity\", position=\"identity\")+ggtitle(\"$s mapping result\")
dev.off()
";
	}
	close RS;
	system ("Rscript cmd.r");
	#system ("rm *tmp");
}

sub sort_by_num {
	my @a=@_;
	my (@num,@leter);
	for my $aa(@a){
		if ($aa=~/\D/g){
			push @leter,$aa;
		}else{
			push @num,$aa;
		}
	}
	@num=sort {$a<=>$b} @num;
	my @new=(@num,@leter);
	return (@new);
}

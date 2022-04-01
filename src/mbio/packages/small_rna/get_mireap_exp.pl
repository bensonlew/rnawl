#! /usr/bin/perl -w
use strict;
use warnings;
use Carp;
use Getopt::Long;
use Config::IniFiles;
use List::Util qw ( sum);

my %opts;
my $VERSION="1.0";
GetOptions( \%opts,"exp=s","mrd=s","count=s","norm=s","config=s","up=i","down=i", "h!");

my $usage = <<"USAGE";
       Program : $0
       Version : $VERSION
       Contact : yuntao.guo\@majorbio.com
       Lastest modify:2016-9-6  for mireap
       Discription:
              Usage :perl $0 [options]
                -exp*           result of quantifier.pl, ex: miRNAs_expressed_all_samples_dir.csv
		-mrd*           result of quantifier.pl, ex: ./expression_analyses/expression_analyses_dir/miRBase.mrd
                -count		output count table, defult: known_miR_count.xls
                -norm           outout normalized exp table (TPM), defult: known_miR_norm.xls
		-config*	config.ini	config file for uniqueFastaqs.pl
                -up	        upstream	upstream input when run miRdeep2,default:2
                -down           downstream	downstream input when run miRdeep2,default:5
                -h				Display this usage information
                * 				must be given Argument      
USAGE

die $usage if ( !( $opts{exp} && $opts{mrd} && $opts{config}) || $opts{h} ); 
$opts{count}=$opts{count}?$opts{count}:'known_miR_count.xls';
$opts{norm}=$opts{norm}?$opts{norm}:'known_miR_norm.xls';
$opts{up}=$opts{up}?$opts{up}:2;
$opts{down}=$opts{down}?$opts{down}:5;

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

my %exp;
my @t;
open EXP,"< $opts{exp}" or die "Error:Cannot open file exp $opts{exp} : $! \n";
while(<EXP>){
	chomp;	
	if(/^#/){
		@t=split("\t",$_);
	}else{
		my @a=split("\t",$_);
		my $n=(scalar(@a)-4)/2;
		for(my $i=4;$i<scalar(@a)-$n;$i++){
			if(exists($exp{$a[0]}{$t[$i]})){
				$exp{$a[0]}{$t[$i]}=$a[$i] if($a[$i]>$exp{$a[0]}{$t[$i]});
			}else{
				$exp{$a[0]}{$t[$i]}=$a[$i];
			}
		}
	}	
}
close EXP;

my %seqs;
my %map_seq;
my $name;
my $exp;
my @pos;
open MRD,"<$opts{mrd}" or die "Error:Cannot open mrd file $opts{mrd} : $! \n";
while(<MRD>){
	chomp;
	if(/^>(\S+)/){
		$name=$1;
		undef $exp;
		undef @pos;
	}elsif(/^exp\s+(\S+)$/){
		$exp=$1;
		while($exp=~/[M53]+/g){
			my $end=pos($exp);
			my $start=$end-length($&)+1;
			my %hash;
			$hash{start}=$start;
			$hash{end}=$end;
			push(@pos,\%hash);
		}		
	}elsif(/^(\w+\_\d+_x\d+)\s+(\S+)\s+\d*$/){
		my $seq_name=$1;
		if(&check_exp_read(\@pos,$2)){
			$seq_name=~/^(\w+)\_\d+\_x(\d+)$/;
			$map_seq{$1}{$seq_name}=$2;
			$map_seq{'total'}{$seq_name}=$2;
			push(@{$seqs{$seq_name}},$name);
		}
	}
}
close MRD;


sub check_exp_read(){
	my $mpos=shift;
	my $seq=shift;
	$seq=~/[^.]\S+[^.]/g;
	my $end=pos($seq);
	my $start=$end-length($&)+1;	
	foreach my $p (@$mpos){
		if($p->{start}-$opts{up} <=$start && $p->{end}+$opts{down} >= $end){
			return 1;
			last;
		}
	}
	return 0;	
}

open C,">$opts{count}" or die "Error:Cannot create file  $opts{count} : $! \n";

my %mapped_seqs;
print C "miRNA";
foreach my $n (@names){
	print C "\t".$titles{$n};
	$mapped_seqs{$n}=sum(values %{$map_seq{$n}});
} 
print C "\n";

foreach my $m (sort keys(%exp)){
	next if sum(values %{$exp{$m}})==0;
	print C $m;
	foreach my $n (@names){
		my $value=int($exp{$m}{$n});
		print C "\t".$value;
	}
	print C "\n";
}
close C;

open N,"> $opts{norm}" or die "Error:Cannot create file $opts{norm} : $! \n";
print N "miRNA";
foreach my $n (@names){
	print N "\t".$titles{$n};
} 
print N "\n";
foreach my $m (sort keys(%exp)){
	next if sum(values %{$exp{$m}})==0;
	print N "$m";
	foreach my $n (@names){
		my $value=$mapped_seqs{$n}==0?0:$exp{$m}{$n} * 1000000 /$mapped_seqs{$n};
		print N "\t". sprintf("%.2f",$value) ;
	}
	print N "\n";

}
close N;

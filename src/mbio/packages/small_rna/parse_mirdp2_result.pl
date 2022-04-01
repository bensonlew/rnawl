#! /usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Config::IniFiles;
use List::Util qw (sum);

my %opts;
GetOptions (\%opts,"exp=s","prediction=s","config=s","count=s","norm=s","mature=s","tr=s","pfa=s","mfa=s","o=s","h!");
&USAGE if (!( $opts{exp} && $opts{prediction} && $opts{config} && $opts{o}) || $opts{h});

$opts{count}=$opts{count}?$opts{count}:"novel_miR_count.xls";
$opts{norm}=$opts{norm}?$opts{norm}:"novel_miR_norm.xls";
$opts{mature}=$opts{mature}?$opts{mature}:"novel_miR_mature_infor.xls";
$opts{tr}=$opts{tr}?$opts{tr}:"novel_miR_pre.tr";
$opts{pfa}=$opts{pfa}?$opts{pfa}:"novel_miR_pre.fa";
$opts{mfa}=$opts{mfa}?$opts{mfa}:"novel_miR_mature.fa";

#step 1 store sample names--3 letters names into %titles
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
		$titles{$f}=$name;
	}
}

#step 2 read known count table to get sample total counts (as the novel miRNA num is generally small,use the novel total count do normalization is not credible)
my %total;
open KNOWN,"<$opts{exp}" || die "can't open known exp file, please check the file name!\n";
my $hh=<KNOWN>;
chomp $hh;
my @hs=split /\t/,$hh;
while (<KNOWN>){
	my @line=split /\t/;
	my $end=1+ scalar(@hs)/2;
	for (my $k=4;$k<=$end;$k++){
		if (exists $total{$hs[$k]}{$line[0]}  and $total{$hs[$k]}{$line[0]}>=$line[$k] ){
			next;
		}else{
			$total{$hs[$k]}{$line[0]} = $line[$k];
		}
	}
}
close KNOWN;

#step3 parse prediction file and get novel miRNA info
my (%pre,%infor_mature,%mature_seq,%pre_seq,%pre_fa,%map_seq,%infor_pre,$exp,$pre_name,$sequence);
open PREDICTION,"<$opts{prediction}" || die "can't open prediction file, please check the file name!\n";
while (<PREDICTION>){
    chomp;
    my @line=split /\t/;
    my $miRNA_id=$line[10];
    my $pre_miRNA_id=$line[3];
    my $miRNA_start=(split /\.\./, $line[4])[0];
    my $miRNA_end=(split /\.\./, $line[4])[1];
    my $premiRNA_start=(split /\.\./, $line[5])[0];
    my $premiRNA_end=(split /\.\./, $line[5])[1];
    my @reads=split /;/,$line[11];
    $pre{$miRNA_id}=$pre_miRNA_id;
    $mature_seq{$miRNA_id}=$line[6];
    $infor_mature{$miRNA_id}=$line[0]."\t".$line[1]."\t".$miRNA_start."\t".$miRNA_end;
    $infor_pre{$pre_miRNA_id}=$line[0]."\t".$premiRNA_start."\t".$premiRNA_end."\t".$line[1]."\t".length($line[7])."\t".$line[8];
    $pre_seq{$pre_miRNA_id}=$line[9];
	$pre_fa{$pre_miRNA_id}=$line[7];
	foreach my $read(@reads){
	    $read=~/^(\S{3})_\d+_x(\d+)$/;
	    my $seq_name=$1;
	    my $count=$2;
	    $map_seq{$miRNA_id}{$seq_name}+=$count;
	    $total{$seq_name}{$miRNA_id}+=$count;
	}
}

# step 4 print output
# out norm and count table
my $header;
my @s=sort (keys %titles);
for my $na(@s){
	$header.="\t$titles{$na}";
}

system("mkdir -p $opts{o}");
open COUNT,">$opts{o}/$opts{count}" || die "can't open count table\n";
open NORM,">$opts{o}/$opts{norm}" || die "can't open norm table\n";

print COUNT "miRNA$header\n";
print NORM "miRNA$header\n";
for my $id (sort keys %map_seq){
	print COUNT $id;
	print NORM $id;
	for my $sample(@s){
		if (exists $map_seq{$id}{$sample}){
			print COUNT "\t$map_seq{$id}{$sample}";
			my $norm_exp = 0 ;
			if(sum (values %{$total{$sample}}) == 0){
			}else{
			    $norm_exp=sprintf ("%.2f",$map_seq{$id}{$sample}*1000000/sum (values %{$total{$sample}}));
			}
			
			print NORM "\t$norm_exp";
		}else {
			print COUNT "\t0";
			print NORM "\t0";
		}
	}
	print COUNT "\n";
	print NORM "\n";

}
close COUNT;
close NORM;


# out pre.structure and pre.fa
open STR,">$opts{o}/$opts{tr}" || die "can't open str file\n";
open PSEQ,">$opts{o}/$opts{pfa}" || die "can't open fa file\n";
for my $p_id (sort keys %pre_seq){
	print STR ">$p_id\n$pre_seq{$p_id}\n";
	print PSEQ ">$p_id\n$pre_fa{$p_id}\n";
}
close STR;
close PSEQ;

# out mature.fa 
open MSEQ,">$opts{o}/$opts{mfa}" || die "can't open out mature fa file\n";
for my $k(sort keys %mature_seq){
	print MSEQ ">$k\n$mature_seq{$k}\n";
}
close MSEQ;

#out infor.xls 
open MA,">$opts{o}/$opts{mature}" || die "can't open out mature infor xls\n";
print MA "M_id\tM_chr\tM_strand\tM_start\tM_end\tP_id\tP_chr\tP_start\tP_end\tP_strand\tP_lenth(nt)\tP_energy(kcal/mol)\n";
for my $m(sort keys %pre){
	print MA "$m\t$infor_mature{$m}\t$pre{$m}\t$infor_pre{$pre{$m}}\n";
}
close MA;


sub USAGE {
	die (qq/
	Program : $0
	Contact : caiping.shi\@majorbio.com
	Lastest modify:2019-9-16  for mirdp2
	Discription:
	Usage :perl $0 [options]
	Opinions:
	Basic
			-exp        known_exp_file      i.e. miRNAs_expressed_all_samples_dir.csv
			-prediction     prediction file     i.e. filter_P_prediction
			-config		config file		i.e. config.ini
			-o		outputdir		i.e. outdir_prefix
        Opinional
			-count		count.matrix		default:novel_miR_count.xls
			-norm		norm.matrix		default:novel_miR_norm.xls
			-mature		mature miRNA infor	default:novel_miR_mature_infor.xls
			-str		2d_structure		default:novel_miR_pre.tr
			-pfa		precursor.fa		default:novel_miR_pre.fa
			-mfa		mature.fa		default:novel_miR_mature.fa
			-h		usage			this usage information
	\n/);
}

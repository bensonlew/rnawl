#!/usr/bin/perl -w
#xiaoyue.wang@majorbio.com
use strict;
use warnings;
use Getopt::Long;

my %opts;
GetOptions(\%opts, "faa=s", "signalp=s", "tmhmm=s", "outfile=s", "help!");
my $usage =<< "USAGE";
Program: $0
Discription: The result of the combination of TMHMM and signalp.
Options:
	-faa      <file>	faa file of gene
	-signalp  <file>	signalp software outfile
	-tmhmm    <file>        tmhmm software outfile
	-outfile  <file>	The output file name
	-help           	print help information
Example:
	perl signalp_tmhmmtxt2xls.pl -faa V1.faa -signalp V1.signalp.txt -tmhmm V1.tmhmm.txt -outfile V1.signalp.tmhmm.xls
	perl signalp_tmhmmtxt2xls.pl -faa V1.faa -signalp V1.signalp.txt -outfile V1.signalp.xls
	perl signalp_tmhmmtxt2xls.pl -faa V1.faa -tmhmm V1.tmhmm.txt -outfile V1.tmhmm.xls
USAGE

die $usage if($opts{help});
die "faa is must have!!!\n\n$usage" if(! defined $opts{faa});
die "signalp and tmhmm must have at least one!!!\n\n$usage" if(! defined $opts{signalp} && ! defined $opts{tmhmm});
die "please give a output file name!!!\n\n$usage" if(! defined $opts{outfile});

open OUT,">$opts{outfile}" || die $!;
my %result;my %gene_len;my @sort_gene;
(-s $opts{faa}) || die "file: $opts{faa} is not exists!\n";
open FAA,"<$opts{faa}" || die "file: $opts{faa} is not exists!\n";
$/=">";<FAA>;
while(<FAA>){
	chomp;
	my @faa_line=split /\n/;
	my $gene_name=(split /\s+/,$faa_line[0])[0];
	my $seq=join("",@faa_line[1..$#faa_line]);
	$seq=~s/\s//g;
	my $seq_len=length($seq);
	$gene_len{$gene_name}=$seq_len;
	push(@sort_gene,$gene_name);
}
close FAA;
$/="\n";

print OUT "##Gene ID\n##Length(aa)\n";
if(defined $opts{signalp}){
	(-s $opts{signalp}) || die "file: $opts{signalp} is not exists!\n";
	open IN1,"<$opts{signalp}" || die "file: $opts{signalp} is not exists!\n";
	print OUT"##SignalP\n##Description: D value is the score that is used to discriminate signal peptides from non-signal peptides;if D > D-cutoff,SignalP=YES\n";
	while(<IN1>){
		chomp;
		if(/^Name=/){
			my @signalp_line=split /\s+/;
			my $signalp_id=shift@signalp_line;
			$signalp_id=~s/Name=//g;
			my $signalp=shift@signalp_line;
			$signalp=~s/SP=//g;
			$signalp=~s/\'//g;
			my $signalp_opt=join(" ",@signalp_line);
			$result{$signalp_id}{signalp}=$signalp."\t".$signalp_opt;
		}
	}
	close IN1;
}

if(defined $opts{tmhmm}){
	(-s $opts{tmhmm}) || die "file: $opts{tmhmm} is not exists!\n";
	open IN2,"<$opts{tmhmm}" || die "file: $opts{tmhmm} is not exists!\n";
	my @eachgene;
	print OUT "##Number of predicted TMHs: The number of predicted transmembrane helices.\n";
	print OUT "##Exp number of AAs in TMHs: The expected number of amino acids intransmembrane helices. If this number is larger than 18 it is very likely to be a transmembrane protein (OR have a signal peptide).\n";
	print OUT "##Exp number, first 60 AAs: The expected number of amino acids in transmembrane helices in the first 60 amino acids of the protein. If this number more than a few, you should be warned that a predicted transmembrane helix in the N-term could be a signal peptide.\n";
	print OUT "##Total prob of N-in: The total probability that the N-term is on the cytoplasmic side of the membrane.\n";
	print OUT "##POSSIBLE N-term signal sequence: a warning that is produced when \"Exp number, first 60 AAs\" is larger than 10.\n";
	print OUT "##Topology: eg. outside:1-316:表示从1-316个氨基酸在膜外\n";
	while(<IN2>){
		chomp;
		if(/^# (\S+) Length:\s*([0-9\.]+)/){
			my $gene_id=$1;my $length=$2;
			if(@eachgene != 0){
				$result{$eachgene[0]}{tmhmm}=(join "\t",@eachgene[2..$#eachgene]);
			}
			@eachgene = undef;
			$eachgene[0]=$gene_id;$eachgene[1]=$length;
		}
		if(/^# (\S+) Number of predicted TMHs:\s*([0-9\.]+)/){
			my $gene_id=$1;my $TMHs_num=$2;
			($gene_id eq $eachgene[0]) ? ($eachgene[2]=$TMHs_num) : ($eachgene[2]="-");
		}
		if(/^# (\S+) Exp number of AAs in TMHs:\s*([0-9\.]+)/){
			my $gene_id=$1;my $AAs_num=$2;
			($gene_id eq $eachgene[0]) ? ($eachgene[3]=$AAs_num) : ($eachgene[3]="-");
		}
		if(/^# (\S+) Exp number, first 60 AAs:\s*([0-9\.]+)/){
			my $gene_id=$1;my $first60_num=$2;
			($gene_id eq $eachgene[0]) ? ($eachgene[4]=$first60_num) : ($eachgene[4]="-");
			my $next1=<IN2>;
			chomp $next1;
			if($next1=~/^# (\S+) Total prob of N-in:\s*([0-9\.]+)/){
				my $gene_id=$1;my $prob=$2;
				($gene_id eq $eachgene[0]) ? ($eachgene[5]=$prob) : ($eachgene[5]="-");
			}
			if($first60_num>10){
				my $next2=<IN2>;chomp $next2;
				if($next2=~/^# (\S+) POSSIBLE N-term signal sequence/){
					my $gene_id=$1;
					($gene_id eq $eachgene[0]) ? ($eachgene[6]="YES") : ($eachgene[6]="-");
				}
			}else{
				$eachgene[6]="NO";
			}
		}
		if(/^[^#]/){
			my @tmp = split /\s+/;
			($tmp[0] eq $eachgene[0]) ? ($eachgene[7] .= $tmp[2].":".$tmp[3]."-".$tmp[4].";") : ($eachgene[7] ="-");
		}
	}
	close IN2;
	if(@eachgene != 0){
		$result{$eachgene[0]}{tmhmm}=(join "\t",@eachgene[2..$#eachgene]);
	}
}

if(defined $opts{signalp} && defined $opts{tmhmm}){
	print OUT "\n\n\nGene ID\tLength(aa)\tSignalP\tDescription\tNumber of predicted TMHs\tExp number of AAs in TMHs\tExp number, first 60 AAs\tTotal prob of N-in\tPOSSIBLE N-term signal sequence\tTopology\n";
}elsif(defined $opts{signalp} && ! defined $opts{tmhmm}){
	print OUT "\n\n\nGene ID\tLength(aa)\tSignalP\tDescription\n";
}elsif(! defined $opts{signalp} && defined $opts{tmhmm}){
	print OUT "\n\n\nGene ID\tLength(aa)\tNumber of predicted TMHs\tExp number of AAs in TMHs\tExp number, first 60 AAs\tTotal prob of N-in\tPOSSIBLE N-term signal sequence\tTopology\n";
}

for my $k1 (@sort_gene){
	print OUT "$k1\t$gene_len{$k1}";
	for my $k2 (sort {$a cmp $b} keys %{$result{$k1}}){
		print OUT "\t$result{$k1}{$k2}";
	}
	print OUT "\n";
}
close OUT;

##_end_

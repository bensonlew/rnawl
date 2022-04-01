#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fIn,$fOut,$fPos,$PID,$MID,$ParentDepth,$OffDepth);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"vcf:s"=>\$fIn,
	"out:s"=>\$fOut,
	"PID:s"=>\$PID,
	"MID:s"=>\$MID,
			) or &USAGE;
&USAGE unless ($fIn and $fOut and $PID and $MID);
$ParentDepth||=-1;
$OffDepth||=-1;
open In,$fIn;
if ($fIn =~ /.gz/) {
	close In;
	open In,"zcat $fIn|";
}
open Out,">$fOut";
my @indi;
my $n=0;
my %stat;
while (<In>) {
	chomp;
	next if ($_ eq "" || /^$/ || /^##/);
	if (/^#/) {
		(undef,undef,undef,undef,undef,undef,undef,undef,undef,@indi)=split(/\t/,$_);
		my @outid;
		for (my $i=0;$i<@indi;$i++) {
			if ($indi[$i] eq $PID) {
				next;
			}
			if ($indi[$i] eq $MID) {
				next;
			}
			push @outid,$indi[$i];
		}
		die "error input PID or MID" if (!defined $PID || !defined $MID);
		print Out join("\t","#MarkerID","Gtype","Vtype","Pdep","Mdep",@outid),"\n";
	}else{
		my ($chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$format,@geno)=split(/\t/,$_);
		my @alt=split(/,/,join(",",$ref,$alt));
		my %length;
		for (my $i=0;$i<@alt;$i++) {
			$length{length($alt[$i])}=1;
		}
		my $vtype="SNP";
		if (scalar keys %length ne 1) {
			$vtype="INDEL";
		}
		my @format=split(/:/,$format);
		my %info;
		for (my $i=0;$i<@indi;$i++) {
			if ($geno[$i] eq "./.") {
				$info{$indi[$i]}{gt}="./.";
				$info{$indi[$i]}{ad}="0";
				$info{$indi[$i]}{dp}="0";
				next;
			}
			my @info=split(/:/,$geno[$i]);
			for (my $j=0;$j<@info;$j++) {
				$info{$indi[$i]}{gt}=$info[$j] if ($format[$j] eq "GT");
				$info{$indi[$i]}{ad}=$info[$j] if ($format[$j] eq "AD");
				$info{$indi[$i]}{dp}=$info[$j] if ($format[$j] eq "DP");
			}
			if ($info{$indi[$i]}{gt} eq "./.") {
				$info{$indi[$i]}{ad}=0;
				$info{$indi[$i]}{dp}=0;
			}elsif (!exists $info{$indi[$i]}{ad} && exists $info{$indi[$i]}{dp}) {
				my @gt=split(/\//,$info{$indi[$i]}{gt});
				my $max=(sort{$a<=>$b} @gt)[-1];
				my @dep;
				if ($gt[0] eq $gt[1]) {
					for (my $j=0;$j <= $max;$j++) {
						if ($j eq $gt[0]) {
							push @dep,$info{$indi[$i]}{dp};
						}else{
							push @dep,0;
						}
					}
				}else{
					for (my $j=0;$j <= $max;$j++) {
						if ($j eq $gt[0] || $j eq $gt[1]) {
							push @dep,$info{$indi[$i]}{dp}/2;
						}else{
							push @dep,0;
						}
					}
				}
				$info{$indi[$i]}{ad}=join(",",@dep);
			}
		}
		$info{$PID}{gt}||="./.";
		$info{$MID}{gt}||="./.";
		if ($info{$PID}{gt} eq "./." ||$info{$MID}{gt} eq "./." ){$stat{miss}++;;next;};#missing parent
		my ($p1,$p2)=split(/\//,$info{$PID}{gt});
		my ($m1,$m2)=split(/\//,$info{$MID}{gt});
		if ($p1 eq $p2 && $m1 eq $m2 && $m1 eq $p2){$stat{aaxaa}++;next;};#aaxaa
		my @pd=split(/\,/,$info{$PID}{ad});
		my @md=split(/\,/,$info{$MID}{ad});
		my $sum1=$pd[$p1];$sum1+=$pd[$p2] if($p1 ne $p2);
		my $sum2=$md[$m1];$sum2+=$md[$m2] if($m1 ne $m2);
		$stat{"$sum1:$sum2"}++;
		if ($sum1 < $ParentDepth || $sum2 < $ParentDepth){$stat{pmdep}++;next;}
		my %geno;
		$geno{$p1}++;
		$geno{$p2}++;
		$geno{$m1}++;
		$geno{$m2}++;
		my %ale;
		my @out;
		push @out,"$chr\_$pos";
		if (scalar keys %geno == 4) { #abxcd
			 %ale=(
				$p1=>"a",
				$p2=>"b",
				$m1=>"c",
				$m2=>"d",
			);
			push @out,"abxcd";
			$stat{abxcd}++;
		}elsif (scalar keys %geno == 3 && $p1 eq $p2 && $m1 ne $m2) {#ccxab
			 %ale=(
				$p1=>"c",
				$m1=>"a",
				$m2=>"b",
			);
			push @out,"ccxab";
			$stat{ccxab}++;
		}elsif (scalar keys %geno == 3 && $p1 ne $p2 && $m1 eq $m2) {#abxcc
			 %ale=(
				$p1=>"a",
				$p2=>"b",
				$m1=>"c",
			);
			push @out,"abxcc";
			$stat{abxcc}++;

		}elsif (scalar keys %geno == 3 && $p1 ne $p2 && $m1 ne $m2) {#efxeg
			my @ale=sort {$geno{$a}<=>$geno{$b}} keys %geno;
			 %ale=(
				$ale[0]=>"f",
				$ale[1]=>"g",
				$ale[2]=>"e",
			);
			push @out,"efxeg";
			$stat{efxeg}++;

		}elsif (scalar keys %geno == 2 && $info{$PID}{gt} eq $info{$MID}{gt}) {#hkxhk
			%ale=(
				$p1=>"h",
				$p2=>"k",
			);
			push @out,"hkxhk";
			$stat{hkxhk}++;

		}elsif (scalar keys %geno == 2 && $p1 eq $p2 && $m1 ne $m2) {#nnxnp
			my @ale=sort {$geno{$a}<=>$geno{$b}} keys %geno;
			 %ale=(
				$ale[0]=>"p",
				$ale[1]=>"n",
			);
			push @out,"nnxnp";
						$stat{nxxnp}++;

		}elsif (scalar keys %geno == 2 && $p1 ne $p2 && $m1 eq $m2) {#lmxll
			my @ale=sort {$geno{$a}<=>$geno{$b}} keys %geno;
			 %ale=(
				$ale[0]=>"m",
				$ale[1]=>"l",
			);
			push @out,"lmxll";
			$stat{lmxll}++;

		}elsif (scalar keys %geno == 2 && $p1 eq $p2 && $m1 eq $m2) {#aaxbb
			 %ale=(
				$p1=>"a",
				$m1=>"b",
			);
			push @out,"aaxbb";
			$stat{aaxbb}++;
			#if($pos eq "1209058"){print Dumper $p1,"\t",$m1,"\n";print Dumper %info;die};
		}
		push @out,$vtype;
		push @out,$sum1;
		push @out,$sum2;
		for (my $i=0;$i<@indi;$i++) {
			next if ($indi[$i] eq $PID || $indi[$i] eq $MID);
			if ($info{$indi[$i]}{gt} eq "./.") {
				push @out,"0,--";
				next;
			}
			my ($g1,$g2)=split(/\//,$info{$indi[$i]}{gt});
			my @dp=split(/\,/,$info{$indi[$i]}{ad});
			my $sum=$dp[$g1];
			$sum+=$dp[$g2] if($g1 ne $g2);
			if ($sum < $OffDepth) {
				push @out,"0,--";
			}else{
				my ($g1,$g2)=split(/\//,$info{$indi[$i]}{gt});
				if (!exists $ale{$g1} || !exists $ale{$g2}) {
					push @out,"0,--";
				}else{
					my $geno=join("",sort ($ale{$g1},$ale{$g2}));
					push @out,"$sum,$geno";
				}
			}
		}
		print Out join("\t",@out),"\n";
	}
	#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT,@INDI)=split;
}
close Out;
close In;
open Out,">$fOut.log";
# print Dumper \%stat;
foreach my $stat (sort keys %stat) {
	print Out join("\t",$stat,$stat{$stat}),"\n";
}
close Out;
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub USAGE {#
        my $usage=<<"USAGE";
Contact:        long.huang\@majorbio.com;
Script:			$Script
Description:

Usage:
  Options:
  -vcf	<file>	input file vcf file
  -out	<file>	output keys of file name
  -PID	<str>	paternal ID
  -MID	<str>	maternale ID
  -h         Help

USAGE
        print $usage;
        exit;
}

#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fIn,$dOut,$select,$pop);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Spreadsheet::ParseExcel;
use Statistics::Distributions;
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"trit:s"=>\$fIn,
	"out:s"=>\$dOut,
	"select:s"=>\$select,
	"pop:s"=>\$pop,
			) or &USAGE;
&USAGE unless ($fIn and $dOut and $select);

mkdir $dOut if (!-d $dOut);
open Out,">$dOut/total.sample.xls";
my $parser   = Spreadsheet::ParseExcel->new();
my $workbook = $parser->Parse($fIn);
if ( !defined $workbook ) {
	die $parser->error(), ".\n";
}

for my $worksheet ( $workbook->worksheets() ) {
	my $name=$worksheet->get_name();
	my @out=();
	my @line;
	if ( $name eq "Sheet1") {
		my ( $row_min, $row_max ) = $worksheet->row_range();
		my ( $col_min, $col_max ) = $worksheet->col_range();
		
		for my $row ( $row_min .. $row_max ) {
			for my $col ( $col_min .. $col_max ) {
				my $cell = $worksheet->get_cell( $row, $col );
				next unless $cell;
				if($row==0 and $col==0){
					print Out "Sample ID\t";
					next;
				}else{
					print Out $cell->value(),"\t";
				}
			}
			print Out "\n";
		}
		close Out;
	}
}

my (@sample,@trits);
my %stat;
my %region;
open Out,">$dOut/annovar.sample.xls";
open File,$select;
while(<File>){
	chomp;
	next if($_ eq ""|| /^$/);
	push @sample,$_;
}
close File;
open In,"$dOut/total.sample.xls";
while(<In>){
	chomp;
	my($id,@undi)=split(/\t/,$_);
	if($id=~/Sample/){
		@trits=@undi;
		print Out $_,"\n";
	}else{
		foreach my $sample(@sample){
			if($id eq $sample){
				print Out $_,"\n";
				for (my $i=0;$i<@trits;$i++) {
					my $trits=$trits[$i];
					$stat{$trits}{$sample}=$undi[$i];
					if($trits=~/(\w+)(\W+)(\d+)/){
						my$tri=$1;
						push @{$region{$tri}{$sample}},$undi[$i];
					}
				}
			}
		}
	}
}
close In;
close Out;

my %hash;
my %info;
mkdir "$dOut/trit" if(!-d "$dOut/trit");
my $nind=scalar@sample;
foreach my $tri(sort keys %region){
	open OUT,">$dOut/trit/$tri.txt";
	if($pop eq "CP"){
		print OUT "ntrt=2\nnind=$nind\nmiss=NaN\nsampleID\t$tri\n";
	}
	my @id=();
	my @tvalue=();
	foreach my$sample(sort keys %{$region{$tri}}){
		my $num=scalar@{$region{$tri}{$sample}};
		my $pva;
		foreach(@{$region{$tri}{$sample}}){
			$pva+= $_;
		}
		if($pop eq "CP"){
			print OUT "$sample\t",sprintf("%.3f",$pva/$num),"\n";
		}
		push @id,$sample;
		push @tvalue,sprintf("%.3f",$pva/$num);
	}
	if($pop ne "CP"){
		print OUT "Genotype\t",join(",",@id),"\n$tri",join(",",@tvalue),"\n";
	}
	close OUT;
}
foreach my $trits(sort keys %stat){
	my @info=();
	open Out,">$dOut/trit/$trits.txt";
	if($pop eq "CP"){
		print Out "ntrt=2\nnind=$nind\nmiss=NaN\nsampleID\t$trits\n";
	}
	my @samplename=();
	my @tritvalue=();
	foreach my $sample(sort keys %{$stat{$trits}}) {
		push @info,$stat{$trits}{$sample};
		if($pop eq "CP"){
			print Out "$sample\t$stat{$trits}{$sample}\n";
		}else{
			push @samplename,$sample;
			push @tritvalue,$stat{$trits}{$sample};
		}
	}
	if($pop ne "CP"){
		print Out "Genotype\t",join(",",@samplename),"\n$trits",join(",",@tritvalue),"\n";
	}
	close Out;
	my $trit;
	if($trits=~/(\w+)(\W+)(\d+)/) {
		$trit=$1;
	}else{
		$trit=$trits;
	}
	my $value=join("\,",@info);
	push @{$hash{$trit}},$value;
}

open ANNO,">$dOut/annovar.trit.xls";
print ANNO "Trit\tdf\tMS\tF Value\tP Value\n";
foreach my$trit(sort keys%hash){
	
	my $out=&S2($trit,\@{$hash{$trit}});
	print ANNO "$trit\t$out\n";
}
close ANNO;


#######################################################################################
#print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub S2 {#
	my($trit,$value)=@_;
	my($g,$df,$ms,$f,$p);
	my $sum=0;
	my @eve;
	my $r=scalar(@$value);
	my %per;
	foreach (@$value){
		my @arry=split(/\,/,$_);
		$g=0;
		for(my$m=0;$m<scalar@arry;$m++){
			if($arry[$m] ne "NaN"){
				push @{$per{$m}},$arry[$m];
				push @eve,$arry[$m];
				$g++ ;
				$sum=$sum + $arry[$m];
			}else{
				next;
			}
		}
	}
	$df = $g - 1;
	my $yii=sprintf("%.4f",($sum/($g*$r)));
	my ($ssg,$ssr);
	if($r eq "1"){
		foreach my$m(keys %per){
			my $yik;
			foreach(@{$per{$m}}){
				$yik=$_;
			}
			$ssg+=($yik - $yii)**2;
		}
		$ssr=0;
		$ms=sprintf("%.4f",$ssg/$df);
		$f="--";
		$p="--";
	}else{
		foreach my$m(keys %per){
			my $yi=0;
			foreach(@{$per{$m}}){
				$yi=$yi + $_;
			}
			$yi=sprintf("%.4f",$yi/(scalar@{$per{$m}}));
			
			foreach(@{$per{$m}}){
				my $y2=($_ - $yi)**2;
				$ssr+= $y2;
			}
			$ssg+=($yi - $yii)**2;
		}
		my $msg=$ssg*$r/$df;
		my $num=($r - 1)*$g;
		my $msr=$ssr/$num;
		$f=$msg/$msr;
		$msg=sprintf("%.4f",$msg);
		$msr=sprintf("%.4f",$msr);
		$f=sprintf("%.4f",$f);
		$ms="$msg\/$msr";
		$p=sprintf("%.4f",(Statistics::Distributions::fprob ($df,$num,$f)));
	}
	my $out=join("\t",$df,$ms,$f,$p);
	return $out;
}

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
  -trit	<file>	input trit file(*.txt)
  -out	<file>	output result file
  -pop	<stri>	pop type
  -select	<file>	trit annovar's sample
  -h         Help

USAGE
        print $usage;
        exit;
}

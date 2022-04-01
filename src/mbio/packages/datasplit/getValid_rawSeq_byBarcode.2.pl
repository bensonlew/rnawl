#!/usr/bin/perl -w
use strict;
use FindBin qw/$Bin/;

die "usage: perl $0 raw.1.fq raw.2.fq barcode.config out.prefix\n" unless(@ARGV == 4);
my ($fq1,$fq2,$splitConfig,$prefix) = @ARGV;

#my (%total_Fbarcode,%total_Rbarcode);
my (%Fbarcode,%Rbarcode,%sample,%barcodeTag);


my (%check_sam,%check_bar);
open INL,$splitConfig or die "$!\n";
while(<INL>){
	chomp;
	next if(/^#/);
	my @temp = split;
	die "error: repeat sample $_\n" if(exists $check_sam{$temp[0]});
	die "error: repeat barcode $_\n" if(exists $check_bar{$temp[2]} or exists $check_bar{$temp[3]});
	$check_sam{$temp[0]} = 1;$check_bar{$temp[2]} = 1;$check_bar{$temp[3]} = 1;
	$barcodeTag{$temp[1]} = 1;
	$Fbarcode{$temp[1]} = $temp[2];
	$Rbarcode{$temp[1]} = $temp[3];
	$sample{$temp[1]} = $temp[0];
}
close INL;

if($fq1 =~ /gz$/){
	open INF1,"gzip -dc $fq1 |" or die "$!\n";
}else{
	open INF1,$fq1 or die "$!\n";
}
if($fq2 =~ /gz$/){
	open INF2,"gzip -dc $fq2 |" or die "$!\n";
}else{
	open INF2,$fq2 or die "$!\n";
}


my ($totalSeqNum,$validSeqNum,$chimericNum,$noRbarcodeNum,$noFbarcodeNum);
$totalSeqNum = 0;$validSeqNum = 0;$chimericNum = 0;$noRbarcodeNum = 0;$noFbarcodeNum = 0;
my %eachSampleValid;
open OUT1,"> $prefix.valid.1.fq" or die "$!\n";
open OUT2,"> $prefix.valid.2.fq" or die "$!\n";
open OUT3,"> $prefix.discard.1.fq" or die "$!\n";
open OUT4,"> $prefix.discard.2.fq" or die "$!\n";
open OUS,"> $prefix.seq2sam.stat" or die "$!\n";
print OUS "#Sequence\tSample\n";
while(<INF1>){
	my $flag = 0;my $flag2 = 0;
	$totalSeqNum++;
	chomp;
	chomp(my $seq1 = <INF1>);
	chomp(my $direction1 = <INF1>);
	chomp(my $quality1 = <INF1>);
	chomp(my $head2 = <INF2>);
	chomp(my $seq2 = <INF2>);
	chomp(my $direction2 = <INF2>);
	chomp(my $quality2 = <INF2>);
	foreach my $tag(sort keys %barcodeTag){
#		my $randBaseNum;
		my ($fRandBaseNum,$rRandBaseNum);
		if($tag =~ /^(\d)(\d)/){
			$fRandBaseNum = $1; $rRandBaseNum = $2;
		}elsif($tag =~ /^(\d)/){
			$fRandBaseNum = $1; $rRandBaseNum = $1;
		}else{
			die "barcode tag error: $tag\n";
		}
		if($seq1 =~ /^(\w{$fRandBaseNum})$Fbarcode{$tag}/){
			my $randSeq1 = $1;
			$flag = 1;
			if($seq2 =~ /^(\w{$rRandBaseNum})$Rbarcode{$tag}/){
				$validSeqNum++;
				my $randSeq2 = $1;
				my $sub_seq1 = substr($seq1,$fRandBaseNum);
				my $sub_quality1 = substr($quality1,$fRandBaseNum);
				my $sub_seq2 = substr($seq2,$rRandBaseNum);
				my $sub_quality2 = substr($quality2,$rRandBaseNum);
				print OUT1 "$_\n$sub_seq1","$randSeq1\n$direction1\n$sub_quality1","#" x $fRandBaseNum,"\n";
				print OUT2 "$head2\n$sub_seq2","$randSeq2\n$direction2\n$sub_quality2","#" x $rRandBaseNum,"\n";
				my $temp_name = (split /\s+/,$_)[0];
				print OUS "$temp_name\t$sample{$tag}\n";
				$eachSampleValid{"$sample{$tag}\t$tag"} ++;
			}else{
				foreach my $Rtag(sort keys %Rbarcode){
					my $rRandBaseNum2;
					if($Rtag =~ /^(\d)(\d)/){
						$rRandBaseNum2 = $2;
					}elsif($Rtag =~ /^(\d)/){
						$rRandBaseNum2 = $1;
					}else{
						die "barcode tag error: $Rtag\n";
					}
					if($seq2 =~ /^(\w{$rRandBaseNum2})$Rbarcode{$Rtag}/){
						$flag2 = 1;
						print OUT3 "$_\tFbarcode=$tag\tRbarcode=$Rtag\n$seq1\n$direction1\n$quality1\n";
						print OUT4 "$head2\tFbarcode=$tag\tRbarcode=$Rtag\n$seq2\n$direction2\n$quality2\n";
						$chimericNum++;
						last;
					}
				}
				unless($flag2){
					$noRbarcodeNum++;
					print OUT3 "$_\tnoRbarcode\n$seq1\n$direction1\n$quality1\n";
					print OUT4 "$head2\tnoRbarcode\n$seq2\n$direction2\n$quality2\n";
				}
			}
			last;
		}elsif($seq2 =~ /^(\w{$fRandBaseNum})$Fbarcode{$tag}/){
			$flag = 1;
			my $randSeq2 = $1;
			if($seq1 =~ /^(\w{$rRandBaseNum})$Rbarcode{$tag}/){
				$validSeqNum++;
				my $randSeq1 = $1;
				my $sub_seq1 = substr($seq1,$rRandBaseNum);
				my $sub_quality1 = substr($quality1,$rRandBaseNum);
				my $sub_seq2 = substr($seq2,$fRandBaseNum);
				my $sub_quality2 = substr($quality2,$fRandBaseNum);
				#print OUT1 "$_\n$sub_seq1","$randSeq1\n$direction1\n$sub_quality1","#" x $rRandBaseNum,"\n";
				print OUT2 "$_\n$sub_seq1","$randSeq1\n$direction1\n$sub_quality1","#" x $rRandBaseNum,"\n";
				#print OUT2 "$head2\n$sub_seq2","$randSeq2\n$direction2\n$sub_quality2","#" x $fRandBaseNum,"\n";
				print OUT1 "$head2\n$sub_seq2","$randSeq2\n$direction2\n$sub_quality2","#" x $fRandBaseNum,"\n";
				my $temp_name = (split /\s+/,$_)[0];
				print OUS "$temp_name\t$sample{$tag}\n";
				$eachSampleValid{"$sample{$tag}\t$tag"} ++;
			}else{
				foreach my $Rtag(sort keys %Rbarcode){
					my $rRandBaseNum2;
					if($Rtag =~ /^(\d)(\d)/){
						$rRandBaseNum2 = $2;
					}elsif($Rtag =~ /^(\d)/){
						$rRandBaseNum2 = $1;
					}else{
						die "barcode tag error: $Rtag\n";
					}
					if($seq1 =~ /^(\w{$rRandBaseNum2})$Rbarcode{$Rtag}/){
						$flag2 = 1;
						print OUT3 "$_\tFbarcode=$tag\tRbarcode=$Rtag\n$seq1\n$direction1\n$quality1\n";
						print OUT4 "$head2\tFbarcode=$tag\tRbarcode=$Rtag\n$seq2\n$direction2\n$quality2\n";
						$chimericNum++;
						last;
					}
				}
				unless($flag2){
					$noRbarcodeNum++;
					print OUT3 "$_\tnoRbarcode\n$seq1\n$direction1\n$quality1\n";
					print OUT4 "$head2\tnoRbarcode\n$seq2\n$direction2\n$quality2\n";
				}
			}
			last;
		}
	}
	next if($flag);
	$noFbarcodeNum++;
	print OUT3 "$_\tnoFbarcode\n$seq1\n$direction1\n$quality1\n";
	print OUT4 "$head2\tnoFbarcode\n$seq2\n$direction2\n$quality2\n";
}
close INF1;
close INF2;
close OUT1;
close OUT2;
close OUT3;
close OUT4;
close OUS;

open OUST,"> $prefix.valid.stat" or die "$!\n";
print OUST "#RawNum\tChimericNum\tNoFbarcodeNum\tNoRbarcodeNum\tValidNum\tRate\n";
print OUST "$totalSeqNum\t$chimericNum\t$noFbarcodeNum\t$noRbarcodeNum\t$validSeqNum\t",$validSeqNum/$totalSeqNum,"\n";
close OUST;

open OUTL, "> $prefix.valid.eachsample.log" or die "$!\n";
print OUTL "#Sample\tBarcode\tValidReads\n";
foreach my $key (sort keys %eachSampleValid){
    print OUTL "$key\t$eachSampleValid{$key}\n";
}
close OUTL;

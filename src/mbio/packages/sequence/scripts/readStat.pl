#!/usr/bin/perl -w
use strict;
use warnings;

die "usage: perl $0 rawStat.list cleanFq.list prefix
or	perl $0 rawStat.list QC.fq.list cleanFq.list prefix\n" unless(@ARGV == 3 or @ARGV == 4);
my ($rawList,$QCFqList,$cleanFqList,$prefix);
if(@ARGV == 4){
	($rawList,$QCFqList,$cleanFqList,$prefix) = @ARGV;
}else{
	($rawList,$cleanFqList,$prefix) = @ARGV;
}

my (%rawRead,%rawBase);
open INR,$rawList or die "$!\n";
while(my $file = <INR>){
	chomp $file;
	my $sam_ = (split /\//,$file)[-1];
	my $sam = (split /\./,$sam_)[0];
	my ($readNum,$baseNum);
	open INF,$file or die "$!\n";
	while(<INF>){
		chomp;
		next if(/^column/);
		my @temp = split;
		$readNum = $temp[-1];
		$baseNum += $temp[1]
	}
	close INF;
	$rawRead{$sam} += $readNum;
	$rawBase{$sam} += $baseNum;
}
close INR;

open OUR,"> $prefix.rawData.stat" or die "$!\n";
print OUR "#Sample\tReadsNum\tBaseNum\tAverageLength\n";
foreach my $sam(sort keys %rawRead){
	print OUR "$sam\t$rawRead{$sam}\t$rawBase{$sam}\t",$rawBase{$sam}/$rawRead{$sam},"\n";
}
close OUR;

my (%QCRead1,%QCRead2,%QCReads,%QCBase1,%QCBase2,%QCBases);
if($QCFqList){
#	my (%QCRead1,%QCRead2,%QCReads,%QCBase1,%QCBase2,%QCBases);
	open INQ,$QCFqList or die "$!\n";
	while(my $file = <INQ>){
		chomp $file;
		my @temp = split /\//,$file;
		my $sam_name = (split /\./, $temp[-1])[0]; 
		print "read $file ...\n";
		open INF,$file or die "$!\n";
		if($temp[-1] =~ /1\.fq$/ or $temp[-1] =~ /1\.fastq$/){
			while(<INF>){
				chomp;
				die "error: fastq format error $_\n" unless(/^@/);
				chomp(my $seq = <INF>);
				<INF>;<INF>;
				$QCRead1{$sam_name}++;
				$QCBase1{$sam_name} += length($seq);
			}
		}elsif($temp[-1] =~ /2\.fq/ or $temp[-1] =~ /2\.fastq/){
			while(<INF>){
				chomp;
				die "error: fastq format error $_\n" unless(/^@/);
				chomp(my $seq = <INF>);
				<INF>;<INF>;
				$QCRead2{$sam_name}++;
				$QCBase2{$sam_name} += length($seq);
			}
		}elsif($temp[-1] =~ /s\.fq/ or $temp[-1] =~ /s\.fastq/){
			while(<INF>){
				chomp;
				die "error: fastq format error $_\n" unless(/^@/);
				chomp(my $seq = <INF>);
				<INF>;<INF>;
				$QCReads{$sam_name}++;
				$QCBases{$sam_name} += length($seq);
			}
		}else{
			die "error: fastq file name must be *.[12s].[fq|fastq] : $file\n";
		}
		close INF;
	}
	close INQ;
}

my (%cleanRead1,%cleanRead2,%cleanReads,%cleanBase1,%cleanBase2,%cleanBases);
open INC,$cleanFqList or die "$!\n";
while(my $file = <INC>){
	chomp $file;
	my @temp = split /\//,$file;
	my $sam_name = (split /\./, $temp[-1])[0];
	print "read $file ...\n";
	open INF,$file or die "$!\n";
	if($temp[-1] =~ /1\.fq$/ or $temp[-1] =~ /1\.fastq$/){
		while(<INF>){
			chomp;
			die "error: fastq format error $_\n" unless(/^@/);
			chomp(my $seq = <INF>);
			<INF>;<INF>;
			$cleanRead1{$sam_name}++;
			$cleanBase1{$sam_name} += length($seq);
		}
	}elsif($temp[-1] =~ /2\.fq/ or $temp[-1] =~ /2\.fastq/){
		while(<INF>){
			chomp;
			die "error: fastq format error $_\n" unless(/^@/);
			chomp(my $seq = <INF>);
			<INF>;<INF>;
			$cleanRead2{$sam_name}++;
			$cleanBase2{$sam_name} += length($seq);
		}
	}elsif($temp[-1] =~ /s\.fq/ or $temp[-1] =~ /s\.fastq/){
		while(<INF>){
			chomp;
			die "error: fastq format error $_\n" unless(/^@/);
			chomp(my $seq = <INF>);
			<INF>;<INF>;
			$cleanReads{$sam_name}++;
			$cleanBases{$sam_name} += length($seq);
		}
	}else{
        next;
		#die "error: fastq file name must be *.[12s].[fq|fastq] : $file\n";
	}
	close INF;
}
close INC;

if($QCFqList){
	open OUQ,"> $prefix.QCData.stat" or die "$!\n";
	print OUQ "#Sample\tReadsNum\tBasesNum\tAverageLength\n";
	foreach my $sam(sort keys %QCRead1){
		#print OUQ "$sam\tPair1\t$QCRead1{$sam}\t$QCBase1{$sam}\t",$QCBase1{$sam}/$QCRead1{$sam},"\n";
		#print OUQ "$sam\tPair2\t$QCRead2{$sam}\t$QCBase2{$sam}\t",$QCBase2{$sam}/$QCRead2{$sam},"\n";
		#print OUQ "$sam\tSingle\t$QCReads{$sam}\t$QCBases{$sam}\t",$QCBases{$sam}/$QCReads{$sam},"\n";
		print OUQ "$sam\t",$QCRead1{$sam} + $QCRead2{$sam} + $QCReads{$sam},"\t",$QCBase1{$sam} + $QCBase2{$sam} + $QCBases{$sam},"\t",($QCBase1{$sam} + $QCBase2{$sam} + $QCBases{$sam})/($QCRead1{$sam} + $QCRead2{$sam} + $QCReads{$sam}),"\n";
	}
}
close OUQ;

open OUC,"> $prefix.cleanData.stat" or die "$!\n";
print OUC "#Sample\tReadsNum\tBasesNum\tAverageLength\n";
foreach my $sam(sort keys %cleanRead1){
	#print OUC "$sam\tPair1\t$cleanRead1{$sam}\t$cleanBase1{$sam}\t",$cleanBase1{$sam}/$cleanRead1{$sam},"\n";
	#print OUC "$sam\tPair2\t$cleanRead2{$sam}\t$cleanBase2{$sam}\t",$cleanBase2{$sam}/$cleanRead2{$sam},"\n";
	#print OUC "$sam\tSingle\t$cleanReads{$sam}\t$cleanBases{$sam}\t",$cleanBases{$sam}/$cleanReads{$sam},"\n";
	print OUC "$sam\t",$cleanRead1{$sam} + $cleanRead2{$sam} + $cleanReads{$sam},"\t",$cleanBase1{$sam} + $cleanBase2{$sam} + $cleanBases{$sam},"\t",($cleanBase1{$sam} + $cleanBase2{$sam} + $cleanBases{$sam})/($cleanRead1{$sam} + $cleanRead2{$sam} + $cleanReads{$sam}),"\n";
}
close OUC;


#open OUS,"> $prefix.QC.stat" or die "$!\n";
#if($QCFqList){
#	print OUS "#Sample\tRawReads\tQCReads\tCleanReads\tRate\tRawBases\tQCBases\tCleanBases\tRate\n";
#	foreach my $sam(sort keys %rawRead){
#		print OUS "$sam\t$rawRead{$sam}\t",$QCRead1{$sam} + $QCRead2{$sam} + $QCReads{$sam},"\t",$cleanRead1{$sam} + $cleanRead2{$sam} + $cleanReads{$sam},"\t",($cleanRead1{$sam} + $cleanRead2{$sam} + $cleanReads{$sam})/$rawRead{$sam},"\t$rawBase{$sam}\t",$QCBase1{$sam} + $QCBase2{$sam} + $QCBases{$sam},"\t",$cleanBase1{$sam} + $cleanBase2{$sam} + $cleanBases{$sam},"\t",($cleanBase1{$sam} + $cleanBase2{$sam} + $cleanBases{$sam})/$rawBase{$sam},"\n";
#	}
#}else{
#	print OUS "#Sample\tRawReads\tCleanReads\tRate\tRawBases\tCleanBases\tRate\n";
#	foreach my $sam(sort keys %rawRead){
#		print OUS "$sam\t$rawRead{$sam}\t",$cleanRead1{$sam} + $cleanRead2{$sam} + $cleanReads{$sam},"\t",($cleanRead1{$sam} + $cleanRead2{$sam} + $cleanReads{$sam})/$rawRead{$sam},"\t$rawBase{$sam}\t",$cleanBase1{$sam} + $cleanBase2{$sam} + $cleanBases{$sam},"\t",($cleanBase1{$sam} + $cleanBase2{$sam} + $cleanBases{$sam})/$rawBase{$sam},"\n";
#	}
#}
#close OUS;

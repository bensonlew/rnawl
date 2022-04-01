#!/usr/bin/perl 
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fIn,$fOut,$fMap,$Pos);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$fIn,
				"m:s"=>\$fMap,
				"o:s"=>\$fOut,
				) or &USAGE;
&USAGE unless ($fIn and $fOut);
my %SORT;
my $nmarker=0;
if (!$Pos) {
	open In,$fMap;
	my $GID;
	while (<In>) {
		chomp;
		next if ($_ eq ""||/^$/ ||/^;/ );
		if (/group/) {
			$GID=(split(/\s+/,$_))[-1];
			next;
		}
		my ($id,$pos)=split(/\s+/,$_);
		push @{$SORT{$GID}},$id;
		$nmarker++;
	}
	close In;
}else{
	open In,$fMap;
	while (<In>) {
		chomp;
		next if ($_ eq ""||/^$/ ||/^;/);
		my ($id,$chr,$pos)=split(/\s+/,$_);
		push @{$SORT{$chr}},$id;
		$nmarker++;
	}
	close In;
}
open In,$fIn;
my %Geno;
my @Head;
my $nind;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/);
	my ($id,@info)=split;
	 if(scalar @info < 3 || /MarkerID/ || $info[0] eq ""){
		next if(/number_of_loci/);
		push @Head,$_;
		next;
	}
	$nind=scalar @info;
	for (my $i=0;$i<@info;$i++) {
		$Geno{"indi$i"}{$id}=$info[$i];
	}
}
close In;
my $win=15;
my $step=1;
my %OGeno;
foreach my $Chr (sort keys %SORT) {
	my @Order=@{$SORT{$Chr}};
	for (my $i=0;$i<@Order;$i++) {
		my $start=$i-$win;
		my $end=$i+$win;
		$start=0 if ($start < 0);
		$end=scalar @Order -1 if ($end > scalar @Order -1);
		if (abs($start - $i) != abs($end - $i)) {
			my $winshift1=abs($start-$i);
			my $winshift2=abs($end-$i);
			my $winshift=$winshift1;
			$winshift=$winshift2 if ($winshift2 < $winshift1);
			$start=$i-$winshift;
			$end=$i+$winshift;
		}
		foreach my $indi (sort keys %Geno) {
			my %stat;
			my $n=0;
			for (my $j=$start;$j<=$end;$j++) {
				
				if (!exists $Geno{$indi}{$Order[$j]}) {
					next;
				}
				next if ($Geno{$indi}{$Order[$j]} eq "U");
				$stat{$Geno{$indi}{$Order[$j]}}++;
				$n++;
			}
			if (scalar %stat == 0) {
				$OGeno{$indi}{$Order[$i]}="U";
			}else{
				my @Geno=sort{$stat{$b}<=>$stat{$a}} keys %stat;
				if (scalar @Geno >1 && $stat{$Geno[0]} == $stat{$Geno[1]}) {
					$OGeno{$indi}{$Order[$i]}="U";
				}else{
					$OGeno{$indi}{$Order[$i]}=$Geno[0];
					#$OGeno{$indi}{$Order[$i]}="U" if($n <10);
				}
			}
		}
	}
}
open Out,">$fOut";
print Out join("\n",@Head[0..9],"number_of_loci $nmarker",@Head[10..11]),"\n";
foreach my $Chr (sort keys %SORT) {
	my @Order=@{$SORT{$Chr}};
	foreach my $id (@Order) {
		my @out;
		for (my $i=0;$i<$nind;$i++) {
			push @out,$OGeno{"indi$i"}{$id};
		}
		print Out join("\t",$id,@out),"\n";
	}
}
close Out;
open Out,">$fOut.detail";
foreach my $Chr (sort keys %SORT) {
	my @Order=@{$SORT{$Chr}};
	foreach my $id (@Order) {
		my @out;
		for (my $i=0;$i<$nind;$i++) {
			if ($Geno{"indi$i"}{$id} eq $OGeno{"indi$i"}{$id}) {
				push @out,$OGeno{"indi$i"}{$id}
			}else{
				push @out,join("->",$Geno{"indi$i"}{$id},$OGeno{"indi$i"}{$id});
			}
		}
		print Out join(" ",$id,@out),"\n";
	}
}
close Out;
#######################################################################################
print "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

sub USAGE {#
	my $usage=<<"USAGE";
Description: 
Version:  $Script
Contact: long.huang

Usage:
  Options:
	-i	<file>	input file 
	-o	<file>	output file
	-m	<file>	intput map file
USAGE
	print $usage;
	exit;
}

#!/usr/bin/perl -w
use warnings;
use strict;

die "usage: perl $0 blast.m8.txt gi_tax.[prot|nucl].xls output\n" unless(@ARGV == 3);
my ($blastM8,$taxFile,$output) = @ARGV;

my (%exist,%gi2tax,%used);

open INB,$blastM8 or die "$!\n";
while(<INB>){
	chomp;
	my @temp = split;
	if($temp[1] =~ /gi\|(\d+)/){
		$exist{$1} = 1;
	}else{
		$exist{$temp[1]}=1;
	}
}
close INB;

open INT,$taxFile or die "$!\n";
while(<INT>){
	chomp;
	next if(/^#/);
	my @temp = split /\t/;
	$gi2tax{$temp[1]} = $temp[3] if(exists $exist{$temp[1]});
}
close INT;

open INB2,$blastM8 or die "$!\n";
open OUT,"> $output" or die "$!\n";
#print OUT "#Sequence\tTaxonomy\n";
my $last = "";
my %sam2gi;
while(<INB2>){
	chomp;
	next if(/^#/);
	my @temp = split;
        my $hit = join("\t",@temp[1..11]);
	$last = $temp[0] unless($last);
	#next if(exists $used{$temp[0]});
	my $gi = $temp[1];
        next if(exists $sam2gi{$temp[0]}{$gi});
	$sam2gi{$temp[0]}{$gi} = 1;
	if(exists $gi2tax{$gi}){
		#print OUT "$temp[0]\t$gi2tax{$gi}\n";
		print OUT "$temp[0]\t$gi2tax{$gi}\t$hit\n";
		$used{$temp[0]} = 1;
	}
	unless($last eq $temp[0]){
		print OUT "$last\tNO-TaxID\t$hit\n" unless(exists $used{$last});
	}	
	$last = $temp[0];
}
print OUT "$last\tNO-TaxID\n" unless(exists $used{$last});
close INB2;
close OUT;

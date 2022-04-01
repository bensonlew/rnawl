#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fIn,$fOut,$popt,$homo,$seg,$mis,$map);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Statistics::Distributions;
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"input:s"=>\$fIn,
	"output:s"=>\$fOut,
	"popt:s"=>\$popt,
	"homo:s"=>\$homo,
	"map:s"=>\$map,
			) or &USAGE;
&USAGE unless ($fIn and $fOut and $popt);
open In,$map;
my %mapped;
my %group;
my $groupid;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/ );
	if (/group/) {
		$groupid=(split(/\s+/,$_))[1];
	}
	my ($id,$dis)=split(/\s+/,$_);
	$mapped{$id}=$dis;
	$group{$id}=$groupid;
}
close In;
close In;
$homo||="aa";
$mis||=0.3;
$seg||=0.05;
open In,$fIn;
open Out,">$fOut";
print Out "#MarkerID\tGroupID\tDis\ttype\tNind\tNmiss\tGeno\tNGeno\tpvalue\tSegretion\n";
my %stat;
	my $t=0;
	my $n=1;

my %region;
while (<In>) {
	chomp;
	next if ($_ eq "" || /^$/ ||/=/||/^#/ || /MarkerID/ );
	last if(/^;/ || /individual names/);
	my ($id,$type,$phase,@info);
	if ($popt eq "CP") {
		 ($id,$type,$phase,@info) =split(/\t/,$_);
		 $type=~s/\>//g;
		 $type=~s/\<//g;
	}else{
		 s/\tX\t/\tab\t/g;
		 s/\tA\t/\taa\t/g;
		 s/\tB\t/\tbb\t/g;
		 s/\tU\t/\t--\t/g;
		 s/\t-\t/\t--\t/g;
		 ($id,@info)=split(/\t/,$_);
		 $type="aaxbb";
	}
	next if (!exists $mapped{$id});
	$stat{total}++;

	if ($type eq "aaxbb" && $popt eq "BC1") {
		for (my $i=0;$i<@info;$i++) {
			if ($info[$i] eq "bb") {
				$info[$i] = "--";
			}
		}
	}
	if ($type eq "abxcc") {
		$type = "lmxll";
		for (my $i=0;$i<@info;$i++) {
			if ($info[$i] eq "ac") {
				$info[$i]="ll"
			}elsif ($info[$i] eq "bc") {
				$info[$i]="lm"
			}elsif ($info[$i] eq "ab"||$info[$i] eq "cc") {
				$info[$i]="--"
			}
		}
	}elsif ($type eq "ccxab") {
		$type = "nnxnp";
		for (my $i=0;$i<@info;$i++) {
			if ($info[$i] eq "ac") {
				$info[$i]="np"
			}elsif ($info[$i] eq "bc") {
				$info[$i]="nn"
			}elsif ($info[$i] eq "ab"|| $info[$i] eq "cc") {
				$info[$i]="--"
			}
		}
	}
	my %gstat;
	my $miss=0;
	for (my $i=0;$i<@info;$i++){
		if ($info[$i] eq "ff" || $info[$i] eq "gg" ||$info[$i] eq "mm" || $info[$i] eq "pp") {
			$info[$i] = "--";
		}elsif (($info[$i] eq "ab"||$info[$i] eq "cd") && $type eq "abxcd") {
			$info[$i] = "--";
		}
		if ($info[$i] eq "--"){
			$miss++;
			next;
		}
		$gstat{$info[$i]}++;
	}
	my $flag;
	my $order;
	my($p)=&SegregationX2($type,\@info,\$flag,\$order);
	print Out "$id\t",$group{$id},"\t",$mapped{$id},"\t$type\t",scalar @info,"\t",$miss,"\t",$order,"\t",$flag,"\n";
	if ($p < 0.05) {
		my $LG=$group{$id};
		if (!exists $region{$LG}) {
			$n++;
			push @{$region{$LG}{$n}},$id;
			$t=0;
		}else{
			if ($t == 0) {
				push @{$region{$LG}{$n}},$id;
			}else{
				$n++;
				push @{$region{$LG}{$n}},$id;
				$t=0;
			}
		}
	}else{
		$t=1;
	}
}
close Out;
close In;
open Region,">$fOut.seg.region";
print Region "#LG\tPos1\tPos2\tMarkerNum\n";
foreach my $LG (sort keys %region) {
	foreach my $n (sort {$a<=>$b} keys %{$region{$LG}}) {
		my @m=@{$region{$LG}{$n}};
		print Region join("\t",$LG,$mapped{$m[0]},$mapped{$m[-1]},scalar @m),"\n";
	}
}
close Region;
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub SegregationX2{#.....
	my($type,$genotype,$flag,$order)=@_;
	my %Data;
	foreach  (@$genotype) {
		$Data{$_}++;
	}
	if($type eq "hkxhk"){
		my ($hh,$hk,$kk,$missingData);
		$hh=exists $Data{"hh"}?$Data{"hh"}:0;
		$hk=exists $Data{"hk"}?$Data{"hk"}:0;
		$kk=exists $Data{"kk"}?$Data{"kk"}:0;
		$missingData=exists $Data{"--"}?$Data{"--"}:0;
		my $all=$hh+$hk+$kk+$missingData;
		my $valid=$hh+$hk+$kk;
		my $genotypeOrder="hh:hk:kk";
		my $theoretical_segregation="1:2:1";
		my $segregation="$hh:$hk:$kk";
		$$order="$genotypeOrder\t$segregation";
		my $Segregation_p=Segregation($theoretical_segregation,$segregation,$valid);
		$$flag=" $Segregation_p\t$genotypeOrder=$segregation";
		return ($Segregation_p);
	}elsif($type eq "efxeg"){
		my ($eg,$ef,$ee,$fg,$missingData);
		$eg=exists $Data{"eg"}?$Data{"eg"}:0;
		$ef=exists $Data{"ef"}?$Data{"ef"}:0;
		$ee=exists $Data{"ee"}?$Data{"ee"}:0;
		$fg=exists $Data{"fg"}?$Data{"fg"}:0;
		$missingData=exists $Data{"--"}?$Data{"--"}:0;
		my $all=$eg+$ef+$ee+$fg+$missingData;
		my $valid=$eg+$ef+$ee+$fg;
		my $genotypeOrder="ee:ef:eg:fg";
		my $theoretical_segregation="1:1:1:1";
		my $segregation="$ee:$ef:$eg:$fg";
		$$order="$genotypeOrder\t$segregation";
		my $Segregation_p=Segregation($theoretical_segregation,$segregation,$valid);
		$$flag=" $Segregation_p\t$genotypeOrder=$segregation";
		return ($Segregation_p);
	}elsif($type eq "abxcd"){
		my ($ac,$ad,$bc,$bd,$missingData);
		$ac=exists $Data{"ac"}?$Data{"ac"}:0;
		$ad=exists $Data{"ad"}?$Data{"ad"}:0;
		$bc=exists $Data{"bc"}?$Data{"bc"}:0;
		$bd=exists $Data{"bd"}?$Data{"bd"}:0;
		$missingData=exists $Data{"--"}?$Data{"--"}:0;
		my $all=$ac+$ad+$bc+$bd+$missingData;
		my $valid=$ac+$ad+$bc+$bd;
		my $genotypeOrder="ac:ad:bc:bd";
		my $theoretical_segregation="1:1:1:1";
		my $segregation="$ac:$ad:$bc:$bd";
		$$order="$genotypeOrder\t$segregation";
		my $Segregation_p=Segregation($theoretical_segregation,$segregation,$valid);
		$$flag=" $Segregation_p\t$genotypeOrder=$segregation";
		return ($Segregation_p);
	}elsif($type eq "lmxll"){
		my ($ll,$lm,$missingData);
		$ll=exists $Data{"ll"}?$Data{"ll"}:0;
		$lm=exists $Data{"lm"}?$Data{"lm"}:0;
		$missingData=exists $Data{"--"}?$Data{"--"}:0;
		my $all=$ll+$lm+$missingData;
		my $valid=$ll+$lm;
		my $genotypeOrder="lm:ll";
		my $theoretical_segregation="1:1";
		my $segregation="$lm:$ll";
		$$order="$genotypeOrder\t$segregation";
		my $Segregation_p=Segregation($theoretical_segregation,$segregation,$valid);
		$$flag=" $Segregation_p\t$genotypeOrder=$segregation";
		return ($Segregation_p);
	}elsif($type eq "nnxnp"){
		my ($nn,$np,$missingData);
		$nn=exists $Data{"nn"}?$Data{"nn"}:0;
		$np=exists $Data{"np"}?$Data{"np"}:0;
		$missingData=exists $Data{"--"}?$Data{"--"}:0;
		my $all=$nn+$np+$missingData;
		my $valid=$nn+$np;
		my $genotypeOrder="nn:np";
		my $theoretical_segregation="1:1";
		my $segregation="$nn:$np";
		$$order="$genotypeOrder\t$segregation";
		my $Segregation_p=Segregation($theoretical_segregation,$segregation,$valid);
		$$flag=" $Segregation_p\t$genotypeOrder=$segregation";
		return ($Segregation_p);
	}elsif($type eq "abxcc"){
		my ($ac,$bc,$missingData);
		$ac=exists $Data{"ac"}?$Data{"ac"}:0;
		$bc=exists $Data{"bc"}?$Data{"bc"}:0;
		$missingData=exists $Data{"--"}?$Data{"--"}:0;
		my $all=$ac+$bc+$missingData;
		my $valid=$ac+$bc;
		my $genotypeOrder="ac:bc";
		my $theoretical_segregation="1:1";
		my $segregation="$ac:$bc";
		$$order="$genotypeOrder\t$segregation";
		my $Segregation_p=Segregation($theoretical_segregation,$segregation,$valid);
		$$flag=" $Segregation_p\t$genotypeOrder=$segregation";
		return ($Segregation_p);
	}elsif($type eq "ccxab"){
		my ($ac,$bc,$missingData);
		$ac=exists $Data{"ac"}?$Data{"ac"}:0;
		$bc=exists $Data{"bc"}?$Data{"bc"}:0;
		$missingData=exists $Data{"--"}?$Data{"--"}:0;
		my $all=$ac+$bc+$missingData;
		my $valid=$ac+$bc;
		my $genotypeOrder="ac:bc";
		my $theoretical_segregation="1:1";
		my $segregation="$ac:$bc";
		$$order="$genotypeOrder\t$segregation";
		my $Segregation_p=Segregation($theoretical_segregation,$segregation,$valid);
		$$flag=" $Segregation_p\t$genotypeOrder=$segregation";
		return ($Segregation_p);
	}elsif($type eq "aaxbb"){
		my ($aa,$bb,$ab,$missingData);
		$aa=exists $Data{"aa"}?$Data{"aa"}:0;
		$bb=exists $Data{"bb"}?$Data{"bb"}:0;
		$ab=exists $Data{"ab"}?$Data{"ab"}:0;
		$missingData=exists $Data{"--"}?$Data{"--"}:0;
		my $all=$aa+$bb+$ab+$missingData;
		my $valid=$aa+$bb+$ab;
		my $genotypeOrder="ab:aa:bb";
		if($popt=~/BC\d+/i){
			$genotypeOrder=$homo=~/aa/i?"aa:ab":"bb:ab";
		}elsif ($popt=~/Ri\d+/i || $popt=~/DH/i) {
			$genotypeOrder="aa:bb";
		}
		my $theoretical_segregation;
		if($popt=~/f2/i || $popt=~/F2/i){
			$theoretical_segregation="2:1:1";
		}elsif ($popt=~/Ri\d+/i || $popt=~/DH/i) {
			$theoretical_segregation="1:1";
		}elsif($popt=~/BC\d+/i){
			$theoretical_segregation="1:1"; 
		}elsif ($popt = "CP") {
			my $segregation="$ab:$aa:$bb";
			$$flag="$segregation";
			$$order="$genotypeOrder\t$segregation";
			return("-");
		}else{
			warn "please input true group!\n";
			exit(0);
		}
		my $segregation="$ab:$aa:$bb";
		if($popt=~/bc\d+/i){
			$segregation=$homo=~/aa/i?"$aa:$ab":"$bb:$ab";
		}elsif ($popt=~/Ri\d+/i || $popt=~/DH/i) {
			$segregation="$aa:$bb";
		}
		$$order="$genotypeOrder\t$segregation";
		my $Segregation_p=Segregation($theoretical_segregation,$segregation,$valid);
		$$flag=" $Segregation_p\t$genotypeOrder=$segregation";
		return ($Segregation_p);
	}else{#.......
		warn "wrong grouptype! $type\n";
		exit(0);
	}

	sub Segregation {#
		my ($theoretical_segregation,$segregation,$all)=@_;
		my @a=split ":",$theoretical_segregation;
		my @b=split ":",$segregation;
		return "0.01" if (scalar @a != scalar @b || $all == 0) ;
		my @theoretical;
		my $a_sum=0;
		$a_sum+=$_ foreach (@a);
		push @theoretical,$_/$a_sum*$all foreach (@a);
		my $df=scalar @a -1;
		my $X2=0;
		if ($df == 1) {
			for (my $i=0;$i<@a ;$i++) {
				$X2+=X2df2($b[$i],$theoretical[$i]);
			}
		}else{
			for (my $i=0;$i<@a ;$i++) {
				$X2+=X2df1($b[$i],$theoretical[$i]);
			}
		}
		my $p=0;
		$p=Statistics::Distributions::chisqrprob($df,$X2);
		return int($p*10000)/10000;
	}

	sub X2df1 {#
		my ($A,$T)=@_;
		return ($A-$T)**2/$T;
	}

	sub X2df2 {#
		my ($A,$T)=@_;
		return (abs($A-$T)-0.5)**2/$T;
	}
}

sub USAGE {#
        my $usage=<<"USAGE";
Contact:        long.huang\@majorbio.com;
Script:			$Script
Description:
	eg:
	perl $Script -i -o 

Usage:
  Options:
  -input	<file>	input file name
  -output	<file>	input keys of file name
  -popt	<str>	population type CP/BCi/Fi/Rix/
  -homo	<str>	if BC homotype,default=aa
  -seg	<num>	p value of segragation,default 0.05
  -map	<file>	map file
  -h         Help

USAGE
        print $usage;
        exit;
}

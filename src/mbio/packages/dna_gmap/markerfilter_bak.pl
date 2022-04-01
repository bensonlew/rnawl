#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fIn,$fOut,$fPos,$Seg,$Miss,$minPMdep,$minOdep,$maxPMdep,$maxOdep,$Pop,$Vtype,$newmarker,$indilist,$mmarker,$mindi);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Statistics::Distributions;

my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"info:s"=>\$fIn,
	"out:s"=>\$fOut,
	"pop:s"=>\$Pop,
	"seg:s"=>\$Seg,
	"mis:s"=>\$Miss,
	"minPMdep:s"=>\$minPMdep,
	"maxPMdep:s"=>\$maxPMdep,
	"minOdep:s"=>\$minOdep,
	"maxOdep:s"=>\$maxOdep,
	"vtype:s"=>\$Vtype,
	"newmarker:s"=>\$newmarker,
	"indilist:s"=>\$indilist,
			) or &USAGE;
&USAGE unless ($fIn and $fOut and $Pop);
my $homo="aa";
$Seg||=0.05;
$Miss||=0.3;
$Pop||="F2";
$minPMdep||=10;
$minOdep||=5;
$maxPMdep||=1000000;
$maxOdep||=1000000;
$minPMdep=~s/_//;
$minOdep=~s/_//;
$maxPMdep=~s/_//;
$maxOdep=~s/_//;
$Vtype||="ALL";
$indilist||="";
my %finalindi;
my $count=1;
	open In,$fIn;
	if ($fIn =~ /.gz/) {
		close In;
		open In,"zcat $fIn|";
	}
	while (<In>) {
		chomp;
		next if ($_ eq "" || /^$/ || /^##/);
		if (/^#/) {
			my (undef,undef,undef,undef,undef,@indi)=split(/\t/,$_);
			for (my $i=0;$i<@indi;$i++) {
				$finalindi{$indi[$i]}++;
			}
			last;
		}
	}
	close In;

if ($indilist ne "") {
	my @indilist=split(/,/,$indilist);
	# my %finalindi;
	foreach my $indi(@indilist){
		$finalindi{$indi}++
	}
	$count++;
}
if (defined $newmarker) {
		$count++;
		open In,$newmarker;
		if ($newmarker =~ /.gz/) {
			close In;
			open In,"zcat $newmarker|";
		}
		while (<In>) {
			chomp;
			next if ($_ eq "" || /^$/ || /^##/);
			if (/^#/) {
				my (undef,undef,@indi)=split(/\t/,$_);
				for (my $i=0;$i<@indi;$i++) {
					$finalindi{$indi[$i]}++;
				}
			}
		}
		close In;
}


open In,$fIn;
if ($fIn =~ /.gz/) {
	close In;
	open In,"zcat $fIn|";
}
open Out,">$fOut.detail.info";
open Marker,">$fOut.marker";
open MarkerStat,">$fOut.markerStat";
my @indi;
my $n=0;
my %stat;
my %Mdep;
my %Idep;
my %Alldep;
my %Type;
while (<In>) {
	chomp;
	next if ($_ eq "" || /^$/ || /^##/);
	if (/^#/) {
		(undef,undef,undef,undef,undef,@indi)=split(/\t/,$_);
		my @outid;
		for (my $i=0;$i<@indi;$i++) {
			next if (!exists $finalindi{$indi[$i]} || $finalindi{$indi[$i]} != $count);
			push @outid,$indi[$i];
		}
		print Out join("\t","#MarkerID","Gtype","Vtype","Pdep","Mdep",@outid),"\n";
		print Marker join("\t","#MarkerID","type",@outid),"\n";
		print MarkerStat "#MarkerID\ttype\tNind\tNmiss\tGeno\tNGeno\tSegretion\n";

	}else{
		my ($markerid,$gtype,$vtype,$pdep,$mdep,@geno)=split(/\t/,$_);
		next if ($Pop eq "CP" && $gtype eq "aaxbb"); #F1
		next if ($Pop ne "CP" && $gtype ne "aaxbb"); #F2
		next if ($pdep ne "--" && ($pdep < $minPMdep || $mdep < $minPMdep)); #亲本深度
		next if ($pdep ne "--" && ($pdep > $maxPMdep || $mdep > $maxPMdep));
		next if ($Vtype ne "ALL" && $Vtype ne $vtype); #变异类型
		my $tdep=0;
		my @info;
		my $miss=0;
		my @Geno;
		for (my $i=0;$i<@geno;$i++) {
			next if (!exists $finalindi{$indi[$i]} || $finalindi{$indi[$i]} != $count);
			my ($dep,$geno)=split(/\,/,$geno[$i]);
			if ($Pop eq "DH" && $geno eq "ab") {
				$geno="--";
				$dep=0;
				$geno[$i]=join(",","0",$geno);
			}elsif ($gtype eq "nnxnp" && $geno eq "pp") {
				$geno="--";
				$dep=0;
				$geno[$i]=join(",","0",$geno);
			}elsif ($gtype eq "lmxll" && $geno eq "mm") {
				$geno="--";
				$dep=0;
				$geno[$i]=join(",","0",$geno);
			}elsif ($gtype eq "abxcc") {
				$gtype="lmxll";
				if ($geno[$i] eq "ac") {
					$geno[$i]="ll"
				}elsif ($geno[$i] eq "bc") {
					$geno[$i]="lm"
				}else {
					$dep=0;
					$geno[$i]="--"
				}
				$geno[$i]=join(",","0",$geno);
			}elsif ($gtype eq "ccxab"){
				$gtype="nnxnp";
				if ($geno[$i] eq "ac") {
					$geno[$i]="np"
				}elsif ($geno[$i] eq "bc") {
					$geno[$i]="nn"
				}else {
					$dep=0;
					$geno[$i]="--"
				}
				$geno[$i]=join(",","0",$geno);
			}
			$tdep+=$dep if ($dep ne "--");;
			$miss++ if ($geno eq "--");
			push @info,$geno;
			push @Geno,$geno[$i];
		}
		my $flag;
		my $order;
		my($p)=&SegregationX2($gtype,\@info,\$flag,\$order);
		next if ($p < $Seg);#偏分离
		next if ($miss/scalar @info > $Miss);#缺失
		next if ($tdep!=0 && $tdep / scalar @info < $minOdep);#子代深度
		next if ($tdep!=0 && $tdep / scalar @info > $maxOdep);#子代深度
		$Type{$gtype}{$vtype}++;
		$Mdep{$markerid}=$tdep / scalar @info > $maxOdep;
		for (my $i=0;$i<@geno;$i++) {
			next if (!exists $finalindi{$indi[$i]} || $finalindi{$indi[$i]} != $count);
			my ($dep,$geno)=split(/\,/,$geno[$i]);
			my $gdep=$dep;
			if ($dep eq "--" && $geno ne "--") {
				$gdep=$minOdep;				
			}elsif ($geno eq "--") {
				$gdep=0;
			}
			$Idep{$indi[$i]}{dep}+=$gdep if ($geno ne "--");
			$Idep{$indi[$i]}{num}++ if ($geno ne "--");
			$Idep{$indi[$i]}{miss}++ if ($geno eq "--");
			$Alldep{$markerid}{$indi[$i]}=$gdep;
		}
		print Marker join("\t",$markerid,$gtype,@info),"\n";
		print Out join("\t",$markerid,$gtype,$vtype,$pdep,$mdep,@Geno),"\n";
		print MarkerStat "$markerid\t$gtype\t",scalar @info,"\t",$miss,"\t",$order,"\t",$flag,"\n";

	}
}
close In;
if (defined $newmarker) {
	open In,$newmarker;
	while (<In>) {
		chomp;
		next if ($_ eq "" || /^$/ || /^##/);
		if (/^#/) {
			(undef,undef,@indi)=split(/\t/,$_);
			my @outid;
			for (my $i=0;$i<@indi;$i++) {
				next if (!exists $finalindi{$indi[$i]} || $finalindi{$indi[$i]} != $count);
				push @outid,$indi[$i];
			}
		}else{
			my ($markerid,$gtype,@geno)=split(/\t/,$_);
			next if ($Pop eq "CP" && $gtype eq "aaxbb"); #F1
			next if ($Pop ne "CP" && $gtype ne "aaxbb"); #F2
			my @info;
			my $miss=0;
			my @Geno;
			for (my $i=0;$i<@geno;$i++) {
				next if (!exists $finalindi{$indi[$i]} || $finalindi{$indi[$i]} != $count);
				my $geno=$geno[$i];
				my $dep="--";
				if ($Pop eq "DH" && $geno eq "ab") {
					$geno="--";
				}elsif ($gtype eq "nnxnp" && $geno eq "pp") {
					$geno="--";
					$dep=0;
				}elsif ($gtype eq "lmxll" && $geno eq "mm") {
					$geno="--";
				}elsif ($gtype eq "abxcc") {
					$gtype="lmxll";
					if ($geno[$i] eq "ac") {
						$geno[$i]="ll"
					}elsif ($geno[$i] eq "bc") {
						$geno[$i]="lm"
					}else {
						$geno[$i]="--"
					}
				}elsif ($gtype eq "ccxab"){
					$gtype="nnxnp";
					if ($geno[$i] eq "ac") {
						$geno[$i]="np"
					}elsif ($geno[$i] eq "bc") {
						$geno[$i]="nn"
					}else {
						$geno[$i]="--"
					}
				}
				$geno[$i]=join(",",$dep,$geno);
				$miss++ if ($geno eq "--");
				push @info,$geno;
				push @Geno,$geno[$i];
			}
			my $flag;
			my $order;
			my($p)=&SegregationX2($gtype,\@info,\$flag,\$order);
			#next if ($p < $Seg);#偏分离
			next if ($miss/scalar @info > $Miss);#缺失
			$Type{$gtype}{upload}++;
			$Mdep{$markerid}=$minOdep;
			for (my $i=0;$i<@geno;$i++) {
				next if (!exists $finalindi{$indi[$i]} || $finalindi{$indi[$i]} != $count);
				my ($dep,$geno)=split(/\,/,$geno[$i]);
				$Idep{$indi[$i]}{dep}+=$minOdep if ($geno ne "--");
				$Idep{$indi[$i]}{num}++ if ($geno ne "--");
				$Idep{$indi[$i]}{miss}++ if ($geno eq "--");
				$Alldep{$markerid}{$indi[$i]}=$minOdep;
			}
			print Marker join("\t",$markerid,$gtype,@info),"\n";
			print Out join("\t",$markerid,$gtype,"--","--","--",@Geno),"\n";
			print MarkerStat "$markerid\t$gtype\t",scalar @info,"\t",$miss,"\t",$order,"\t",$flag,"\n";
		}
	}
	close In;
}
close Out;
close Marker;
close MarkerStat;
my %Adep;
foreach my $indi (sort keys %Idep) {
	$Adep{$indi}=$Idep{$indi}{dep}/$Idep{$indi}{num};
}
open Matrix,">$fOut.matrix";
print Matrix join("\t","#MarkerID",sort {$Adep{$b}<=>$Adep{$a}}keys %Adep),"\n";
foreach my $markerid (sort {$Mdep{$b}<=>$Mdep{$a}} keys %Mdep) {
	my @out;
	push @out,$markerid;
	foreach my $indiid (sort {$Adep{$b}<=>$Adep{$a}}keys %Adep) {
		push @out,$Alldep{$markerid}{$indiid};
	}
	print Matrix join("\t",@out),"\n";
}
close Matrix;
open Stat,">$fOut.indi.stat";
print Stat "#Indi\tAverageDepth\tMissRatio\n";
foreach my $indi(sort keys %Adep){
	$Idep{$indi}{num}||=0;
	$Idep{$indi}{miss}||=0;
	print Stat join("\t",$indi,$Adep{$indi},$Idep{$indi}{miss}/($Idep{$indi}{num}+$Idep{$indi}{miss})),"\n";
}
close Stat;
open Out,">$fOut.gtype.stat";
foreach my $type  (sort keys %Type) {
	$Type{$type}{SNP}||=0;
	$Type{$type}{INDEL}||=0;
	print Out $type,"\t",$Type{$type}{SNP},"\t",$Type{$type}{INDEL},"\n";
}
close Out;
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
		if($Pop=~/BC\d+/i){
			$genotypeOrder=$homo=~/aa/i?"aa:ab":"bb:ab";
		}elsif ($Pop=~/Ri\d+/i || $Pop=~/DH/i) {
			$genotypeOrder="aa:bb";
		}
		my $theoretical_segregation;
		if($Pop=~/f2/i || $Pop=~/F2/i){
			$theoretical_segregation="2:1:1";
		}elsif ($Pop=~/Ri\d+/i || $Pop=~/DH/i) {
			$theoretical_segregation="1:1";
		}elsif($Pop=~/BC\d+/i){
			$theoretical_segregation="1:1"; 
		}elsif ($Pop = "CP") {
			my $segregation="$ab:$aa:$bb";
			$$flag="$segregation";
			$$order="$genotypeOrder\t$segregation";
			return("-");
		}else{
			warn "please input true group!\n";
			exit(0);
		}
		my $segregation="$ab:$aa:$bb";
		if($Pop=~/bc\d+/i){
			$segregation=$homo=~/aa/i?"$aa:$ab":"$bb:$ab";
		}elsif ($Pop=~/Ri\d+/i || $Pop=~/DH/i) {
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

Usage:
  Options:
	"info:s"=>\$fIn,   #marker detail
	"out:s"=>\$fOut,   #output file prefix
	"pop:s"=>\$Pop,    #population type
	"seg:s"=>\$Seg,    #segmental 0.05
	"mis:s"=>\$Miss,   # miss 0.3
	"minPMdep:s"=>\$minPMdep, #minPMdep 10
	"maxPMdep:s"=>\$maxPMdep, #maxPMdep 100000
	"minOdep:s"=>\$minOdep, #minOdep 5
	"maxOdep:s"=>\$maxOdep, #maxOdep 100000
	"vtype:s"=>\$Vtype, #Vtype SNP/INDEL/ALL
	"newmarker:s"=>\$newmarker, #custom marker file
	"indilist:s"=>\$indilist, # individual list
  -h         Help

USAGE
        print $usage;
        exit;
}

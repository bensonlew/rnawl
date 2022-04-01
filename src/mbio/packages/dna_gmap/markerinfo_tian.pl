#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fIn,$mark,$fOut,$popt);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Statistics::Distributions;
my $version="1.0.0";
GetOptions(
    "help|?" =>\&USAGE,
    "vcf:s"=>\$fIn,
    "mark:s"=>\$mark,
    "out:s"=>\$fOut,
    "popt:s"=>\$popt,
            ) or &USAGE;
&USAGE unless ($fIn and $mark and $fOut);
$popt||="F2";
my (%stat,%region,%info);
my (@sample,@marker);
my $homo||="aa";
open In,$mark;
my %sample;
while (<In>){
    next if ($_ eq ""|| /^$/);
    chomp;
    if($_=~/^#/){
        (undef,undef,@sample)=split(/\t/,$_);
        foreach my $sample (@sample) {
            $sample{$sample}=1;
        }
    }else{
        #my($id,@undi)=split(/\t/,$_);
        my($id,$type,@info)=split(/\t/,$_);
        #$info{$id}{$type}=join("\t",@info);
        push @marker,$id;
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
                    $info[$i]="ll";
                }elsif ($info[$i] eq "bc") {
                    $info[$i]="lm";
                }elsif($info[$i] eq "ab"||$info[$i] eq "cc") {
                    $info[$i]="--";
                }
            }
        }elsif($type eq "ccxab") {
            $type = "nnxnp";
            for (my $i=0;$i<@info;$i++) {
                if ($info[$i] eq "ac") {
                    $info[$i]="np";
                }elsif($info[$i] eq "bc"){
                    $info[$i]="nn";
                }elsif($info[$i] eq "ab" || $info[$i] eq "cc"){
                    $info[$i]="--";
                }
            }
        }
        my $miss=0;
        for (my $i=0;$i<@info;$i++){
            if ($info[$i] eq "ff" || $info[$i] eq "gg" ||$info[$i] eq "mm" || $info[$i] eq "pp") {
                $info[$i] = "--";
            }elsif (($info[$i] eq "ab"||$info[$i] eq "cd") && $type eq "abxcd") {
                $info[$i] = "--";
            }
            if ($info[$i] eq "--"){
                $miss++ ;
                next;
            }
        }
        my $flag;
        my $order;
        my $X2||=0;
        my $df||=0;
        $miss=sprintf("%.2f",$miss/(scalar@sample));
        #my ($p)=&SegregationX2($type,\@info,\$X2,\$df,\$flag,\$order);
        my @result=&SegregationX2($type,\@info,\$X2,\$df,\$flag,\$order);
        $info{$id}{type}=join("\t",$miss,@result);
        
    }
}
close In;
my (%dep);
open In,$fIn;
if($fIn=~/gz$/){
    close In;
    open In,"gunzip -c $fIn|";
}
my @Indi;
my %vtype;
while (<In>) {
    chomp;
    next if ($_ eq ""|| /^$/ || /^##/);
    if (/^#/) {
        my ($chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$format,@indi)=split(/\t/,$_);
        @Indi=@indi;
    }else{
        my ($chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$format,@samples)=split(/\t/,$_);
        next if (!exists $info{$id});
        my @format=split(/\:/,$format);
        my $markerdep=0;
        my $marker=$id; 
        my @alles=split(",",join(",",$ref,$alt));
        my %len;
        foreach my $alles(@alles){
            $len{length($alles)}++;
        }
        my $vtype="SNP";
        $vtype="INDEL" if(scalar keys %len >1);
        $vtype{$id}=$vtype;
        for (my $i=0;$i<@Indi;$i++) {
                    my @info=split(/\:/,$samples[$i]);
                    my $samples=$Indi[$i];
                    next if (!exists $sample{$samples});
                            for (my $j=0;$j<@info;$j++) {
                                if ($format[$j] eq "DP") {
                                    #$per{$sample}{$marker}=$info[$j];
                                    #$stat{$marker}{$sample}{dep}=$info[$j];
                                    $markerdep=$markerdep + $info[$j];
                                }
                            }

                }
                $markerdep=$markerdep/(scalar @sample);
                $stat{$marker}{depth}=$markerdep;
    }
}
close In;

#$miss,$fiag,$order;
open Out,">$fOut";
print Out "#Locus\tChr\tPosition\tVariant Type\tAverage Depth\tMiss Ratio\tX2\tDf\tSignif\tClassification\n";
    #my($miss,$fiag,$order)=split/\t/,$info{$id}{type};
    foreach my $marker(sort keys%stat){
            my($chr,$pos)=split(/\_/,$marker);
            #print Out join("\t",$marker,$chr,$pos,sprintf("%.2f",$stat{$marker}{depth}),$miss,$fiag,$order),"\n";
            print Out join("\t",$marker,$chr,$pos,$vtype{$marker},sprintf("%.2f",$stat{$marker}{depth}),$info{$marker}{type}),"\n";
        }

close Out;


#######################################################################################
#print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub SegregationX2{#.....
        my($type,$genotype,$X2,$df,$flag,$order)=@_;
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
                #my $Segregation_p=Segregation($theoretical_segregation,$segregation,$valid);
                #$$flag="$Segregation_p\t$genotypeOrder=$segregation";
                #return ($Segregation_p);
                my @result=Segregation($theoretical_segregation,$segregation,$valid);
                $$flag="$genotypeOrder=$segregation";
                push@result,$$flag;
                return @result;
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
                #my $Segregation_p=Segregation($theoretical_segregation,$segregation,$valid);
                #$$flag=" $Segregation_p\t$genotypeOrder=$segregation";
                #return ($Segregation_p);
                my @result=Segregation($theoretical_segregation,$segregation,$valid);
                $$flag="$genotypeOrder=$segregation";
                push@result,$$flag;
                return @result;
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
                #my $Segregation_p=Segregation($theoretical_segregation,$segregation,$valid);
                #$$flag=" $Segregation_p\t$genotypeOrder=$segregation";
                #return ($Segregation_p);
                my @result=Segregation($theoretical_segregation,$segregation,$valid);
                $$flag="$genotypeOrder=$segregation";
                push@result,$$flag;
                return @result;
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
                #my $Segregation_p=Segregation($theoretical_segregation,$segregation,$valid);
                #$$flag=" $Segregation_p\t$genotypeOrder=$segregation";
                #return ($Segregation_p);
                my @result=Segregation($theoretical_segregation,$segregation,$valid);
                $$flag="$genotypeOrder=$segregation";
                push@result,$$flag;
                return @result;
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
                #my $Segregation_p=Segregation($theoretical_segregation,$segregation,$valid);
                #$$flag=" $Segregation_p\t$genotypeOrder=$segregation";
                #return ($Segregation_p);
                my @result=Segregation($theoretical_segregation,$segregation,$valid);
                $$flag="$genotypeOrder=$segregation";
                push@result,$$flag;
                return @result;
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
                #my $Segregation_p=Segregation($theoretical_segregation,$segregation,$valid);
                #$$flag=" $Segregation_p\t$genotypeOrder=$segregation";
                #return ($Segregation_p);
                my @result=Segregation($theoretical_segregation,$segregation,$valid);
                $$flag="$genotypeOrder=$segregation";
                push@result,$$flag;
                return @result;
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
                #my $Segregation_p=Segregation($theoretical_segregation,$segregation,$valid);
                #$$flag=" $Segregation_p\t$genotypeOrder=$segregation";
                #return ($Segregation_p);
                my @result=Segregation($theoretical_segregation,$segregation,$valid);
                $$flag="$genotypeOrder=$segregation";
                push@result,$$flag;
                return @result;
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
                my @result=Segregation($theoretical_segregation,$segregation,$valid);
                #$$flag="$Segregation_p\t$genotypeOrder=$segregation";
                $$flag="$genotypeOrder=$segregation";
                push@result,$$flag;
                return @result;
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
                #&Statistics::Distributions::chisqrprob
                $p=Statistics::Distributions::chisqrprob($df,$X2);
                #return int($p*10000/10000);
                $p=sprintf("%.4f",$p);
                $X2=sprintf("%.2f",$X2);
                my @result;
                push @result,join("\t",$X2,$df,$p);
                return @result;
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
Contact:        chongqing.shi\@majorbio.com;
Script:         $Script
Description:
    fq thanslate to fa format
    eg:
    perl $Script -i -o -k -c

Usage:
  Options:
  -vcf  <file>  input pop.vcf file
  -mark <file>  input filtered.marker
  -out  <file>  output result file
  -popt <stri>  pop type
  -h         Help

USAGE
        print $usage;
        exit;
}

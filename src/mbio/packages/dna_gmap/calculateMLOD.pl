#!/usr/bin/perl -w
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
my ($fIn,$fOut,$only,$simple);
GetOptions(
                "help|?" =>\&USAGE,
                "i:s"=>\$fIn,
                "o:s"=>\$fOut,
                "s"=>\$simple,
                "p:s"=>\$only,
                ) or &USAGE;
&USAGE unless ($fIn and $fOut);

my $dirname = dirname($fOut);
mkdir $dirname unless (-d $dirname) ;


#------------------------------------------------------------------
# global value  
#------------------------------------------------------------------

my %all = ();  ## genotype hash 

#------------------------------------------------------------------
# get data 
#------------------------------------------------------------------
#$debug=1;##<
open (IN,$fIn) or die "can not open input file!";
$/="\n";
my $k=0;
##>readline <IN>;##<
while (<IN>) {
    chomp;
    next if (/Type/ || /type/ || /^;/ || /^$/ || /individual/ || /Individual/ || /^#/ || /MarkerID/);##<  
    my @in=split /\s+/,$_;
    if ($in[1]=~/lmxll/ | /nnxnp/){  
        $all{$k}{"genotype_num"}=1;
        $all{$k}{"marker"}=$_;
    }elsif($in[1]=~/efxeg/ | /abxcd/){
        $all{$k}{"genotype_num"}=3;
        $all{$k}{"marker"}=$_;
    }else{
        $all{$k}{"genotype_num"}=2;
        $all{$k}{"marker"}=$_;
    }
    $k++;
}
$only||=$k-1;
close IN;
#------------------------------------------------------------------
# calculate mLOD and output 
#------------------------------------------------------------------

################ Method introduced from joinmap official website : http://www.kyazma.nl/index.php/mc.JoinMap/sc.FAQ/#Q9
#   For each pair of loci a contingency table is produced of the genotypes. The dimensions of the table depend on the population type and the segregation types of the loci (unknown genotypes are ignored). From this table the G statistic (-2 times the logarithm of the likelihood ratio for the Poisson distribution, see e.g. Fienberg, 1979, The analysis of cross-classified categorical data, MIT Press) is calculated to test for independence; the expected number (E) in each cell is calculated from the R-total (R), the column-total (C) and the grand-total (T): 
#
#   E = R*C/T .
#
#
#   The G statistic then is a summation (SUM) over all cells (O is the observed number, ln() is the natural logarithm): 
#
#   G = 2 * SUM [ O*ln(O/E) ] .
#
#
#   The G statistic (Gd) has an approximate chi-square distribution with the number of Rs in the table minus 1 multiplied by the number of columns minus 1 as the degrees of freedom (d). When the loci have different numbers of genotypes in their segregation, this would present a problem, because this number affects the degrees of freedom in the G test and in the testing of linkage one would need to take account of the degrees of freedom. In order to remove this problem and to ensure the comparability of data coming from different population or segregation types the G statistic with d degrees of freedom, Gd, is transformed approximately to a G statistic, G1, that would have been obtained if there was just a single degree of freedom (as if in a backcross). This approximate transformation is an empirically determined formula (exp() is the exponential function): 
#
#   e = exp( -Gd/(2*(d-1)) ) ,
#
#   G1 = ((4-e)*e - 3)*(d-1) + Gd . 
#
#
#   Because in genetics one is used to LOD scores, which are likelihood ratio statistics using the 10-base logarithm instead of the natural logarithm multiplied by -2, the modified LOD score (mLOD) is simply derived from G1: 
#
#   mLOD = G1 / (2*ln(10)) .
#
#   It can be shown that for the case of two loci each segregating in two genotypes (e.g. population types BC1 or DH1) when there is no segregation distortion this mLOD equals the 'normal' LOD score. The modified LOD score is not sensitive to segregation distortion, in contrast to the 'normal' LOD score.
#
#####################################################

open (OUT,">$fOut") or die "open out file error";
print OUT "#Loci_1\tLoci_2\tMLOD\n";
for (my $i=0; $i<$only; $i++) {
    for (my $j=$i+1; $j<$k; $j++) {
        my (%pairs,%R,%C);
        my($marker1,$genotype1,@indi1)=split /\s+/,$all{$i}{"marker"};
        my($marker2,$genotype2,@indi2)=split /\s+/,$all{$j}{"marker"};
        # die "different number of individuals between $marker1 and $marker2\n" if (@indi1 != @indi2) ;
        my $T = 0;
        for (my $p=0;$p<@indi2 ;$p++) {
            if ($indi1[$p] eq "--" or $indi2[$p] eq "--") { 
                next;
            }else{
                $R{"$indi1[$p]"}++;
                $C{"$indi2[$p]"}++;
                $pairs{"$indi1[$p]x$indi2[$p]"}++;
                $T++;
            }

        }
        my ($bG,$G);
        foreach my $key1 (keys %pairs) {
            my ($r,$c)=split /x/, $key1;
            my $Rv=$R{"$r"};##<
            my $Cv=$C{"$c"};##<
            my $E=$Rv*$Cv/$T;##<
            $bG+=$pairs{$key1}*log(($pairs{$key1})/($E)) if $E!=0;##<
        }
        if (!defined $bG) {
            $G=0;
        }else{
            $G=2*$bG;
        }
        my ($G1,$e);
        if (scalar keys %R < 2 or scalar keys %C < 2) {
            print "$marker1 and $marker2,lack of freedom to calculate mLOD\n";
            next;
        }
        if (scalar keys %R == 2 and scalar keys %C == 2) {
             $G1=$G;
        }else{
            my $df = (scalar(keys %R) -1) * (scalar(keys %C) -1);
            $e = exp(-$G /(2 * ($df -1)));
            $G1 = ((4-$e)*$e-3)*($df-1)+$G;
        }
        my $mlod=$G1/(2*log(10));
        $mlod=int($mlod) if($simple);
        print OUT"$marker1\t$marker2\t$mlod\n";
    }
}
close OUT;
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------

sub GetTime {
    my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
    return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}


sub USAGE {#
    my $usage=<<"USAGE";
Program:$0
Version: $version
Contact: huangl <long.huang\@majorbio.com>
Description:

    calculate modified LOD score of pairwise loci.
    
    perl $0 -i xxx.genotype -o XXX.mLOD
        
Usage:
  Options:
  -i <file> input file,genotype file,forced 
  -o <fiel> output file,hold result of mlod
  -s    input   simple file each mlod is integer to free malloc
  -h         Help

USAGE
    print $usage;
    exit;
}

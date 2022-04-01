#!/usr/bin/perl -w

use strict;
use Getopt::Long;
my ($windl,$step,$outdir,$help);
GetOptions(
        "windl:i"=>\$windl,
	"step:i"=>\$step,
        "outdir:s"=>\$outdir,
        "help"=>\$help
);
($help || @ARGV!=2) && die"Name: gc_depth.pl
Description: to draw gc_depth distribution figure with soap.coverage result
Usage: perl gc_depth.pl <in.fasta> <soap.coverf>
  <in.fasta>           input genomics file in fasta form
  <soap.coverf>        soap.coverage -depthsingle output result
  -windl  <num>        the length of window, default=10000
  -step   <num>        the length of step, default=500
  -outdir <dir>        output dir, default ./
  -help                output help information to screen

Example:
 nohup perl gc_depth.pl arowana.fa arowana.coverage.detail -outdir GC_depth -windl 5000 &
 
Notic:
 1 -x(y)_mun set can be ' ', then the scale will set accroding to the real data.
 2 draw.sh will creat at current dir, you can change some optins and run it to rewrite the figure.
 3 800M genomics use RES 800m and running time 15min
 4 results at outdir:
   gc_distribution, depth_distribution, gc_depth.wind\n";
## check error  ##
foreach(@ARGV){(-s $_) || die"Error: can't find file $_\n";}
my ($fasta,$coverf) = @ARGV;
## get options values  ##
$outdir ||= '.';
(-d $outdir) || mkdir"$outdir";
$windl ||= 10000;
$step ||= 500;
## Main process for stat ##
my (%lenh,%depth);
gc_depth($coverf,$windl,\%lenh,\%depth,1);
my $outf = "$outdir/gc_depth";
open GC,">$outf.wind";
select GC;
gc_depth($fasta,$windl,\%lenh,\%depth);
close GC;
(-s "$outf.wind") || die"Note: have not result, maybe you can decrease -windl $windl\n";
foreach(keys %depth){delete $depth{$_};}
foreach(keys %lenh){delete $lenh{$_};}
#`rm draw.r`;
## get gc and depth distribution file ##
my $get_distribution = 'awk \'BEGIN{ARGC=2}(x){if(n==$1){c=c" "$4;d=d" "$6}else{print n"\n"c >> ARGV[2];print n"\n"d >> ARGV[3];n=$1;c=$4;d=$6}}';
$get_distribution .= '(!x){n=$1;c=$4;d=$6;x=1}END{print n"\n"c >> ARGV[2];print n"\n"d >> ARGV[3]}\'';
my ($gcf,$depthf) = ("$outdir/gc_distribution","$outdir/depth_distribution");
`>$gcf;>$depthf;$get_distribution $outf.wind $gcf $depthf`;
#=============#
sub gc_depth
#=============#
#usage: gc_depth($infa,$windl,\%lenh,\%depth,$cover)
{
  my ($infa,$windl,$lenh,$depth,$cover) = @_;
  open IN,$infa;
  $/=">";<IN>;$/="\n";
  while(<IN>){
    /^(\S+)/ || next;
    my $id = $1;
    $/=">";
    chomp(my $seql = <IN>);
    $/="\n";
    my @seq;
    my $len;
    if($cover){
    	@seq = split/\s+/,$seql;
    	$len = @seq;
    	($len < $windl) && next;
    	$lenh->{$id} = $len;
    }else{
    	$len = ($lenh->{$id} || 0);
    	$len || next;
    	$seql =~ s/\s+//g;
    }
    my ($s, $e) = (1, $windl);
    until($len < $windl){
    	if($cover){
      		my @sub_seq = @seq[$s-1..$e-1];
     		depth_cacul(\@sub_seq,"$id $s",$depth);#sub1.1
    	}else{
    		my $sub_seq = substr($seql,$s-1,$windl);
      		gc_cacul($sub_seq,$id,$s,$e,($depth->{"$id $s"} || 0));#sub1.2
      	}
      $s += $step;
      $e += $step;
      $len -= $step;
    }
  }
  close IN;
}
#sub1.1
#==============#
sub depth_cacul
#==============#
{
	my ($seq,$id,$depth) = @_;
	my $dep = 0;
	foreach(@{$seq}){
		$_ || next;
		$dep += $_;
	}
	$dep || return(0);
	$depth->{$id} = $dep;
}

#sub1.2
#===========#
sub gc_cacul
#===========#
{
  $_ = shift;
  my ($id,$s,$e,$depth) = @_;
  my $gc = (s/[GC]//ig);
  my $at = (s/[AT]//ig);
  my $toal = $gc + $at;
  $toal || return(0);
  $gc /= $toal;
  $depth /= $toal;
  printf("%s\t%d\t%d\t%.4f\t%d\t%.4f\n",$id,$s,$e,$gc*100,$toal,$depth);
}

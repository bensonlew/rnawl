#!/usr/bin/perl -w
#Author: hao.gao\@majorbio.com
#Version: 1.0, date: 2011-5-30

#by Modify: xiaoyue.wang\@majorbio.com
#Version: 2.0, data: 2016-11-9
#The content of the modified: add -step, GC-Depth Image upgrade
#last modified by liulinmeng,20180511; change png file generation method CairoPNG

use strict;
use Getopt::Long;
my ($windl,$step,$outdir,$x_mun,$y_mun,$help);
GetOptions(
        "windl:i"=>\$windl,
	"step:i"=>\$step,
        "outdir:s"=>\$outdir,
        "x_mun:s"=>\$x_mun,
        "y_mun:s"=>\$y_mun,
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
  -x_mun  <str>        the x-axis scale min,unit,number, default 0,20,5
  -y_mun  <str>        the y-axis scale min,unit,number, default 0,150,4
  -help                output help information to screen

Example:
 nohup perl gc_depth.pl arowana.fa arowana.coverage.detail -outdir GC_depth -windl 5000 &
 
Notic:
 1 -x(y)_mun set can be ' ', then the scale will set accroding to the real data.
 2 draw.sh will creat at current dir, you can change some optins and run it to rewrite the figure.
 3 800M genomics use RES 800m and running time 15min
 4 results at outdir:
   gc_distribution, depth_distribution, gc_depth.wind and gc_depth.svg, gc_depth.png.\n";
## check error  ##
foreach(@ARGV){(-s $_) || die"Error: can't find file $_\n";}
my ($fasta,$coverf) = @ARGV;
## get options values  ##
$outdir ||= '.';
(-d $outdir) || mkdir"$outdir";
$windl ||= 10000;
$step ||= 500;
$x_mun ||= '0,20,5';
$y_mun ||= '0,150,4';
my @x_mun;my @y_mun;
@x_mun = split /,/,$x_mun;
@y_mun = split /,/,$y_mun;
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
## draw figure ##
open DR,">$outdir/depth_gc.cmd.r";
print DR "library(Cairo)
mydata <- read.table(\"$outf.wind\" ,header=T, sep=\"\\t\")
mydata <- mydata[,c(4,6)]\nnames(mydata) <- c(\"GC\",\"Depth\")
mycols <- densCols(mydata, colramp=colorRampPalette(c(\"black\", \"white\")), nbin = 5000)
mydata\$dens <- col2rgb(mycols)[1,] + 1L
cols <- colorRampPalette(c(\"gray\", \"orange\", \"red\"), space = \"Lab\")(256)
mydata\$col <- cols[mydata\$dens]
submydata <- mydata[(mydata\$GC >= $x_mun[0] & mydata\$GC <= $x_mun[0]+$x_mun[1]*$x_mun[2]) & (mydata\$Depth >= $y_mun[0] & mydata\$Depth <= $y_mun[0]+$y_mun[1]*$y_mun[2]),]
svg(\"$outf.svg\", width = 8, height = 8)
layout(matrix(c(2,0,1,3),2,2,byrow=TRUE), c(3,1), c(1,3), TRUE)
par(mar=c(6,6.3,1,1.2))
plot(Depth ~ GC, data=submydata[order(submydata\$dens),], col=col, ylab=\"Average depth (X)\", xlab=\"GC content (%)\", cex.lab = 1.4, xaxt=\"n\", yaxt=\"n\", pch = 20, cex=0.3, ylim = c($y_mun[0],$y_mun[1]*$y_mun[2]+$y_mun[0]), xlim = c($x_mun[0],$x_mun[0]+$x_mun[1]*$x_mun[2]))
axis(1, at=seq($x_mun[0],$x_mun[0]+$x_mun[1]*$x_mun[2],$x_mun[1]), cex.axis = 1.3)
axis(1, at=seq($x_mun[0],$x_mun[0]+$x_mun[1]*$x_mun[2],$x_mun[1]/4), labels=NA, cex.axis = 0.5)
axis(2, at=seq($y_mun[0],$y_mun[0]+$y_mun[1]*$y_mun[2],$y_mun[1]), cex.axis = 1.3)
axis(2, at=seq($y_mun[0],$y_mun[0]+$y_mun[1]*$y_mun[2],$y_mun[1]/4), labels=NA, cex.axis = 0.5)
par(mar=c(0,6.3,1,1.2))\nGChist <- hist(submydata\$GC, plot=FALSE ,breaks=seq($x_mun[0],$x_mun[0]+$x_mun[1]*$x_mun[2],$x_mun[1]*$x_mun[2]/125))
barplot(GChist\$counts, axes=T, space=0, col = \"plum\")
par(mar=c(6,0,1,1.2))\nDephist <- hist(submydata\$Depth, plot=FALSE ,breaks=seq($y_mun[0],$y_mun[0]+$y_mun[1]*$y_mun[2],$y_mun[1]*$y_mun[2]/125))
barplot(Dephist\$counts, axes=T, space=0, horiz=TRUE, col = \"plum\")
dev.off()
#png(\"$outf.png\", width = 21, height = 20, units = \"cm\", res = 600)
CairoPNG(file =\"$outf.png\", width = 21, height = 20, units = \"cm\", res = 600)
layout(matrix(c(2,0,1,3),2,2,byrow=TRUE), c(3,1), c(1,3), TRUE)
par(mar=c(6,6.3,1,1.2))
plot(Depth ~ GC, data=submydata[order(submydata\$dens),], col=col, ylab=\"Average depth (X)\", xlab=\"GC content (%)\", cex.lab = 1.4, xaxt=\"n\", yaxt=\"n\", pch = 20, cex=0.3, ylim = c($y_mun[0],$y_mun[1]*$y_mun[2]+$y_mun[0]), xlim = c($x_mun[0],$x_mun[0]+$x_mun[1]*$x_mun[2]))
axis(1, at=seq($x_mun[0],$x_mun[0]+$x_mun[1]*$x_mun[2],$x_mun[1]), cex.axis = 1.3)
axis(1, at=seq($x_mun[0],$x_mun[0]+$x_mun[1]*$x_mun[2],$x_mun[1]/4), labels=NA, cex.axis = 0.5)
axis(2, at=seq($y_mun[0],$y_mun[0]+$y_mun[1]*$y_mun[2],$y_mun[1]), cex.axis = 1.3)
axis(2, at=seq($y_mun[0],$y_mun[0]+$y_mun[1]*$y_mun[2],$y_mun[1]/4), labels=NA, cex.axis = 0.5)
par(mar=c(0,6.3,1,1.2))\nGChist <- hist(submydata\$GC, plot=FALSE ,breaks=seq($x_mun[0],$x_mun[0]+$x_mun[1]*$x_mun[2],$x_mun[1]*$x_mun[2]/125))
barplot(GChist\$counts, axes=T, space=0, col = \"plum\")
par(mar=c(6,0,1,1.2))\nDephist <- hist(submydata\$Depth, plot=FALSE ,breaks=seq($y_mun[0],$y_mun[0]+$y_mun[1]*$y_mun[2],$y_mun[1]*$y_mun[2]/125))
barplot(Dephist\$counts, axes=T, space=0, horiz=TRUE, col = \"plum\")
dev.off()\n";
#`R --restore --no-save < draw.r`;
#`rm draw.r`;
## get gc and depth distribution file ##
my $get_distribution = 'awk \'BEGIN{ARGC=2}(x){if(n==$1){c=c" "$4;d=d" "$6}else{print n"\n"c >> ARGV[2];print n"\n"d >> ARGV[3];n=$1;c=$4;d=$6}}';
$get_distribution .= '(!x){n=$1;c=$4;d=$6;x=1}END{print n"\n"c >> ARGV[2];print n"\n"d >> ARGV[3]}\'';
my ($gcf,$depthf) = ("$outdir/gc_distribution","$outdir/depth_distribution");
`>$gcf;>$depthf;$get_distribution $outf.wind $gcf $depthf`;
##====  process end  ====##

##=============================================##
##               SUB FUNCTION                  ##
##=============================================##
#sub1
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

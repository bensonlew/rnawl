#!/usr/bin/perl -w
use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use File::Basename qw(basename dirname);

use FindBin qw($Bin $Script);
my ($orginLnc,$chain,$out);
GetOptions(
	"fa:s"=>\$orginLnc,
	"chain:s"=>\$chain,
	"help|?" =>\&USAGE,
	"out:s"=>\$out)or &USAGE;
&USAGE unless ($orginLnc and $chain and $out);


mkdir $out if (!-d $out);
$out = ABSOLUTE_DIR($out);
$chain = ABSOLUTE_DIR($chain);
$orginLnc = ABSOLUTE_DIR($orginLnc);


open IN ,"<$orginLnc";
open OUT,">$out/lncRNA.bed";
my $i = 1;
my %lnc;
while (<IN>) {
	chomp;
	next if ($_ eq ""|| /^$/);
#	my @line = split /\t/;
	$lnc{$i} = $_;
	print OUT "$_\t$i\n";
	$i++;
}
close IN;
close OUT;
#print Dumper %lnc;
#die;

`liftOver $out/lncRNA.bed $chain $out/lift.out $out/unlift.out`;

#open SH ,">$out/lift.sh";
#print SH "liftOver $out/lncRNA.bed $chain $out/lift.out $out/unlift.out";
#close SH;
#my $queue||="DNA";
#my $job="qsub-slurm.pl --Queue $queue --Resource mem=1G --CPU 1  $out/lift.sh";
#`$job`;

JUGE: if (-e "$out/lift.out") {
	my $j=1;
	open IN, "<$out/lift.out";
	open OUT, ">$out/liftover.out";
	print OUT "id\tchr_target\tstart\tend\tchr_query\tstart\tend\n";
	while (<IN>) {
		chomp;
		next if ($_ eq ""|| /^$/);
		my @line = split /\t/;
		if (exists $lnc{$line[-1]}) {
			my $orgin = $lnc{$line[-1]};
			my @trans = split (/\t/,$orgin);
			pop @trans;
	#		print $lnc{$line[-1]};
#			die;
			my $num = pop @line;
			my $gene = pop @line;
			 
			print OUT $gene,"\t",join("\t",@line),"\t",join("\t",@trans),"\n";
		}else{
			next;
		}
	}
}else{
	sleep (2);
	last JUGE;
}
close IN;
close OUT;
`rm $out/lncRNA.bed`;

#######################################################################################
sub ABSOLUTE_DIR #$pavfile=&ABSOLUTE_DIR($pavfile);
{
	my $cur_dir=`pwd`;chomp($cur_dir);
	my ($in)=@_;
	my $return="";
	if(-f $in){
		my $dir=dirname($in);
		my $file=basename($in);
		chdir $dir;$dir=`pwd`;chomp $dir;
		$return="$dir/$file";
	}elsif(-d $in){
		chdir $in;$return=`pwd`;chomp $return;
	}else{
		warn "Warning just for file and dir \n$in";
		exit;
	}
	chdir $cur_dir;
	return $return;
}

sub USAGE {#
        my $usage=<<"USAGE";
Contact:        licui.li\@majorbio.com;
Description: Trans coordinate of different genomes with liftOver 
Usage:	perl CoordTrans.pl -fa <file> -chain <chain file>  -out <output file>
  Options:
  -fa	<file>	input, lncRNA with bed format
  -chain	<chain file>	input, chain file for liftOver
  -out	<output> output, new file with lift successful
  -h         Help

USAGE
        print $usage;
        exit;
}



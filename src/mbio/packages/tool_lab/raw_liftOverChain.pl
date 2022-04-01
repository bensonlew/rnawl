#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($lnRNA,$query,$target,$output,$queue);
my (@targetId,@lavlist,@falist,@psllist,@pslList,@chainlist);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use File::Path qw(make_path);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"target:s"=>\$query,
	"ref:s"=>\$target,
	"output:s"=>\$output,
	"queue:s"=>\$queue,
			) or &USAGE;
&USAGE unless ( $query and $target and $output);

mkdir $output if (!-d $output);
$output=ABSOLUTE_DIR($output);
$query=ABSOLUTE_DIR($query);
$target=ABSOLUTE_DIR($target);

`mkdir $output/target` ;
open IN, "<$target";
my $nSeq = 0;
my $chr;
$/ = ">";
while (<IN>) {
	chomp;
	next if ($_ eq ""||/^$/);
	$nSeq++;
	my ($id, @seq) = split(/\n/,$_);
	$id=(split(/\s+/,$id))[0];
	$chr = "chr".$nSeq.".fa";
	push @falist, $chr ;
	open FA, ">$output/target/$chr";
	print FA ">$id\n",join("\n",@seq),"\n";
	close FA;
	my $chainId =$id.".chain";
	push @targetId,$chainId;
	
}
close IN;
$/ = "\n";

`mkdir "$output/lastz"` ;
`mkdir "$output/Info"`;
open SH,">$output/01.Info.sh";
print SH "faSplit -lift=$output/Info/target.lft size $target -oneFile 5000000  -extra=10000 $output/Info/target","&&";
print SH "faSplit -lift=$output/Info/query.lft size $query -oneFile 5000000  -extra=10000 $output/Info/query","&&";
print SH "faToTwoBit $target  $output/Info/target.2bit","&&";
print SH "faToTwoBit $query  $output/Info/query.2bit","&&";
print SH "twoBitInfo $output/Info/target.2bit $output/Info/target.chrom.sizes","&&";
print SH "twoBitInfo $output/Info/query.2bit $output/Info/query.chrom.sizes","\n";
close SH;

open SH, ">$output/02.lastz.sh";
foreach my $fa (@falist) {
	print SH "c-1.04.00 $output/target/$fa $output/Info/query.fa --masking=50 --hspthresh=2200 --ydrop=3400 --gappedthresh=4000 --inner=2000 --ambiguous=iupac --output=$output/lastz/$fa.query.lav","\n";
	my $fa_query_lav = $fa.".query.lav";
	push @lavlist,$fa_query_lav;
}
close SH;

`mkdir $output/psl`;
`mkdir $output/psl/psl`;
`mkdir $output/chain`;
`mkdir $output/chain/chain`;
`mkdir $output/net`;
`mkdir $output/subset`;

open SH,">$output/03.chain.sh";

if (@lavlist){
	foreach my $lav (@lavlist) {
		print SH "lavToPsl $output/lastz/$lav $output/psl/$lav.psl","&&","\n";
		my $lav_psl = $lav.".psl";
		push @pslList,$lav_psl;
	}
}
if (@pslList) {
	foreach my $psl (@pslList) {
		print SH "liftUp -pslQ $output/psl/psl/query.$psl $output/Info/query.lft warn $output/psl/$psl","&&","\n";
		my $query_psl = "query.".$psl;
		push @psllist,$query_psl;
	}
}
if (@psllist) {
	foreach my $psllist (@psllist) {
		print SH "axtChain -linearGap=medium -psl $output/psl/psl/$psllist $output/Info/target.2bit  $output/Info/query.2bit $output/chain/$psllist.chain","&&","\n";
		my $psllist_chain = $psllist.".chain";
		push @chainlist,$psllist_chain;
	}
}
if (@chainlist) {
	print SH "chainMergeSort $output/chain/*.chain | chainSplit $output/chain/chain/ stdin && cat $output/chain/chain/*.chain > $output/chain/chain/target.chain","&&","\n";
}
close SH;

open SH ,">$output/04.net.sh";

if (@targetId) {
	foreach my $targetId(@targetId) {
		print SH "chainNet $output/chain/chain/$targetId $output/Info/target.chrom.sizes $output/Info/query.chrom.sizes $output/net/$targetId.net /dev/null","&&","\n";
	}
	print SH "cat $output/net/*.net > $output/net/target.net","&&","\n";

	print SH "netChainSubset $output/net/target.net $output/chain/chain/target.chain $output/subset/target.chain","&&","\n";

	print SH "less $output/subset/target.chain |grep -v \"\#\"> $output/subset/lift.chain","&&","\n";

	#	print SH "liftOver $lnRNA $output/subset/lift.chain $output/lift.bed $output/unlift.bed","\n";
}

close SH;

$queue||="DNA";
my $proc||=20;
my $job="qsub-slurm.pl --Queue $queue --Resource mem=20G --CPU 1 --maxjob $proc $output/01.Info.sh";
 print $job;
 `$job`;
$job="qsub-slurm.pl --Queue $queue --Resource mem=100G --CPU 1 --maxjob $proc $output/02.lastz.sh";
 print $job;
 `$job`;
$job="qsub-slurm.pl --Queue $queue --Resource mem=20G --CPU 1 --One --maxjob $proc $output/03.chain.sh";
 print $job;
 `$job`;
$job="qsub-slurm.pl --Queue $queue --Resource mem=20G --CPU 1 --One --maxjob $proc $output/04.net.sh";
 print $job;
 `$job`;

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
Script:			$Script
Description:
Usage: perl liftOverChain.pl  -ref <file> -target <file> -output <dir>
  Options:
  -ref	<file>	input Old genomic sequence
  -target	<file>	input New	 genomic sequence
  -output	<dir>	output dir 
  -h         Help

USAGE
        print $usage;
        exit;
}

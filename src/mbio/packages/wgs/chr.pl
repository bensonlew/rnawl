#!/usr/bin/env perl
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fIn,$fOut);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$fIn,
	"o:s"=>\$fOut,
			) or &USAGE;
&USAGE unless ($fIn and $fOut);
open In,$fIn;
my %chr;
while (<In>) {
	chomp;
	next if ($_ eq "" || /^$/ || !/^\@SQ/);
	my ($SQ,$SN,$LN,undef)=split(/\s+/,$_);
	my $id=(split(/\:/,$SN))[-1];
	my $len=(split(/\:/,$LN))[-1];
	$chr{$id}=$len;
}
close In;
open Out,">$fOut";
if (exists $chr{chr1}) {
	my %order;
	foreach my $id (keys %chr) {
		if ($id =~ /chr(\d+)/) {
			$order{$1}=$id;
		}
	}
	foreach my $id (sort {$a<=>$b} keys %order) {
		print Out $order{$id},"\t",$chr{$order{$id}},"\n";
	}
}else{
	my $n=0;
	foreach my $id (sort{$chr{$b}<=>$chr{$a}}keys %chr) {
		print Out $id,"\t",$chr{$id},"\n";
		$n++;
		last if ($n == 20);
	}
}
close Out;
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
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
Contact:        long.huang\@majorbio.com;
Script:			$Script
Description:

Usage:
  Options:
  -i	<file>	input dict name
  -o	<file>	out chr list name
  -h         Help

USAGE
        print $usage;
        exit;
}

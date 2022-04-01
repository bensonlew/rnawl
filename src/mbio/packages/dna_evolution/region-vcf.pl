#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fIn,$fOut,$select);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$fIn,
	"r:s"=>\$select,
	"o:s"=>\$fOut,
			) or &USAGE;
&USAGE unless ($fIn and $fOut and $select);
open In,$select;
my %region;
while (<In>) {
	chomp;
	next if ($_ eq ""|| /^$/ || /#/);
	my ($chr,$start,$end)=split(/\s+/,$_);
	$region{$chr}{join("\t",sort{$a<=>$b}($start,$end))}=1;
}
close In;
open In,$fIn;
my $head;
my %info;
my %stat;
my %einfo;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/);
	if (/^#/) {
		$head=$_;
	}else{
		my ($chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$ann,undef)=split(/\s+/,$_);
		foreach my $region (sort keys %{$region{$chr}}) {
			my ($pos1,$pos2)=split(/\s+/,$region);
			if ($pos >= $pos1 && $pos <= $pos2) {
				push @{$info{$chr}{$region}},$_;
				if($info=~/ANN=([^\;]*)/){$ann=$1;my @ann=split(/\|/,$ann);$ann=join("|",$ann[1],$ann[2],$ann[3],$ann[4])}
				my ($fun,$impact,$genename,$geneid)=split(/\|/,$ann);
				push @{$einfo{$chr}{$region}},$_ if($impact eq "HIGH" || $impact eq "MODERATE");
				$impact||="unknown";
				my @alt=split(/\,/,$alt);
				my $len=0;
				for (my $i=0;$i<@alt;$i++) {
					if ($len < length($alt[$i])){
						$len=length($alt[$i]);
					}
				}
				if ($len ne length($ref)) {
					$stat{$chr}{$region}{indel}++;
					$stat{$chr}{$region}{effindel}++ if($impact eq "HIGH" || $impact eq "MODERATE");
				}else{
					$stat{$chr}{$region}{snp}++;
					$stat{$chr}{$region}{effsnp}++ if($impact eq "HIGH" || $impact eq "MODERATE");
				}
			}
		}
	}
}
close In;
open Out,">$fOut.total";
print Out "#\@chr\tpos1\tpos2\tsnp\teffsnp\n";
print Out "$head\n";
foreach my $chr (sort keys %region) {
	foreach my $region (sort keys %{$region{$chr}}) {
		$stat{$chr}{$region}{snp}||=0;
		$stat{$chr}{$region}{effsnp}||=0;
		print Out join("\t","\@$chr",$region,$stat{$chr}{$region}{snp},$stat{$chr}{$region}{effsnp}),"\n";
		print Out join("\n",@{$info{$chr}{$region}}),"\n";
	}
}
close Out;
open Out,">$fOut.eff";
print Out "#\@chr\tpos1\tpos2\n";
print Out "$head\n";
foreach my $chr (sort keys %region) {
	foreach my $region (sort keys %{$region{$chr}}) {
		print Out join("\t","\@$chr",$region),"\n";
		print Out join("\n",@{$einfo{$chr}{$region}}),"\n";
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
	fq thanslate to fa format
	eg:
	perl $Script -i -o -k -c

Usage:
  Options:
  -i	<file>	input file name
  -o	<file>	output file name
  -r	<file>	input region file
  -h         Help

USAGE
        print $usage;
        exit;
}

#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fin,$fout,$select);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$fin,
	"o:s"=>\$fout,
	"r:s"=>\$select,
			) or &USAGE;
&USAGE unless ($fout);
open In,$select;
my %region;
while (<In>) {
	chomp;
	next if ($_ eq ""|| /^$/ || /#/);
	my ($chr,$start,$end)=split(/\s+/,$_);
	$region{$chr}{join("\t",sort{$a<=>$b}($start,$end))}=1;
}
close In;
open In,$fin;
my $head;
my %info;
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
			}
		}
	}
}
close In;
#print Dumper \%info;
my %filehandle;
open OUT,">$fout/vcf.list";
foreach my $chr (sort keys %info) {
	foreach my $region (sort keys %{$info{$chr}}) {
		my ($pos1,$pos2)=split/\s+/,$region;
		my $file = join("_",$chr,$pos1,$pos2);
		if (!exists $filehandle{$file}) {
			open $filehandle{$file},">$fout/$file.LD.vcf";
			print {$filehandle{$file}} "$head\n";
		}
		print {$filehandle{$file}} join("\n",@{$info{$chr}{$region}}),"\n";
		print OUT "$fout/$file.LD.vcf\n";
	}
}
close OUT;
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub ABSOLUTE_DIR 
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
Contact:        minghao.zhang\@majorbio.com;
Script:			$Script
Description:
	fq thanslate to fa format
	eg:
	perl $Script -i -o -k -c

Usage:
  Options:
	-i		input region.variant.total.vcf		       
    -o		output dir
	-r		input select.region
    -h      Help

USAGE
        print $usage;
        exit;
}

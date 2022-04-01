#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($wt,$mut,$output,$region);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"wt:s"=>\$wt,
	"mut:s"=>\$mut,
	"region:s"=>\$region,
	"output:s"=>\$output,
			) or &USAGE;
&USAGE unless ($mut and $region);
open In,$region;
my %region;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/);
	my ($chr,$start,$end,$regionid,$seq,undef)=split(/\s+/,$_);
	$region{$regionid}{position}="$chr:$start-$end";
	$region{$regionid}{seq}="$seq";
}
close In;
my %reference;
if (defined $wt) {
	my @wtfile=glob("$wt/*/*/Alleles_frequency_table.txt");
	#print Dumper @wtfile;
	foreach my $wtfile (@wtfile) {
		my $wtregion=(split(/\//,$wtfile))[-2];
		$wtregion=(split("on_",$wtregion))[-1];
		next if (!exists $region{$wtregion});
		open In,$wtfile;
		while (<In>) {
			chomp;
			next if ($_ eq ""||/^$/||/#/);
			my ($alle,$ref,$NHEJ,$UNMODIFIED,$HDR,undef,undef,undef,$reads,$reads2)=split;
			$reference{$wtregion}{$alle}=1;
		}
		close In;
	}
}
my @mutfile=glob("$mut/*/*/Alleles_frequency_table.txt");
#print Dumper @mutfile;
my %mut;
my %mutstat;
foreach my $mutfile (@mutfile) {
		my $mutregion=(split(/\//,$mutfile))[-2];
		$mutregion=(split("on_",$mutregion))[-1];
		open In,$mutfile;
		while (<In>) {
			chomp;
			next if ($_ eq ""||/^$/||/#/);
			my ($alle,$ref,$NHEJ,$UNMODIFIED,$HDR,undef,undef,undef,$reads,$reads2)=split;
			if (!defined $wt) {
				$reference{$mutregion}{$ref}=1;
			}
			if (exists $reference{$mutregion}{$alle}){
				$mut{$mutregion}{$alle}{reads}+=$reads;
				$mut{$mutregion}{$alle}{type}="unmodified";
			}else{
				$mut{$mutregion}{$alle}{reads}+=$reads;
				$mut{$mutregion}{$alle}{type}="modified";
			}
			$mutstat{$mutregion}{$mut{$mutregion}{$alle}{type}}++;
		}
		close In;
}
open Out,">$output.detail";
print Out join("\t","#postion","amplication","alle reads","alle dep","alle type","wt alle"),"\n";
foreach my $region (keys %mut) {
	foreach my $ale (sort keys %{$mut{$region}}) {
		print Out join("\t",$region{$region}{position},$region{$region}{seq},$ale,$mut{$region}{$ale}{reads},$mut{$region}{$ale}{type},join(":",sort keys %{$reference{$region}})),"\n";
	}
}
close Out;
open Out,">$output.stat";
print Out join("\t","#postion","amplication","modified","unmodifed","modified ratio"),"\n";
foreach my $region (keys %mut) {
	$mutstat{$region}{modified}||=0;
	$mutstat{$region}{unmodified}||=0;
	print Out join("\t",$region{$region}{position},$region{$region}{seq},$mutstat{$region}{modified},$mutstat{$region}{unmodified},$mutstat{$region}{modified}/($mutstat{$region}{modified}+$mutstat{$region}{unmodified})),"\n";
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
  -wt	<file>	input wt sample file dir
  -mut	<file>	input mut sample file dir
  -region	<file>	input region file
  -output	<file>	input gff file name
  -h         Help

USAGE
        print $usage;
        exit;
}

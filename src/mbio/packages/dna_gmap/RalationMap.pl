#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Statistics::RankCorrelation;
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fMap,$fAlign,$fKey,$dOut);
GetOptions(
				"help|?" =>\&USAGE,
				"m:s"=>\$fMap,
				"k:s"=>\$fKey,
				"o:s"=>\$dOut
				) or &USAGE;
&USAGE unless ($fMap and $fKey );

#die ` pod2text $0 ` unless ($fMap and $fAlign and $fChroLength and $fKey) ;

$dOut||="./";
mkdir $dOut unless (-d $dOut) ;
my $ConvertSVGtoPNG="$Bin/svg2xxx_release/svg2xxx";


#-------------------------------------------------------------------
# Global value
#-------------------------------------------------------------------

my (%Map,%Marker,%Alignment,%chroLength) = ();

#-------------------------------------------------------------------
# Get Map info
#-------------------------------------------------------------------

open (IN,"<",$fMap) or die $!;
$/="\n";
my ($group,$order);
my %spearman;
while (<IN>) {
	chomp;
	s/\r//g;
	next if (/^$/ || /^;/ || /^\#/) ;
	if (/^group\s+(\d+)/) {
		$group=$1;
		$order=0;
	}else{
		my ($marker,$cm,undef)=split(/\s+/,$_);
		$Map{$group}{$marker}{"order"}=$order;
		$Map{$group}{$marker}{"cm"}=$cm;
		$Marker{$marker}{"G"}{"group"}=$group;
		$Marker{$marker}{"G"}{"order"}=$order;
		$Marker{$marker}{"G"}{"cm"}=$cm;
		my ($chr,$pos)=split(/\_/,$marker);
		$chr=~s/\D//g;
		$Marker{$marker}{"C"}{$chr}=$pos;
		$Alignment{$chr}{$marker}=$pos;
		$chroLength{$chr}=$pos if (!exists $chroLength{$chr} || $chroLength{$chr} < $pos);
		$order++;
		$spearman{$group}{$chr}++;
	}
}
close (IN) ;
open Out,">$dOut/$fKey.spearman.xls";
print Out "#LGID\tMarkerNum\tMaxChr\tSpearman\n";
foreach my $chr (sort keys %spearman) {
	my $group=(sort{$spearman{$chr}{$b}<=>$spearman{$chr}{$a}} keys %{$spearman{$chr}})[0];
	my @marker1=sort{$Map{$chr}{$a}<=>$Map{$chr}{$b}}keys %{$Map{$chr}};
	my @marker2=sort{$Alignment{$group}{$a}<=>$Alignment{$group}{$b}} keys %{$Alignment{$group}};
	my (@m1,@m2);
	for (my $i=0;$i<@marker1;$i++) {
		push @m1,$i;
		push @m2,grep{$marker2[$_] eq $marker1[$i]} 0..$#marker2;
	}		
	my $c = Statistics::RankCorrelation->new( \@m1, \@m2,   );
	my $n=$c->spearman;
	print Out $chr,"\t",scalar @marker1,"\t",$chr,"\t",$n,"\n";
}
close Out;

sub USAGE {#
	my $usage=<<"USAGE";
Program: $0
Version: $version
Modified by HONGONG 
Contact: Li Tiancheng <litc\@biomarker.com.cn> <ltc_gs\@qq.com>

	map file: 
	  group X
	  MarkerID	MapDistance

Usage:
  Options:
  -m <file>  map file, JoinMap map format, forced.

USAGE
	print $usage;
	exit;
}

#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($gff , $anno , $fOut ,$select);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$select,
	"a:s"=>\$anno,
	"o:s"=>\$fOut,
			) or &USAGE;
&USAGE unless ($select and $anno and $fOut);
open In,$select;
my %region;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/ ||/#/);
	s/\"//g;
	my ($chr,$pos1,$pos2,undef)=split(/\s+/,$_);
	next if ($chr eq "chr");
	my $regioned=0;
	foreach my $region (sort keys %{$region{$chr}}) {
		my ($pos3,$pos4)=split(/\s+/,$region);
		#print $pos1,"\t",$pos2,"\t",$pos3,"\t",$pos4,"\n";
		if ($pos1 >= $pos3 && $pos1 <= $pos4) {
			my ($p1,$p2,$p3,$p4)=sort {$a<=>$b} ($pos1,$pos2,$pos3,$pos4);
			my $newregion=join("\t",$p1,$p4);
			delete $region{$chr}{$region};
			$region{$chr}{$newregion}++;
			$regioned=1;
		}
#		else{$regioned=0;}
	}
	if ($regioned == 0) {
		$region{$chr}{join("\t",$pos1,$pos2)}++;
	}
}
close In;
open In,$anno;
my %kdetail;
my %gdetail;
my %enrich;
my %edetail;
my %info;
my %stat;
my %astat;
my $head;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/);
	if (/^#/) {
		$head=$_;
		next;
	}
	my ($Gene_name,$Gene_id,$Transcript_id,$Bio_Type,$chr,$Pos1,$Pos2,$High,$Moderate,$Low,$Modifier,$nrid,$nranno,$uniid,$unianno,$koid,$koanno,$goid,$goanno,$egid,$expre)=split(/\t/,$_);
	my $regioned=0;
	foreach my $region (sort keys %{$region{$chr}}) {
		my ($pos3,$pos4)=split(/\t/,$region);
		if (($Pos1 > $pos3 && $Pos1 <$pos4) ||($Pos2 > $pos3 && $Pos2 < $pos4) || ($pos3 > $Pos1 && $pos3 < $Pos2) || ($pos4 > $Pos1 && $pos4 < $Pos2)) {
			push @{$info{$chr}{$region}},$_;
			$stat{$chr}{$region}{total}++;
			$stat{$chr}{$region}{totalnr}++ if($nrid ne "--");
			$stat{$chr}{$region}{totaluni}++ if($uniid ne "--");
			$stat{$chr}{$region}{totalkegg}++ if($koid ne "--");
			$stat{$chr}{$region}{totalgo}++ if($goid ne "--");
			$stat{$chr}{$region}{totaleggnog}++ if($egid ne "--");
			$stat{$chr}{$region}{eff}++ if($High+$Moderate > 0);
			$stat{$chr}{$region}{effnr}++ if($nrid ne "--" && $High+$Moderate > 0);
			$stat{$chr}{$region}{effuni}++ if($uniid ne "--" && $High+$Moderate > 0 );
			$stat{$chr}{$region}{effkegg}++ if($koid ne "--" && $High+$Moderate > 0);
			$stat{$chr}{$region}{effgo}++ if($goid ne "--" && $High+$Moderate > 0);
			$stat{$chr}{$region}{effeggnog}++ if($egid ne "--" && $High+$Moderate > 0);
			$regioned=1;
		}
	}
	$astat{total}++;
	$astat{eff}++ if($regioned == 1);

	my @koid=split(/,/,$koid);
	my @kdetail=split(/:/,$koanno);
	for (my $i=0;$i<@koid;$i++) {
		next if ($koid[$i] eq "--");
		$enrich{$koid[$i]}{total}++;
		$enrich{$koid[$i]}{enrich}++ if($regioned ==1);
		$kdetail{$koid[$i]}=$kdetail[$i];
	}
	my @goid=split(/,/,$goid);
	my @gdetail=split(/:/,$goanno);
	for (my $i=0;$i<@goid;$i++) {
		next if ($goid[$i] eq "--");
		$enrich{$goid[$i]}{total}++;
		$enrich{$goid[$i]}{enrich}++ if($regioned ==1);
		$gdetail{$goid[$i]}=$gdetail[$i];
	}
	my ($eid,$eanno)=split(/:/,$egid);
	$eanno||="";
	my @eid=split(/,/,$eid);
	my @edetail=split(/;/,$eanno);
	for (my $i=0;$i<@eid;$i++) {
		next if ($eid[$i] eq "--");
		$enrich{$eid[$i]}{total}++;
		$enrich{$eid[$i]}{enrich}++ if($regioned ==1);
		$edetail{$eid[$i]}=$edetail[$i];
	}
}
close In;
open Out,">$fOut.total";
print Out "#\@chr\tpos1\tpos2\ttotal\teff\ttotalnr\ttotaluni\ttotalkegg\ttotalgo\ttotaleggnog\teffnr\teffuni\teffkegg\teffgo\teffeggnog\n";
print Out "$head\n";
foreach my $chr (sort keys %region) {
	foreach my $region (sort keys %{$region{$chr}}) {
		$stat{$chr}{$region}{totalnr}||=0;
		$stat{$chr}{$region}{totaluni}||=0;
		$stat{$chr}{$region}{totalkegg}||=0;
		$stat{$chr}{$region}{totalgo}||=0;
		$stat{$chr}{$region}{totaleggnog}||=0;
		$stat{$chr}{$region}{effnr}||=0;
		$stat{$chr}{$region}{effuni}||=0;
		$stat{$chr}{$region}{effkegg}||=0;
		$stat{$chr}{$region}{effgo}||=0;
		$stat{$chr}{$region}{effeggnog}||=0;
		$stat{$chr}{$region}{total}||=0;
		$stat{$chr}{$region}{eff}||=0;
		print Out join("\t","\@$chr",$region,$stat{$chr}{$region}{total},$stat{$chr}{$region}{eff},$stat{$chr}{$region}{totalnr},$stat{$chr}{$region}{totaluni},$stat{$chr}{$region}{totalkegg},$stat{$chr}{$region}{totalgo},$stat{$chr}{$region}{totaleggnog},$stat{$chr}{$region}{effnr},$stat{$chr}{$region}{effuni},$stat{$chr}{$region}{effkegg},$stat{$chr}{$region}{effgo},$stat{$chr}{$region}{effeggnog}),"\n";
		next if (!defined $info{$chr}{$region});
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
		next if (!defined $info{$chr}{$region});
		foreach my $line (@{$info{$chr}{$region}}) {
			my @split=split(/\t/,$line);
			print Out $line,"\n" if($split[7]+$split[8] > 0);
		}
	}
}
close Out;

open Out,">$fOut.kegg.stat";
foreach my $koid (sort keys %kdetail) {
	$enrich{$koid}{enrich}||=0;
	$enrich{$koid}{total}||=0;
	$astat{totalkegg}||=0;
	$astat{effkegg}||=0;
	print Out join("\t",$koid,$kdetail{$koid},$enrich{$koid}{enrich},$enrich{$koid}{total},$astat{eff},$astat{total}),"\n";
}
close Out;
open Out,">$fOut.go.stat";
foreach my $goid (sort keys %gdetail) {
	$enrich{$goid}{enrich}||=0;
	$enrich{$goid}{total}||=0;
	$astat{totalgo}||=0;
	$astat{effgo}||=0;
	print Out join("\t",$goid,$gdetail{$goid},$enrich{$goid}{enrich},$enrich{$goid}{total},$astat{eff},$astat{total}),"\n";
}
close Out;
open Out,">$fOut.eggnog.stat";
foreach my $eggnog (sort keys %edetail) {
	$enrich{$eggnog}{enrich}||=0;
	$enrich{$eggnog}{total}||=0;
	$astat{totaleggnog}||=0;
	$astat{effeggnog}||=0;
	print Out join("\t",$eggnog,$edetail{$eggnog},$enrich{$eggnog}{enrich},$enrich{$eggnog}{total},$astat{eff},$astat{total}),"\n";
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
  -select	<file>	input select file name
  -anno	<file>	input anno file name
  -out	<key>	output keys of file name
  -h         Help

USAGE
        print $usage;
        exit;
}

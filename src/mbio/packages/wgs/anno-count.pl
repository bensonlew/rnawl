#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($snp , $indel , $anno , $fOut );
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"snp:s"=>\$snp,
	"indel:s"=>\$indel,
	"anno:s"=>\$anno,
	"out:s"=>\$fOut,
			) or &USAGE;
&USAGE unless ($snp and $indel and $anno and $fOut );
open In,$snp;
my %info;
my %stat;
my %eff;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/ || /^#/);
	my ($gname,$gid,$tid,$biotype,$high,$low,$moderate,$modifier,undef)=split(/\t/,$_);
	my $info=join("\t",$gname,$gid,$tid,$biotype);
	$info{$gname}=$info;
	$info{$gid}=$info;
	$info{$tid}=$info;
	$eff{$info}++ if($high+$moderate > 0);
	$stat{$info}{high}+=$high;
	$stat{$info}{low}+=$low;
	$stat{$info}{middle}+=$moderate;
	$stat{$info}{unknow}+=$modifier;
}
close In;
open In,$indel;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/ || /^#/);
	my ($gname,$gid,$tid,$biotype,$high,$low,$moderate,$modifier,undef)=split(/\t/,$_);
	my $info=join("\t",$gname,$gid,$tid,$biotype);
	$info{$gname}=$info;
	$info{$gid}=$info;
	$info{$tid}=$info;
	$eff{$info}++ if($high+$moderate > 0);
	$stat{$info}{high}+=$high;
	$stat{$info}{low}+=$low;
	$stat{$info}{middle}+=$moderate;
	$stat{$info}{unknow}+=$modifier;
}
close In;
open In,$anno;
open Out,">$fOut.summary";
my %kdetail;
my %gdetail;
my %enrich;
my %edetail;
my %pdetail;
my %fun;
print Out join("\t","#Gene_name","Gene_id","Transcript_id","Bio_Type","Chr","Pos1","Pos2","High","Moderate","Low","Modifier","NR-ID","NR-ANNO","Uni-ID","Uni-ANNO","KEGG-ID","KEGG-ANNO","GO-ID","GO-ANNO","EggNOG-ID","EggNOG-ANNO","PfamAccession","PfamAnno","InterProAccession","Anno"),"\n";
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/);
	if (/^#/) {
		next;
	}else{
		my ($id,$nrid,$nranno,$uniid,$unianno,$koid,$koanno,$goid,$goanno,$egid,$eganno,$pfid,$pfanno,$inid,$inanno)=split(/\t/,$_);
		my @ids=split(/:/,$id);
		if (!defined $ids[2] || !defined $ids[3] || !defined $ids[4]) {
			next;
		}
		my $pos=join("\t",$ids[2],$ids[3],$ids[4]);
		my $gname=(split(/\|/,$ids[1]))[0];
		my $tranid=(split(/\|/,$ids[1]))[-2];
		if($tranid eq "--"){
			$tranid = $ids[0];
		}
		$info{$ids[0]}||=join("\t",$gname,$ids[0],$ids[1],"--");
		my @info=split("\t",$info{$ids[0]});
		if (!defined $info[3]) {
			$info[3] = "--";
		}
		my $infos=join("\t",$gname,$ids[0],$ids[1],$info[3]);
		my $info=$info{$ids[0]};
		$stat{$info}{high}||=0;
		$stat{$info}{low}||=0;
		$stat{$info}{middle}||=0;
		$stat{$info}{unknow}||=0;
		$fun{total}{nr}++ if($nrid ne "--");
		$fun{eff}{nr}++ if($stat{$info}{high}+$stat{$info}{middle} > 0 && $nrid ne "--");
		$fun{total}{uni}++ if($uniid ne "--");
		$fun{eff}{uni}++ if($stat{$info}{high}+$stat{$info}{middle} > 0 && $uniid ne "--");
		$fun{total}{kegg}++ if($koid ne "--");
		$fun{eff}{kegg}++ if($stat{$info}{high}+$stat{$info}{middle} > 0 && $koid ne "--");
		$fun{total}{go}++ if($goid ne "--");
		$fun{eff}{go}++ if($stat{$info}{high}+$stat{$info}{middle} > 0 && $goid ne "--");
		$fun{total}{eggnog}++ if($egid ne "--");
		$fun{eff}{eggnog}++ if($stat{$info}{high}+$stat{$info}{middle} > 0 && $egid ne "--");
		$fun{total}{pfam}++ if($pfid ne "--");
		$fun{eff}{pfam}++ if($stat{$info}{high}+$stat{$info}{middle} > 0 && $pfid ne "--");
		print Out join("\t",$infos,$pos,$stat{$info}{high},$stat{$info}{middle},$stat{$info}{low},$stat{$info}{unknow},$nrid,$nranno,$uniid,$unianno,$koid,$koanno,$goid,$goanno,$egid,$eganno,$pfid,$pfanno,$inid,$inanno),"\n";
		my @koid=split(/,/,$koid);
		my @kdetail=split(/:/,$koanno);
		for (my $i=0;$i<@koid;$i++) {
			next if ($koid[$i] eq "--");
			$enrich{$koid[$i]}{total}++;
			if (exists $enrich{$koid[$i]}{total_gene}){
				$enrich{$koid[$i]}{total_gene} = "$enrich{$koid[$i]}{total_gene};$tranid";
			}else{
				$enrich{$koid[$i]}{total_gene} = $tranid;
			}
			$enrich{$koid[$i]}{enrich}++ if($stat{$info}{high}+$stat{$info}{middle} > 0);
			$kdetail{$koid[$i]}=$kdetail[$i];
		}
		my @goid=split(/,/,$goid);
		my @gdetail=split(/:/,$goanno);
		for (my $i=0;$i<@goid;$i++) {
			next if ($goid[$i] eq "--");
			$enrich{$goid[$i]}{total}++;
			if (exists $enrich{$goid[$i]}{total_gene}){
				$enrich{$goid[$i]}{total_gene} = "$enrich{$goid[$i]}{total_gene};$tranid";
			}else{
				$enrich{$goid[$i]}{total_gene} = $tranid;
			}
			$enrich{$goid[$i]}{enrich}++ if($stat{$info}{high}+$stat{$info}{middle} > 0);
			$gdetail{$goid[$i]}=$gdetail[$i];
		}
		my ($eid,$eanno)=split(/:/,$egid);
		if ($egid eq "--") {
			next;
		}
		my @eid=split(/,/,$eid);
		my @edetail=split(/;/,$eanno);
		for (my $i=0;$i<@eid;$i++) {
			next if ($eid[$i] eq "--");
			$enrich{$eid[$i]}{total}++;
			if (exists $enrich{$eid[$i]}{total_gene}){
				$enrich{$eid[$i]}{total_gene} = "$enrich{$eid[$i]}{total_gene};$tranid";
			}else{
				$enrich{$eid[$i]}{total_gene} = $tranid;
			}
			$enrich{$eid[$i]}{enrich}++ if($stat{$info}{high}+$stat{$info}{middle} > 0);
			$edetail{$eid[$i]}=$edetail[$i];
		}
		my @pid=split(/,/,$pfid);
		my @pdetail=split(/:/,$pfanno);
		for (my $i=0;$i<@pid;$i++) {
			next if ($pid[$i] eq "--");
			$enrich{$pid[$i]}{total}++;
			if (exists $enrich{$pid[$i]}{total_gene}){
				$enrich{$pid[$i]}{total_gene} = "$enrich{$pid[$i]}{total_gene};$tranid";
			}else{
				$enrich{$pid[$i]}{total_gene} = $tranid;
			}
			$enrich{$pid[$i]}{enrich}++ if($stat{$info}{high}+$stat{$info}{middle} > 0);
			$pdetail{$pid[$i]}=$pdetail[$i];
		}
	}
}
close Out;
close In;
open Out,">$fOut.kegg.stat";
print Out join("\t","#koid","ko_detail","eff_variant","all_gene","total_eff","total_gene","all_gene_list"),"\n";
foreach my $koid (sort keys %kdetail) {
	$enrich{$koid}{enrich}||=0;
	$enrich{$koid}{total}||=0;
	print Out join("\t",$koid,$kdetail{$koid},$enrich{$koid}{enrich},$enrich{$koid}{total},scalar keys %eff,scalar keys %stat,$enrich{$koid}{total_gene}),"\n";
}
close Out;
open Out,">$fOut.go.stat";
print Out join("\t","#goid","go_detail","eff_variant","all_gene","total_eff","total_gene","all_gene_list"),"\n";
foreach my $goid (sort keys %gdetail) {
	$enrich{$goid}{enrich}||=0;
	$enrich{$goid}{total}||=0;
	print Out join("\t",$goid,$gdetail{$goid},$enrich{$goid}{enrich},$enrich{$goid}{total},scalar keys %eff,scalar keys %stat,$enrich{$goid}{total_gene}),"\n";
}
close Out;
open Out,">$fOut.eggnog.stat";
print Out join("\t","#eggnogid","eggnog_detail","eff_variant","all_gene","total_eff","total_gene","all_gene_list"),"\n";
foreach my $eggnog (sort keys %edetail) {
	$enrich{$eggnog}{enrich}||=0;
	$enrich{$eggnog}{total}||=0;
	print Out join("\t",$eggnog,$edetail{$eggnog},$enrich{$eggnog}{enrich},$enrich{$eggnog}{total},scalar keys %eff,scalar keys %stat,$enrich{$eggnog}{total_gene}),"\n";
}
close Out;
open Out,">$fOut.pfam.stat";
print Out join("\t","#Pfam Accession","Pfam Annotation","eff_variant","all_gene","total_eff","total_gene","all_gene_list"),"\n";
foreach my $pfid (sort keys %pdetail) {
	$enrich{$pfid}{enrich}||=0;
	$enrich{$pfid}{total}||=0;
	$enrich{$pfid}{total_gene}||=0;
	print Out join("\t",$pfid,$pdetail{$pfid},$enrich{$pfid}{enrich},$enrich{$pfid}{total},scalar keys %eff,scalar keys %stat,$enrich{$pfid}{total_gene}),"\n";
}
close Out;
open Out,">$fOut.stat.csv";
$fun{total}{nr}||=0;
$fun{eff}{nr}||=0;
$fun{total}{uni}||=0;
$fun{eff}{uni}||=0;
$fun{total}{kegg}||=0;
$fun{eff}{kegg}||=0;
$fun{total}{go}||=0;
$fun{eff}{go}||=0;
$fun{total}{eggnog}||=0;
$fun{eff}{eggnog}||=0;
$fun{total}{pfam}||=0;
$fun{eff}{pfam}||=0;
print Out join("\t","#type","NR","Uniprot","KEGG","GO","EGGNOG","Pfam"),"\n";
print Out join("\t","total",$fun{total}{nr},$fun{total}{uni},$fun{total}{kegg},$fun{total}{go},$fun{total}{eggnog},$fun{total}{pfam}),"\n";
print Out join("\t","eff",$fun{eff}{nr},$fun{eff}{uni},$fun{eff}{kegg},$fun{eff}{go},$fun{eff}{eggnog},$fun{eff}{pfam}),"\n";
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
  -snp	<file>	input snp file name
  -indel	<file>	input indel file name
  -anno	<file>	input anno file name
  -out	<key>	output keys of file name
  -h         Help

USAGE
        print $usage;
        exit;
}

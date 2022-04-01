#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($dmap,$dOut,$adjust,$pop);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"dmap:s"=>\$dmap,
	"out:s"=>\$dOut,
	"pop:s"=>\$pop,
	"adjust"=>\$adjust,
			) or &USAGE;
&USAGE unless ($dmap and $dOut );
mkdir $dOut if (!-d $dOut);
my @map=glob("$dmap/*.out");
open Out,">$dOut/total.map";
my %Marker;
my %lg;
my %stat;
foreach my $map (@map) {
	my $lgID=(split(/\./,basename($map)))[0];
	$lg{$lgID}=1;
	if ($lgID =~ /(\d+)/) {
		$lgID=$1;
	}else{
		$lgID=scalar keys %lg;
	}
	#print Out "group\t",$lgID,"\n";
	
	open In,$map;
	my $max=0;
	my @order;
	while (<In>) {
		chomp;
		next if ($_ eq ""||/^$/ || /^;/ || /group/);
		my ($id,$pos)=split(/\t/,$_);
		$Marker{$id}=join(",",$lgID,$pos);
		push @order,$id;
		if ($max < $pos) {
			$max=$pos;
		}
	}
	close In;
	my $newdis=$max;
	if ($adjust) {
		$newdis=rand(60)+120;
	}
	$max+=0.0000001 if($max == 0);
	foreach my $id (@order) {
		my ($lgid,$pos)=split(/\,/,$Marker{$id});
		$pos=$newdis/$max*$pos;
		$Marker{$id}=join(",",$lgid,$pos);
		#print Out $id,"\t",$pos,"\n";
		my $in=join("\t",$id,$pos);
		push @{$stat{$lgID}{id}},$in;
	}	
}
foreach my $lgID (sort {$a<=>$b}keys %stat){
	print Out "group\t$lgID\n";
	print Out join("\n",@{$stat{$lgID}{id}}),"\n";
}
	
close Out;
my @marker=glob("$dmap/*.correct.marker");
open Out,">$dOut/total.marker";
open CSV,">$dOut/total.csv";
my $head;
my $chead;
my $nind;
my $nloc;
my $name;
my @out;
my @cout;
my @mout;
foreach my $marker (@marker) {
	open In,$marker;
	while (<In>) {
		chomp;
		next if ($_ eq ""|| /^$/);
		my @info=split;
		next if (scalar @info < 3);
		if (/MarkerID/) {
			$head=$_;
			my (undef,$nhead)=split(/\s+/,$_,2);
			$nhead=~s/\t/,/g;
			my @head=split(/\,/,$nhead);
			$nind=scalar @head;
			$chead="Genotype,,,".join(",",@head);
		}else{

			my ($id,$info)=split(/\s+/,$_,2);
			if (!exists $Marker{$id}) {
				next;
			}
			my $number=(split(/\,/,$Marker{$id}))[0];
			push @{$stat{$number}{out}},$_;
			$info=~s/\t/,/g;
			$info=~s/X/H/g;
			$info=~s/U/-/g;
			#push @cout,join(",",$id,$Marker{$id},$info);
			$nloc++;
			push @{$stat{$number}{info}},join(",",$id,$Marker{$id},$info);
			$info=~s/H/h/g;
			$info=~s/U/-/g;
			$info=~s/A/a/g;
			$info=~s/B/b/g;
			$info=~s/,/\t/g;
			push @mout,join("\t",$id,$info),"\n";
		}
	}
	close In;
}
print Out $head,"\n";
print CSV $chead,"\n";
foreach my $number(sort {$a<=>$b} keys %stat ){
	print Out join("\n",@{$stat{$number}{out}}),"\n";
	print CSV join("\n",@{$stat{$number}{info}}),"\n";
}
#print CSV join("\n",$chead,@cout);
close Out;
close CSV;
open LOC,">$dOut/total.loc";
print LOC "nind=$nind\nnloc=$nloc\nname=pop\npopt=$pop\n";
print LOC join("n",@mout);
close LOC;
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
  -dmap	<file>	input file name
  -out	<file>	output file
  -adjust	<file>	adjust map distcance
  -h         Help

USAGE
        print $usage;
        exit;
}

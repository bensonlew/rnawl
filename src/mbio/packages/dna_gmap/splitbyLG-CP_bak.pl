#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Cwd qw(abs_path);
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fIn,$dOut,$fLG,$fKey,$type);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$fIn,
				"l:s"=>\$fLG,
				"d:s"=>\$dOut,
				"t:s"=>\$type,
				) or &USAGE;
&USAGE unless ($fIn and $dOut and $fLG);
open In,$fIn;
my %info;
my $head;
my $nind;
my %type;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/);
	if (/^#/) {
		$head=$_;
		next;
	}
	my ($id,$type,$info)=split(/\s+/,$_,3);
	$info{$id}=$info;
	$type{$id}=$type;
	my @ind=split(/\t/,$info);
	$nind=scalar @ind;
}
close In;
open In,$fLG;
open List,">$dOut/pri.marker.list";
$/=">";
while (<In>) {
	chomp;
	next if ($_ eq "" ||/^$/);
	my ($id,$marker)=split(/\n/,$_);
	$id=(split(/\s+/,$id))[0];
	open Out,">$dOut/$id.pri.marker";
	open Map,">$dOut/$id.pri.map";
	print Map "group $id\n";
	print List "$id\t$dOut/$id.pri.marker\n";
	my @marker=split(/\s+/,$marker);
	my %pos;
	for (my $i=0;$i<@marker;$i++) {
		$pos{$marker[$i]}=$i;
		if ($marker[$i] =~ /\_/) {
			$pos{$marker[$i]}=(split(/\_/,$marker[$i]))[-1];
			if ($pos{$marker[$i]}=~/\-/) {
				$pos{$marker[$i]}=(split(/\-/,$pos{$marker[$i]}))[0];
			}

		}
	}
	my @out;
	my $nloc=scalar @marker;
	foreach my $m (sort{$pos{$a}<=>$pos{$b}} keys %pos) {
		if (!exists $info{$m}) {
			next;
		}
		push @out,join(" ",$m,"<".$type{$m}.">","{--}",$info{$m});
		print Map join("\t",$m,$pos{$m}),"\n";
	}
#	print Out $head,"\n";
	print Out join("\n",@out),"\n";
	close Out;
	close Map;
}
close In;
close List;
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
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

sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

sub USAGE {#
	my $usage=<<"USAGE";
Program: $0
Version: $version
Contact: huangl <long.huang\@majorbio.com> 

Usage:  ��genotype�ļ�������Ⱥ�ָ�
  Options:
  -help			USAGE,
  -i	genotype file�� forced
  -l	linkage lg file
  -o	output dir
  
   
USAGE
	print $usage;
	exit;
}


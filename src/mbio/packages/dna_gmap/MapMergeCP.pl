#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($dmap,$dOut,$marker,$adjust);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"dmap:s"=>\$dmap,
	"marker:s"=>\$marker,
	"out:s"=>\$dOut,
	"adjust"=>\$adjust,
			) or &USAGE;
&USAGE unless ($dmap and $dOut and $marker);
mkdir $dOut if (!-d $dOut);
open In,$marker;
my @Indi;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/);
	if (/MarkerID/) {
		my (undef,undef,@indi)=split;
		push @Indi,@indi;
		last;
	}
}
close In;
my @map=glob("$dmap/*.sexAver.map");
open Out,">$dOut/total.sexAver.map";
#open Draw,">$dOut/total.sexAver.map.draw";
my %info;
foreach my $map (@map) {
	my $lgID=(split(/\./,basename($map)))[0];
	$lgID=~s/\D+//g;
	my $nloc=`wc -l $map`;
	$nloc=(split(/\s+/,$nloc))[0];
	next if ($nloc <= 2);
	print Out "group\t$lgID\n";
	open In,$map;
	my $max=0;
	while (<In>) {
		chomp;
		next if ($_ eq ""||/^$/ || /^;/ || /group/);
		my ($id,$pos)=split(/\s+/,$_);
		$info{$id}{lgID}=$lgID;
		$info{$id}{pos}=$pos;
		$max=$pos if ($pos > $max);
	}
	close In;
	my $newdis=$max;
	if ($adjust) {
		$newdis=rand(80)+120;
	}
	my $n=0;
	my $pos0;
	foreach my $id (sort {$info{$a}{pos}<=>$info{$b}{pos}} keys %info) {
		next if ($info{$id}{lgID} ne $lgID);
		$info{$id}{pos}=$info{$id}{pos}*$newdis/$max;
		if ($n == 0) {
			$pos0=$info{$id}{pos};
			$n++;
		}
		$info{$id}{pos}=$info{$id}{pos}-$pos0;
		print Out $id,"\t",$info{$id}{pos},"\n";
	}

}
close Out;
@map=glob("$dmap/*.male.map");
open Out,">$dOut/total.male.map";
my %males;

foreach my $map (@map) {
	my $lgID=(split(/\./,basename($map)))[0];
	$lgID=~s/\D+//g;
	my $nloc=`wc -l $map`;
	$nloc=(split(/\s+/,$nloc))[0];
	next if ($nloc <= 2);
	print Out "group\t$lgID\n";
	open In,$map;
	my $max=0;
	my %male;
	while (<In>) {
		chomp;
		next if ($_ eq ""||/^$/ || /^;/ || /group/);
		my ($id,$pos)=split(/\s+/,$_);
		$male{$id}{lgID}=$lgID;
		$male{$id}{pos}=$pos;
		$max=$pos if ($pos > $max);
	}
	close In;
	my $newdis=$max;
	if ($adjust) {
		$newdis=rand(80)+120;
	}
	my $n=0;
	my $pos0;
	foreach my $id (sort {$male{$a}{pos}<=>$male{$b}{pos}} keys %male) {
		next if ($info{$id}{lgID} ne $lgID);
		$male{$id}{pos}=$male{$id}{pos}*$newdis/$max;
		if ($n == 0) {
			$pos0=$male{$id}{pos};
			$n++;
		}
		$male{$id}{pos}=$male{$id}{pos}-$pos0;
		$males{$id}{pos}=$male{$id}{pos};
		print Out $id,"\t",$male{$id}{pos},"\n";
	}

}
close Out;
@map=glob("$dmap/*.female.map");
open Out,">$dOut/total.female.map";
	my %females;

foreach my $map (@map) {
	my $lgID=(split(/\./,basename($map)))[0];
	$lgID=~s/\D+//g;
	my $nloc=`wc -l $map`;
	$nloc=(split(/\s+/,$nloc))[0];
	next if ($nloc <= 2);
	print Out "group\t$lgID\n";
	open In,$map;
	my $max=0;
	my %female;

	while (<In>) {
		chomp;
		next if ($_ eq ""||/^$/ || /^;/ || /group/);
		my ($id,$pos)=split(/\s+/,$_);
		$female{$id}{lgID}=$lgID;
		$female{$id}{pos}=$pos;
		$max=$pos if ($pos > $max);
	}
	close In;
	my $newdis=$max;
	if ($adjust) {
		$newdis=rand(100)+120;
	}
	my $n=0;
	my $pos0;
	foreach my $id (sort {$female{$a}{pos}<=>$female{$b}{pos}} keys %female) {
		next if ($info{$id}{lgID} ne $lgID);
		$female{$id}{pos}=$female{$id}{pos}*$newdis/$max;
		if ($n == 0) {
			$pos0=$female{$id}{pos};
			$n++;
		}
		$female{$id}{pos}=$female{$id}{pos}-$pos0;
		$females{$id}{pos}=$female{$id}{pos};
		print Out $id,"\t",$female{$id}{pos},"\n";
	}
	
}
close Out;

my @marker=glob("$dmap/*.correct.loc");
open Out,">$dOut/total.loc";
my $nloc=0;
my $nind=0;
my @out;
foreach my $marker (@marker) {
	next if ($marker  =~ /pri/);
	open In,$marker;
	while (<In>) {
		chomp;
		next if ($_ eq ""||/^$/);
		my @info=split(/\s+/,$_);
		next if (scalar @info < 3);
		if (/=/) {
			next;
		}else{
			if ($nind != scalar @info -3 && $nind != 0) {
				die "error $_";
			}else{
				$nind=scalar @info-3;
			}
			$nloc++;
			$info[2] = "{00}" if ($info[2] eq "{\-\-}");
			push @out,join("\t",@info);
		}
	}
	close In;
}
open Out,">$dOut/total.sexAver.loc";
print Out join("\n","nloc = $nloc","nind = $nind","popt = CP","name = Total",@out),"\n\n";
#print Out "individual names\n",join("\n",@Indi);
close Out;
open Out,">$dOut/total.final.loc";
print Out join("\n","nloc = $nloc","nind = $nind","popt = CP","name = Total",@out),"\n\n";
print Out "individual names\n",join("\n",@Indi);
close Out;
 @marker=glob("$dmap/*.correct.phase");
 @out=();
open Out,">$dOut/total.sexAver.phase";
open Male,">$dOut/total.male.phase";
open Female,">$dOut/total.female.phase";
print Out join(",","Genotype","","",@Indi),"\n";
print Male join(",","Genotype","","",@Indi),"\n";
print Female join(",","Genotype","","",@Indi),"\n";
foreach my $marker (@marker) {
	open In,$marker;
	while (<In>) {
		chomp;
		next if ($_ eq ""||/^$/);
		my ($id,$type,$phase,@info)=split(/\s+/,$_);
	#	next if ($type eq "hkxhk");
		my @out1;
		my @out2;
		push @out1,join(",",$id,$info{$id}{lgID},$info{$id}{pos});
		push @out2,join(",",$id,$info{$id}{lgID},$info{$id}{pos});
		my @outm;
		if (exists $males{$id}) {
			push @outm,join(",",$id,$info{$id}{lgID},$males{$id}{pos});
		}
		my @outf;
		if (exists $females{$id}) {
			push @outf,join(",",$id,$info{$id}{lgID},$females{$id}{pos});
		}
		for (my $i=0;$i<@info;$i++) {
			my ($p1,$p2)=split(//,$info[$i]);
			if ($type =~ /nnxnp/) {
				if ($p1 eq "-") {
					push @out1,"-";
					push @out2,"-";
				}else{
					push @out1,0;
					push @out2,$p2;
				}
			}elsif($type =~ /lmxll/){
				if ($p1 eq "-") {
					push @out1,"-";
					push @out2,"-";
				}else{
					push @out1,$p1;
					push @out2,0;
				}
			}else{
				push @out1,$p1;
				push @out2,$p2;
			}
		}
		print Out join(",",@out1),"\n";
		print Out join(",",@out2),"\n";
		if (exists $males{$id}) {
			print Male join(",",@outm,@out2[1..$#out2]),"\n";
		}
		if (exists $females{$id}) {
			print Female join(",",@outf,@out1[1..$#out1]),"\n";
		}
	}
	close In;
}
close Out;
close Female;
close Male;
open Out,">$dOut/total.qtl";
print Out join("\n","nloc = $nloc","nind = $nind","popt = CP","name = Total",@out),"\n\n";
#print Out "individual names\n",join("\n",@Indi);
close Out;
open Out,">$dOut/total.trt";
print Out "ntrt=1\n";
print Out "nind=$nind\n";
print Out "miss=*\n";
print Out "Genotype\n";
for (my $i=0;$i<$nind;$i++) {
	print Out "$Indi[$i]\n";
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
sub determineHaploSource {#
	my ($type,$linkPhase,$progenyGeno)=@_;

	return "--" if ($progenyGeno eq "--") ;

	my $haploSource='';;
	my (%parentAllel,@gameteCombination)=();

	my %haploIndex=(
		"0"=>"11",
		"1"=>"12",
		"2"=>"21",
		"3"=>"22",
	);

	my ($p1,$p2,$m1,$m2)= $type =~/(.)(.)x(.)(.)/ ;
	@{$parentAllel{"P"}}=($p1,$p2);
	@{$parentAllel{"M"}}=($m1,$m2);
	my ($PLinkPhase, $MLinkPhase) = split //,$linkPhase ;
	if ($PLinkPhase eq '1'){
		@{$parentAllel{"P"}} = reverse @{$parentAllel{"P"}};
	}
	if ($MLinkPhase eq '1'){
		@{$parentAllel{"M"}} = reverse @{$parentAllel{"M"}} ;
	}

	foreach my $pGamete (@{$parentAllel{"P"}}) {
		foreach my $mGamete (@{$parentAllel{"M"}}) {
			push @gameteCombination,$pGamete.$mGamete;
		}
	}

	for (my $j=0;$j<@gameteCombination ;$j++) {
		
		my @haplo=split //,$gameteCombination[$j];
		my $haplo=join("",sort {$a cmp $b} @haplo);

		if ($haplo eq $progenyGeno) {

			my ($allel1,$allel2) = $progenyGeno =~/(.)(.)/;

			if ($p1 eq $p2) {

				my ($lp) = $haploIndex{$j} =~/.(.)/;
				$haploSource = "0".$lp;
				last;
				
			}elsif ($m1 eq $m2) {

				my ($lp) = $haploIndex{$j} =~/(.)./;
				$haploSource = $lp."0";
				last;

			}elsif (join("",sort {$a cmp $b} ($p1,$p2)) eq join("",sort {$a cmp $b} ($m1,$m2)) and $allel1 ne $allel2) {
				
				$haploSource = "--";
				last;

			}else{
				
				$haploSource = $haploIndex{$j};
			}
		}
		
	}
	return $haploSource;

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
  -marker	<file>	input marker file
  -h         Help

USAGE
        print $usage;
        exit;
}

#!/usr/bin/env perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($map,$fOut);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$map,
	"o:s"=>\$fOut,
			) or &USAGE;
&USAGE unless ($map and $fOut);
open In,$map;
my %Marker;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/ ||/^;/) ;
	if (/group/) {
		next;
	}else{
		my ($id,$female,$male,$sexAver)=split(/\s+/,$_);
		if ($female ne "NA") {
			$Marker{female}{$id}=$female;
		}
		if ($male ne "NA") {
			$Marker{male}{$id}=$male;
		}
		if ($sexAver ne "NA") {
			$Marker{sexAver}{$id}=$sexAver;
		}
	}
}
close In;
open Out,">$fOut.sexAver.map";
my $lgid=basename("$fOut.sexAver.map");
$lgid=~s/\.sexAver\.map//g;
print Out "group $lgid\n";
foreach my $id (sort{$Marker{sexAver}{$a}<=>$Marker{sexAver}{$b}} keys %{$Marker{sexAver}}) {
	print Out $id,"\t",$Marker{sexAver}{$id},"\n";
}
close Out;
open Out,">$fOut.male.map";
print Out "group $lgid\n";
foreach my $id (sort{$Marker{male}{$a}<=>$Marker{male}{$b}} keys %{$Marker{male}}) {
	print Out $id,"\t",$Marker{male}{$id},"\n";
}
close Out;
open Out,">$fOut.female.map";
print Out "group $lgid\n";
foreach my $id (sort{$Marker{female}{$a}<=>$Marker{female}{$b}} keys %{$Marker{female}}) {
	print Out $id,"\t",$Marker{female}{$id},"\n";
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
  -i	<file>	input cross link file 
  -o	<file>	output file
  -h         Help

USAGE
        print $usage;
        exit;
}

#!/usr/bin/env perl
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fIn,$fOut);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$fIn,
	"o:s"=>\$fOut,
			) or &USAGE;
&USAGE unless ($fIn and $fOut);
open In,$fIn;
open Out,">$fOut";
my %info;
while (<In>) {
	chomp;
	next if ($_ eq "" || /^$/ );
	my ($chr,$source,$type,$start,$end,@info)=split(/\t/,$_);
	next if ($type eq "region");
	my $info=join("\t",$chr,$start,$end);
	if ($type ne "exon" && $type ne "CDS") {
		if (!exists $info{$info}) {
			if ($type eq "gene") {
				print Out $_,"\n";
				$info{$info}=1;
			}else{
				print Out join("\t",$chr,$source,"gene",$start,$end,@info),"\n";
				print Out $_,"\n";
				$info{$info}=1;
			}
		}else{
			print Out $_,"\n";
		}
	}else{
		print Out $_,"\n";
	}
}	
close In;
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
  -i	<file>	input dict name
  -o	<file>	out chr list name
  -h         Help

USAGE
        print $usage;
        exit;
}

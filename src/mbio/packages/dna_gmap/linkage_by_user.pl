#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fIn,$fOut,$scheme);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Storable;

my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$fIn,
	"o:s"=>\$fOut,
	"s:s"=>\$scheme,
			) or &USAGE;
&USAGE unless ($fIn and $fOut and $scheme);
my $Tree;
open In,$scheme;
my %lg;
my $read;
while ($read=<In>) {
	chomp $read;
	next if ($read eq ""||$read =~ /^$/);
	$lg{$read}++;
}
close In;
$Tree = retrieve("$fIn");
my %out;
sub_tree($Tree,\%lg,\%out);
open Out,">$fOut";
my $lgid=0;
foreach my $lg (sort keys %lg) {
	if (scalar @{$out{$lg}} == 0) {
		print STDERR "$lg is not exists! plead check !";
		next;
	}
	foreach my $l (@{$out{$lg}}) {
		$lgid++;
		print Out ">LG$lgid\t$lg\n";
		print Out $l,"\n";
	}
}
close Out;
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub sub_tree{
	my($tree,$lg,$out)=@_;
	my @queue=();
	push @queue,$tree->{'root'};
	for(my $i=0;$i<@queue;$i++) {
		my $str=$queue[$i]->{'lod'}. "/".$queue[$i]->{'nloc'};
		if (exists $$lg{$str}) {
			push @{$$out{$str}},join("\t",@{$queue[$i]->{'locus'}});
		}
		if (defined $queue[$i]->{'child'} && scalar @{$queue[$i]->{'child'}}!=0) {
			if ($i >= scalar @queue -1) {
				@queue = (@queue[0..$i],@{$queue[$i]->{'child'}});
			}else{
				@queue = (@queue[0..$i],@{$queue[$i]->{'child'}},@queue[$i+1..$#queue]);
			}
		}
	}
	close Out;
}

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
  -i	<file>	input file name
  -o	<file>	output file
  -s	<file>	group scheme list
  -h         Help

USAGE
        print $usage;
        exit;
}

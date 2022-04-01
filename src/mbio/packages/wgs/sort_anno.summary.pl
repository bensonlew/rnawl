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
my %chr;
my %hash;
$fIn=ABSOLUTE_DIR($fIn);
my($gene,$trans,$chr,$start,$end);
open Out,">$fOut";
while (<In>) {
	chomp;
    print $_,"\n";
    next if($_ eq "" || /^$/);
	if ($_ =~ /^#/){print Out "$_\n";next;}
    my $mm = $_;
	my@a=split(/\|/,(split(/\t/,$_))[0]);
    next if($a[0] eq ':--');
	my @b=split(/\:/,$a[0]);
	$gene=$b[0] ;
	next if($gene eq "");
	$trans=$b[1];
	my @c=split(/\:/,$a[-1]);
	$chr=$c[-3];
	$start=$c[-2];
	$end=$c[-1];
	my $nr;
	$trans=$gene if($trans eq "--");
    $trans=~/([a-zA-Z]*)([0-9]*)/;
    die "\nchr/scaID check ERROR" if($1 eq "" and $2 eq "");
    $hash{$1}{$2}{$chr}{$start}{$end}=$mm;
}
foreach my $n (sort keys %hash){
	foreach my $m (sort {$a<=>$b} keys %{$hash{$n}}){
		foreach my $i (sort keys %{$hash{$n}{$m}}){
        	foreach my $j (sort {$a<=>$b} keys %{$hash{$n}{$m}{$i}}){
            	foreach my $k (sort {$a<=>$b} keys %{$hash{$n}{$m}{$i}{$j}}){
                	print Out $hash{$n}{$m}{$i}{$j}{$k},"\n";
				}
			}
		}
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
Contact:        qingmei.cui\@majorbio.com;
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

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

my %hash;
$fIn=ABSOLUTE_DIR($fIn);
open Out,">$fOut";
while (<In>) {
	chomp;
    next if($_ eq "" || /^$/);
	if ($_ =~ /^#/){print Out "$_\n";next;}
    my $mm = $_;
	my($a,@arr)=split(/\t/,$_);
    my @b = split/_/,$a;
	$b[0]=~/([a-z]*)([0-9]*)/;
    die "chr/scaID check ERROR" if($1 eq "" and $2 eq "");
    $hash{$1}{$2}{$b[1]}{$mm}=1;
}
foreach my $n (sort keys %hash){
        foreach my $m (sort {$a<=>$b} keys %{$hash{$n}}){
			foreach my $pos (sort {$a<=>$b} keys %{$hash{$n}{$m}}){
				foreach my $mm (keys %{$hash{$n}{$m}{$pos}}){
					print Out $mm,"\n";}
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

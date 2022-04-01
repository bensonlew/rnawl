#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fIn,$fOut,$select);
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
if ($fIn =~ /\.gz$/) {
	close In;
	open In,"zcat $fIn|";
}
my @Indi;
my %Stat;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/ ||/^##/);
	if (/^#/) {
		my ($chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$format,@indi)=split(/\t/,$_);
		@Indi=@indi;
	}else{
		my ($chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$format,@sample)=split(/\t/,$_);
		my @alt=split(/,/,join(",",$ref,$alt));
		my %ale;
		my %len;
		for (my $i=0;$i<@alt;$i++) {
			$ale{$alt[$i]}=$i;
			$len{length($alt[$i])}=1;
		}
		my $snptype="SNP";
		if (scalar keys %len > 1) {
			$snptype="INDEL";
		}
		my %ann;
		if ($info =~/ANN=([^;]*)/) {
			my $ann=$1;
			my @ann=split(/\,/,$ann);
			for (my $i=0;$i<@ann;$i++) {
				my ($ale,$eff,$type)=split(/\|/,$ann[$i]);
				my @eff=split(/\&/,$eff);
				foreach my $eff (@eff) {
					$ann{$ale{$ale}}{$eff}++;
				}
				$ann{$ale{$ale}}{$type}++;
			}
		}
		my @format=split(/:/,$format);
		for (my $i=0;$i<@Indi;$i++) {
			my @info=split(/:/,$sample[$i]);
			for (my $j=0;$j<@info;$j++) {
				if ($format[$j] eq "GT") {
					next if ($info[$j] eq "./." || $info[$j] eq "0/0");
					my ($g1,$g2)=split(/\//,$info[$j]);
					if ($g1 eq $g2) {
						foreach my $type (sort keys %{$ann{$g1}}) {
							$Stat{$snptype}{$type}{$Indi[$i]}++;
						}
					}else{
						foreach my $type (sort keys %{$ann{$g1}}) {
							$Stat{$snptype}{$type}{$Indi[$i]}++;
						}
						foreach my $type (sort keys %{$ann{$g2}}) {
							$Stat{$snptype}{$type}{$Indi[$i]}++;
						}
					}
				}
			}
		}
	}
}
close In;
mkdir $fOut if (!-d $fOut);
open Out,">$fOut/snp.stat";
print Out join("\t","sampleID",sort keys %{$Stat{"SNP"}}),"\n";
foreach my $Indi (@Indi) {
	my @out;
	push @out,$Indi;
	foreach my $type (sort keys %{$Stat{"SNP"}}) {
		$Stat{"SNP"}{$type}{$Indi}||=0;
		push @out,$Stat{"SNP"}{$type}{$Indi};
	}
	print Out join("\t",@out),"\n";
}
close Out;
open Out,">$fOut/indel.stat";
print Out join("\t","sampleID",sort keys %{$Stat{"INDEL"}}),"\n";
foreach my $Indi (@Indi) {
	my @out;
	push @out,$Indi;
	foreach my $type (sort keys %{$Stat{"INDEL"}}) {
		$Stat{"INDEL"}{$type}{$Indi}||=0;
		push @out,$Stat{"INDEL"}{$type}{$Indi};
	}
	print Out join("\t",@out),"\n";
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
  -i	<file>	input vcf name
  -o	<file>	output file name
  -h         Help

USAGE
        print $usage;
        exit;
}

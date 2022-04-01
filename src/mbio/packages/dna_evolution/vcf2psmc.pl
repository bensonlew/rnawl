#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fIn,$fOut,$pop);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"i:s"=>\$fIn,
	"p:s"=>\$pop,
	"o:s"=>\$fOut
			) or &USAGE;
&USAGE unless ($fIn and $fOut);
mkdir $fOut if (!-d $fOut);
my %pop;
open In,$pop;
while (<In>) {
	chomp;
	next if ($_ eq "" || /^$/ ||/^#/);
	my ($id,$popid)=split(/\s+/,$_);
	$pop{$id}=$popid;
}
close In;
open In,$fIn;
my %seq;
my @Indi;
my %BASE=(
	"AA"=>"A","GG"=>"G","CC"=>"C","TT"=>"T",
	"AT"=>"W","AG"=>"R","AC"=>"M",
	"CG"=>"S","CT"=>"Y",
	"GT"=>"K"
);
while (<In>) {
	chomp;
	next if ($_ eq "" || /^$/ ||/^##/);
	if (/^#/) {
		(undef,undef,undef,undef,undef,undef,undef,undef,undef,@Indi)=split(/\t/,$_);
	}else{
		my ($CHROM,$POS,$ID,$REF,$ALT,$QUAL,$FILTER,$INFO,$FORMAT,@indi)=split(/\t/,$_);
		my @ale=split(/\,/,join(",",$REF,$ALT));
		for (my $i=0;$i<@indi;$i++) {
			my $geno=(split(/\:/,$indi[$i]))[0];
			my ($g1,$g2)=split(/\//,$geno);
			if ($g1 eq ".") {
					if (!exists $seq{$Indi[$i]}{$CHROM}) {
						$seq{$Indi[$i]}{$CHROM}="N";
					}else{
						$seq{$Indi[$i]}{$CHROM}.="N";
					}
			}else{
				if (!exists $BASE{join("",sort($ale[$g1],$ale[$g2]))}){
					if (!exists $seq{$Indi[$i]}{$CHROM}) {
						$seq{$Indi[$i]}{$CHROM}="N";
					}else{
						$seq{$Indi[$i]}{$CHROM}.="N";
					}
					next;
				};
				if (!exists $seq{$Indi[$i]}{$CHROM}) {
					$seq{$Indi[$i]}{$CHROM}=$BASE{join("",sort($ale[$g1],$ale[$g2]))};
				}else{
					$seq{$Indi[$i]}{$CHROM}.=$BASE{join("",sort($ale[$g1],$ale[$g2]))};
				}
			}
		}
	}
}
close In;
my $head=0;
open SH,">$fOut/work1.sh";
my %poppsmc;
foreach my $id (sort keys %seq) {
	print SH "psmc -N25 -t15 -r5 -p \"4+25*2+4+6\" -o $fOut/$id.psmc $fOut/$id.fasta\n";
	push @{$poppsmc{$pop{$id}}},"$fOut/$id.psmc";
	open FA,">$fOut/$id.fasta";
	foreach my $chr (sort keys %{$seq{$id}}) {
		print FA ">$chr\n$seq{$id}{$chr}\n";
	}
	close FA;
}
close SH;
open SH,">$fOut/work2.sh";
my @id;
my @psmcfile;
foreach my $popid (sort keys %poppsmc) {
	print SH "cat ".join(" ",@{$poppsmc{$popid}}) ."> $fOut/$popid.psmc && ";
	print SH "psmc_plot.pl -R -p $fOut/$popid $fOut/$popid.psmc\n";
	push @id,$popid;
	push @psmcfile,"$fOut/$popid.psmc";
}
print SH "psmc_plot.pl -R -M \"".join(",",@id)."\" -p $fOut/combined.psmc ".join(" ",@psmcfile),"\n";
close SH;
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
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
  -o	<file>	output dir
  -p	<file>	population file
  -h         Help

USAGE
        print $usage;
        exit;
}

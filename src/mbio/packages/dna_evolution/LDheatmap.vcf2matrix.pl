#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fin,$fout,$min);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	#"i:s"=>\$fIn,
	#"i2:s"=>\$fIn,
	"i:s"=>\$fin,
	"o:s"=>\$fout,
	"m:s"=>\$min,
			) or &USAGE;
&USAGE unless ($fout);
open IN,$fin;
my @indi;
my %seq;
my %head;
while (<IN>) {
	chomp;
	next if ($_ eq "" || /^$/ || /^##/);
	my ($chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$format,$sam)=split/\s+/,$_,10;
	if (/^#/) {
		@indi = split/\s+/,$sam;
		foreach my $sample (@indi) {
			$seq{$sample}{seq}="";
			$head{id}="";
			$head{pos}="";
		}
	}else{
		my $ID = $chr."_".$pos;
		my @gale = split/\,/,join(",",$ref,$alt);
		my $len =scalar @gale;
		next if ($len ne 2);
		my @geno = split/\s+/,$sam;
		my @format = split/\:/,$format;
		my @id;
		my @pos;
		for (my $i=0;$i<@indi ;$i++) {
			#$seq{$indi[$i]}{pos}.="$pos\t";
			my @ale = split(/\:/,$geno[$i]);
			for (my $j=0;$j<@format;$j++) {
				if ($format[$j] eq "GT") {
					my ($g1,$g2) = split/\//,$ale[$j];
					#print $g1,"\t",$g2;die;
					if ($g1 eq ".") {
						$seq{$indi[$i]}{seq}.="<NA>\t";
					}elsif ($g1 eq 0 && $g2 eq 0) {
						$seq{$indi[$i]}{seq}.="$ref/$ref\t";
					}elsif ($g1 eq 0 && $g2 eq 1) {
						$seq{$indi[$i]}{seq}.="$ref/$alt\t";
					}elsif ($g1 eq 1 && $g2 eq 1) {
						$seq{$indi[$i]}{seq}.="$alt/$alt\t";
					}else{
						$seq{$indi[$i]}{seq}.="<NA>\t";
					}
				}
			}
		}
		$head{pos}.="$pos\t";
		$head{id}.="$ID\t";
	}
}
close IN;
#print Dumper \%seq;die;
open OUT,">$fout";
print OUT "$head{id}\n";
print OUT "$head{pos}\n";
foreach my $seq (sort keys %seq) {
	print OUT "$seq{$seq}{seq}\n";
}
close OUT;
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub USAGE {#
        my $usage=<<"USAGE";
Contact:        minghao.zhang\@majorbio.com;
Script:			$Script
Description:
	fq thanslate to fa format
	eg:
	perl $Script -i -o -k -c

Usage:
  Options:
	"i:s"=>\$fin,                                                                                                                                                                        
    "o:s"=>\$fout,
	"m:s"=>\$min,
  -h         Help

USAGE
        print $usage;
        exit;
}

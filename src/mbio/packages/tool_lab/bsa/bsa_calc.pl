#!/usr/bin/env perl 
use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Data::Dumper;
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fIn,$fOut,$P1,$P2,$B1,$B2,$popt);
GetOptions(
				"help|?" =>\&USAGE,
				"table:s"=>\$fIn,
				"out:s"=>\$fOut,
				"wp:s"=>\$P1,
				"mp:s"=>\$P2,
				"wb:s"=>\$B1,
				"mb:s"=>\$B2,
				"popt:s"=>\$popt,
				) or &USAGE;
&USAGE unless ($fIn and $fOut and $B2);
$popt ||="F2";
open In,$fIn;
open Out,">$fOut";
if (!defined $B1) {
	print Out "#chr\tpos\tdelta\tdepth\n";
}else{
	print Out "#chr\tpos\tindex1\tindex2\tdelta\tn1\tn2\tn3\tn4\n";
}
my @head;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/);
	if (/#/ ||/CHROM/) {
		my ($chr,$pos,$ref,$alt,$type,@sample)=split(/\s+/,$_);
		push @head,@sample;
		pop(@head);
	}else{
		my ($chr,$pos,$ref,$alt,$type,@sample)=split(/\s+/,$_);
		my %info;
		for (my $i=0;$i<@head;$i++) {
			my ($id,$info)=(split(/\./,$head[$i]));
			$info{$id}{$info}=$sample[$i];
		}
		my %ad;
		foreach my $sample (keys %info) {
			my @ad=split(/\,/,$info{$sample}{AD});
			for (my $i=0;$i<@ad;$i++) {
				$ad{$sample}{$ref}=$ad[0];
				$ad{$sample}{$alt}=$ad[1];
			}
		}
		if ($popt eq "F2") {
			if (!defined $B1) {
				my $index;
				if (defined $P1 && !defined $P2) {
					my @gtP=split(/\//,$info{$P1}{GT});
					if ($gtP[0] eq $ref) {
						$index=$ad{$B2}{$alt}/$info{$B2}{DP};
					}else{
						$index=$ad{$B2}{$ref}/$info{$B2}{DP};
					}
				}elsif (defined $P2 && !defined $P1) {
					my @gtP=split(/\//,$info{$P2}{GT});
					$index=$ad{$B2}{$gtP[0]}/$info{$B2}{DP};
				}elsif (defined $P1 && defined $P2) {
					my @gtP=split(/\//,$info{$P2}{GT});
					$index=$ad{$B2}{$gtP[0]}/$info{$B2}{DP};
				}elsif (!defined $P2 && !defined $P1) {
					$index=$ad{$B2}{$alt}/$info{$B2}{DP};
				}
				print Out join("\t",$chr,$pos,$index,$info{$B2}{DP}),"\n";
			}else{
				my ($index1,$index2,$delta);
				if (defined $P1 && !defined $P2) {
					my @gtP=split(/\//,$info{$P1}{GT});
					my $numerator=$alt if ($gtP[0] eq $ref);
					$index1=$ad{$B2}{$numerator}/$info{$B2}{DP};
					$index2=$ad{$B1}{$numerator}/$info{$B1}{DP};
					$delta=$index1-$index2;
				}elsif (defined $P2 && !defined $P1) {
					my @gtP=split(/\//,$info{$P2}{GT});
					my $numerator=$gtP[0];
					$index1=$ad{$B2}{$numerator}/$info{$B2}{DP};
					$index2=$ad{$B1}{$numerator}/$info{$B1}{DP};
					$delta=$index1-$index2;
				}elsif (defined $P1 && defined $P2) {
					my @gtP=split(/\//,$info{$P2}{GT});
					my $numerator=$gtP[0];
					$index1=$ad{$B2}{$numerator}/$info{$B2}{DP};
					$index2=$ad{$B1}{$numerator}/$info{$B1}{DP};
					$delta=$index1-$index2;
				}elsif (!defined $P2 && !defined $P1) {
					$index1=$ad{$B2}{$alt}/$info{$B2}{DP};
					$index2=$ad{$B1}{$alt}/$info{$B1}{DP};
					$delta=abs($index1-$index2);
				}
				print Out join("\t",$chr,$pos,$index1,$index2,$delta,$ad{$B1}{$ref},$ad{$B2}{$ref},$ad{$B1}{$alt},$ad{$B2}{$alt}),"\n";
			}
		}else{#F1
			if (!defined $B1) {
				my $index=$ad{$B2}{$alt}/$info{$B2}{DP};
				print Out join("\t",$chr,$pos,$index,$info{$B2}{DP}),"\n";
			}else{
				my ($index1,$index2,$delta);
				$index1=$ad{$B2}{$alt}/$info{$B2}{DP};
				$index2=$ad{$B1}{$alt}/$info{$B1}{DP};
				$delta=abs($index1-$index2);
				print Out join("\t",$chr,$pos,$index1,$index2,$delta,$ad{$B1}{$ref},$ad{$B2}{$ref},$ad{$B1}{$alt},$ad{$B2}{$alt}),"\n";
			}
		}
	}
}
close Out;
close In;
#######################################################################################
print "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

sub USAGE {#
	my $usage=<<"USAGE";
Description: 
Version:  $Script
Contact: long.huang

Usage:
  Options:

	-table input table for calculate
	-out	output filename for result
	-wp	wild parent
	-mp	mutant parend
	-wb	wild bulk
	-mb	mut bulk
	-popt	population type

USAGE
	print $usage;
	exit;
}

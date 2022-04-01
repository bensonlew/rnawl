#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($vcf,$mark,$fOut);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"vcf:s"=>\$vcf,
	"mark:s"=>\$mark,
	"out:s"=>\$fOut
			) or &USAGE;
&USAGE unless ($vcf and $mark);
my %mark;
my (@sample,@marker);
open In,$mark;
while (<In>){
	next if ($_ eq ""|| /^$/);
	chomp;
	if($_=~/^#Mark/){
		(undef,undef,@sample)=split(/\t/,$_);
	}else{
		my($id,$type,@undi)=split(/\t/,$_);
		my $marker=$id;
		$mark{$marker}=$type;
		push @marker,$id;
	}
}
close In;

my (%stat,%dep,%per);
open In,$vcf;
if($vcf=~/gz$/){
	close In;
	open In,"gunzip -c $vcf|";
}
my @Indi;
while (<In>) {
	chomp;
	next if ($_ eq ""|| /^$/ || /^##/);
	if (/^#/) {
		my ($chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$format,@indi)=split(/\t/,$_);
		@Indi=@indi;
	}else{
		my ($chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$format,@samples)=split(/\t/,$_);
		my @format=split(/\:/,$format);
		my $markerdep=0;
		#foreach my $marker(@marker){
		foreach my $marker(keys %mark){
			if($id eq $marker){
				for (my $i=0;$i<@Indi;$i++) {
					my @info=split(/\:/,$samples[$i]);
					my $samples=$Indi[$i];
					foreach my $sample(@sample){
						if($samples eq $sample){
							for (my $j=0;$j<@info;$j++) {
								if ($format[$j] eq "DP") {
									$per{$sample}{$marker}=$info[$j];
									$stat{$marker}{$sample}{dep}=$info[$j];
									$markerdep=$markerdep + $info[$j];
								}
							}
						}
					}
				}
			my $mark=$marker;
			$dep{$markerdep}{$mark}=$mark{$marker};
			}
		}
	}
}                      
close In;

my(@out1,@out2,@out3);
my $marknum=scalar@marker ;
my %persam;
foreach my $sample(keys %per){
	my $samdep=0;
	foreach my $marker(keys %{$per{$sample}}){
		$samdep=$samdep + $per{$sample}{$marker};
	}
	my $depth=sprintf("%.2f",$samdep / (scalar@marker)) ;
	$persam{$depth}{$sample}=$samdep;
}

foreach my $depth(sort{$b<=>$a}keys%persam){
	foreach my $sample(keys %{$persam{$depth}}){
		push @out1,$sample;
		push @out2,$depth;
	}
}
open Out,">$fOut";
print Out join("\t","#SampleID","","",@out1),"\n",join("\t","#MarkerID","Type","Depth",@out2),"\n";

foreach my $markerdep (sort{$b<=>$a} keys %dep){
	foreach my $mark (keys %{$dep{$markerdep}}){
		my @info=();
		for(my $i=0;$i<@out1;$i++){
			foreach my $marker(%stat){
				if($marker eq $mark){
					foreach my $sample(%{$stat{$marker}}){
						if($sample eq $out1[$i]){
							push @info,$stat{$marker}{$sample}{dep};
						}
					}
				}
			}
		}
		push @out3,join("\t",$mark,$dep{$markerdep}{$mark},$markerdep,join("\t",@info));
	}	
}
print Out join("\n",@out3);
close Out;

#######################################################################################
#print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub USAGE {#
        my $usage=<<"USAGE";
Contact:        caixia.tian\@majorbio.com;
Script:			$Script
Description:
	fq thanslate to fa format
	eg:
	perl $Script -i -o -k -c

Usage:
  Options:
  -vcf	<file>	input pop.vcf file
  -mark	<file>  input filtered.marker
  -out	<file>	output result file
  -h         Help

USAGE
        print $usage;
        exit;
}

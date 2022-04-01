#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fIn,$mark,$x,$y,$fOut);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"depth:s"=>\$fIn,
	"mark:s"=>\$mark,
	"x:s"=>\$x,
	"y:s"=>\$y,
	"out:s"=>\$fOut
			) or &USAGE;
&USAGE unless ($fIn and $mark and $fOut);

my ($sample1,$sample2)=split(/\,/,$x);
my ($start,$end)=split(/\,/,$y);
my (@tsample,@tmarker);
my %info;
open In,$fIn;
while (<In>){
	next if ($_ eq ""|| /^$/|| /^#MarkerID/);
	chomp;
	if($_=~/SampleID/){
		(undef,undef,@tsample)=split(/\t/,$_);
	}else{
		my($id,$type,$tdepth,@pdepth)=split(/\t/,$_);
		my $markerid=$id;
		$info{$markerid}=$type;
		push @tmarker,$id;
	}
}
close In;

my $sample=join("\t",@tsample);
$sample=(split(/$sample1/,$sample))[1];
$sample=(split(/$sample2/,$sample))[0];
$sample=$sample1.$sample.$sample2 ;

my $marker=join("\t",@tmarker);
$marker=(split(/$start/,$marker))[1];
$marker=(split(/$end/,$marker))[0];
$marker=$start.$marker.$end ;

#print $sample,"\n",$marker,"\n";
my @sample=split(/\t/,$sample);
my @marker=split(/\t/,$marker);

my %stat;
my @Indi;
open In,$mark;
while (<In>) {
	chomp;
	next if ($_ eq ""|| /^$/);
	if (/^#MarkerID/) {
		(undef,@Indi)=split(/\t/,$_);
	}else{
		my ($id,@indi)=split(/\t/,$_);
		for (my $i=0;$i<@Indi;$i++) {
			my $samplen=$Indi[$i];
			$stat{$id}{$samplen}{type}=$indi[$i];
		}
	}
}
close In;

my (@samples,@markers);
open Out,">$fOut";
print Out "#MarkID\tType\t",join("\t",@sample),"\n";
foreach my $marker(@marker){
	foreach my $markerid(keys %info){
		if($marker eq $markerid){
			print Out $marker,"\t",$info{$markerid},"\t";
			foreach my $id(keys %stat){
				if($marker eq $id ){
					#print Out $marker,"\t",$info{$markerid},"\t";
					for(my$j=0;$j<scalar@sample;$j++){
						foreach my $samplen(keys %{$stat{$id}}){
							if($samplen eq $sample[$j]){
								#print Out $sample[$j],"\,",$stat{$id}{$samplen}{type},"\t";
								print Out $stat{$id}{$samplen}{type},"\t";
							}
						}
					}
				}
			}
			print Out "\n";
		}
	}
}
close Out;

#######################################################################################
#print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub USAGE {#
        my $usage=<<"USAGE";
Contact:        chongqing.shi\@majorbio.com;
Script:			$Script
Description:
	fq thanslate to fa format
	eg:
	perl $Script -i -o -k -c

Usage:
  Options:
  -mark	<file>	input total.marker
  -depth	<file>	depth matrix	
  -x	<str>  sample1,sample2
  -y	<str>	marker1,marker2
  -out	<file>	output result file
  -h         Help

USAGE
        print $usage;
        exit;
}

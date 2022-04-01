#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($fIn,$group,$fOut);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
    "g:s"=>\$group,
	"i:s"=>\$fIn,
	"o:s"=>\$fOut,
			) or &USAGE;
&USAGE unless ($fIn and $group and $fOut);
my %groups;
my %list;
open IN,$group;
open OUT,">$fOut/fa.list";
while(<IN>){
    $_=~s/[\n\r]//g;
    my ($samples,$gro)=split/\t/,$_;
    $groups{$samples}=$gro;
    $list{$gro}=1;
    }
close IN;
foreach my $fa (sort keys %list){
    print OUT "$fa\t$fOut\/$fa.fa\n";
    }
close OUT;
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
		foreach my $id (@Indi) {
            if(exists $groups{$id}){
               $seq{$groups{$id}}{$id}="";
            }
		}
	}else{
		my ($CHROM,$POS,$ID,$REF,$ALT,$QUAL,$FILTER,$INFO,$FORMAT,@indi)=split(/\t/,$_);
		my @ale=split(/\,/,join(",",$REF,$ALT));   
		for (my $i=0;$i<@indi;$i++) {
            if(exists $groups{$Indi[$i]}){
                my $geno=(split(/\:/,$indi[$i]))[0];
                my ($g1,$g2)=split(/\//,$geno);
                if ($g1 eq ".") {
                    $seq{$groups{$Indi[$i]}}{$Indi[$i]}.="N";
			}else{
				if (!exists $BASE{join("",sort($ale[$g1],$ale[$g2]))}){
					$seq{$groups{$Indi[$i]}}{$Indi[$i]}.="N";
					next;
				}
				if (!exists $seq{$groups{$Indi[$i]}}{$Indi[$i]}) {
					print $Indi[$i],"\n";die;
				}
				$seq{$groups{$Indi[$i]}}{$Indi[$i]}.=$BASE{join("",sort($ale[$g1],$ale[$g2]))};
			}
		}
       }
   }
}
close In;

foreach my $id (sort keys %seq) {
    open OUT,">$fOut/$id.fa";
    foreach my $seqs (keys %{$seq{$id}}){
        print OUT  ">$seqs\n$seq{$id}{$seqs}\n";
    }
    close OUT;
}
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub USAGE {#
        my $usage=<<"USAGE";
Contact:        nanshan.yang\@majorbio.com;
Script:			$Script
Description:
	
	eg:
	perl $Script -g -i -o 

Usage:
  Options:
  -g    <file>  input group list
  -i	<file>	input file name
  -o	<dir>	
  -h         Help

USAGE
        print $usage;
        exit;
}

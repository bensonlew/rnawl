#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
my %opts;
GetOptions (\%opts,"list=s", "o=s", "s=s");
my $usage = <<"USAGE";
        Program : $0
	Discription:combine idba contigs into one file
	Usage:perl $0 [options]
		-list	contig list
                -s  put sample_name into file,type Y or N,default: "N"
		-o	output_file name
	example:combine_contig.pl -list contig.list -s N -o output/contig.fa
USAGE
$opts{s} = "N" if (!$opts{s});
die $usage if (!$opts{list});
die $usage if (!$opts{o});
die $usage if !($opts{s} eq "Y" or $opts{s} eq "N");

my $i = 0;
open LIST, "$opts{list}" or die "Couldn't find file $opts{list}!\n";
open OUT, ">$opts{o}" or die "Couldn't open file $opts{o}!\n";
while (<LIST>) {
	chomp;
    my @tmp_list = split '\t',$_;
    my $contig_file = $tmp_list[0];
    my $sample_name = 0;
    if ($#tmp_list == 1) {
        $sample_name = $tmp_list[1];
    }else {
        die "No sample name in contig list!\n"if ($opts{s} eq "Y");
    }
	open FA, "$contig_file" or die "Couldn't find file $_!\n";
	while (<FA>) {
	    if (/^>/) {
		$i ++;
		my @tmp = split " ",$_;
                $tmp[0] =~ s/>//;
                if ($opts{s} eq 'N'){
		    print OUT ">contig_$i $tmp[1] $tmp[2]\n";
                } else {
                    print OUT ">$sample_name".'_'."$tmp[0] $tmp[1] $tmp[2]\n";
                }
	    }else {
			print OUT $_;
	    }
	}
	close FA;
}
close LIST;
close OUT;

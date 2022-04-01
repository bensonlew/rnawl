#!/usr/bin/perl -w
use strict;
use warnings;
use Carp;


use Getopt::Long qw(:config no_ignore_case bundling);

my %opts;
GetOptions(\%opts,"i=s","o=s","type=i","gene_type=s","thick_type=s");


my $usage = <<__EOUSAGE__;

#################################################################################### 
# Usage:
#	perl $0 --i ref.gtf --o ref.bed
# Required:
#
#  --i		gtf file
#
####################################################################################


__EOUSAGE__

    ;


die $usage  if( !$opts{i} );


open (INFILE,$opts{i}) || die "Can not open $opts{i}\n";

my %exon = ();
my $exon_num = 1;

while (<INFILE>) {
    my $org_line = $_;
    my @line = split(/\t/,$_);
    if (! ($line[8]) ){
	next;
    }
    if ( $line[8] =~ /transcript_id "([^"]*)"/){
	my $id = $1;
	if($line[2] eq "exon"){
	    $exon{$id}{$exon_num}{start} = $line[3];
	    $exon{$id}{$exon_num}{end} = $line[4];
            $exon_num ++ ;
        }
    }
}

close INFILE;

open (INFILE, $opts{i}) || die "Can not open $opts{i}\n";
while (<INFILE>) {
    my $org_line = $_;
    my @line = split(/\t/,$_);
    if (! ($line[8]) ){
        print $org_line;
	next;
    }
    if ( $line[8] =~ /transcript_id "([^"]*)"/){
	my $id = $1;
	if($line[2] eq "exon"){
	    print $org_line;
        }elsif($line[2] eq "cds" || $line[2] eq "CDS"){
	    my $start = $line[3];
	    my $end = $line[4];
	    my  $has_exon = 0;
	    foreach(keys %{$exon{$id}}){
		my $num = $_;
		if ($start >= $exon{$id}{$num}{start} && $end <= $exon{$id}{$num}{end}){
		    $has_exon = 1;
		    last;
                }
            }
	    if ($has_exon == 0){
		$line[2] = "exon";
                print join("\t", @line);
	    }
	    print $org_line;
	}else{
            print $org_line;
        }
    }
}


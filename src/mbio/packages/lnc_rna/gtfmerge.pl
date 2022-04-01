#!/usr/bin/perl -w
use strict;
use warnings;
use Carp;
use Getopt::Std;
use Getopt::Long;



my %opts;

GetOptions(\%opts,"i=s","o=s","ref=s","compare=s","add!");

my $usage = <<__EOUSAGE__;

#################################################################################### 
# Usage:
#	perl $0 -ref ref.gtf -i cufflink.gtf -o merge.gtf -compare compare.id
# Required:
#
#  -i	    gtf file
#  -o	    gtf file        
#  -ref     ref gtf files
#  -compare cuffcmp.merged.gtf.tmap
#  -add     add ref transcripts 
#
####################################################################################


__EOUSAGE__

;

die $usage if( !$opts{ref} || !$opts{i} || !$opts{compare} );
my  $output= $opts{o}?$opts{o}:$opts{i}.".merge.gtf";



my %ref = ();
my %mapping = ();
my %ref_trans2gene = ();

open (MAPPING,$opts{compare}) || die "Can not open $opts{compare}\n";
<MAPPING>;
while(<MAPPING>){
	chomp;
	my @line = split(/\t/, $_);
	$mapping{$line[4]}{type} = $line[2];
	$mapping{$line[4]}{ref_gene_id} = $line[0];
	$mapping{$line[4]}{ref_id} = $line[1];	
}
close MAPPING;


open (REF,$opts{ref}) || die "Can not open $opts{ref}\n";

my $line_num = 1;
while(<REF>){
	chomp;
	my $line_all = $_;
	my @line = split(/\t/, $line_all);	
	my ($transcript, $gene) = ("", "");
	if(! exists $line[8]){
		next;
	}
	if($line[2] eq "gene" || $line[2] eq "transcript"){
		next;
	}
	if( $line[8] =~ /transcript_id "([^"]*)"/ ){
		$transcript = $1;
	}
	if( $line[8] =~ /gene_id "([^"]*)"/ ){
		$gene = $1;
	}else{
		$gene = $transcript;
	}
	$ref{$transcript}{$line_num} = $line_all;
	$ref_trans2gene{$transcript} = $gene;
	$line_num ++;
}

close REF;

open (INFILE,$opts{i}) || die "Can not open $opts{i}\n";
open (OUT,"> $output") || die "Can not create outfile\n";


my %trans_out =();
my %known_out =();

while (<INFILE>) {
	chomp;
	my $line_all = $_;
	my @line = split(/\t/,$_);
	my ($transcript, $gene) = ("", "");
	if( ! exists $line[8]){
		next;
	}
	if( $line[8] =~ /transcript_id "([^"]*)"/ ){
		$transcript = $1;
	}
	if( $line[8] =~ /gene_id "([^"]*)"/ ){
		$gene = $1;
	}
	if(exists $mapping{$transcript}){
		if($mapping{$transcript}{type} eq "="){
			if(exists $trans_out{$transcript}){
			}else{
				my $transcript_old = $mapping{$transcript}{ref_id};
				#print $transcript_old."***\n";

				my %tran_gtf =  %{$ref{$transcript_old}};
				$trans_out{$transcript_old} = 1;
				if( exists $trans_out{$transcript_old} ){
					next;
				}
				
				foreach(keys %tran_gtf){
					print OUT $tran_gtf{$_}." class_code \"=\";\n";
				}
			}
		}else{
			if( exists $mapping{$transcript}{ref_id} && $mapping{$transcript}{ref_id} ne "-" && exists $ref_trans2gene{$mapping{$transcript}{ref_id}} ){
				$gene = $ref_trans2gene{$mapping{$transcript}{ref_id}};
				$line_all =~ s/gene_id "([^"]*)"/gene_id "$gene"/g;
				print OUT $line_all."\n";
				$trans_out{$transcript} = 1;
			}else{
			    print OUT $line_all."\n";
			}		
		}
	}	
}

if( $opts{add} ){
	foreach(keys %ref){
		my  $transcript = $_;
		if(exists $trans_out{$transcript}){
		}else{
			my %tran_gtf =  %{$ref{$transcript}};
			$trans_out{$transcript} = 1;
			foreach(keys %tran_gtf){
				print OUT $tran_gtf{$_}."\n";
			}
		}
	}
}

close INFILE;
close OUT;

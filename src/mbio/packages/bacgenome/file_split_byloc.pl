#! /usr/bin/perl
use strict;
use warnings;

if(@ARGV!=3) {
    print STDERR "extract_seq_bygenenu.pl seq_file sample outputdir\n";
    exit;
}
my ($seq_file,$sample,$outputdir)=@ARGV;

my ($nanalysis,$location,%hash);
if ((split /\//,$seq_file)[-1] =~ /^CDS.*/){
	$nanalysis = "_CDS.";
	}elsif ((split /\//,$seq_file)[-1] =~ /^TRF.*/){
	$nanalysis = "_TRF.";
	}elsif ((split /\//,$seq_file)[-1] =~ /^rRNA.*/){
	$nanalysis = "_rRNA.";
	}elsif ((split /\//,$seq_file)[-1] =~ /^tRNA.*/){
	$nanalysis = "_tRNA.";
}
my $type = (split /\./,$seq_file)[-1];
open(FNN,"<$seq_file") or die;
if ($type  =~ /faa/ or $type  =~ /fnn/){
	$/ = ">",<FNN>;
	while(<FNN>){
		$location = (split /_/,(split / /,$_)[3])[0];
		$_ =~  s/>//;
		#$location =~ s/^Chr/chromosome/;
		#$location =~ s/^p/plasmid/;
		if (not exists $hash{$location}){
			$hash{$location} = ">".$_;
		}else{
			$hash{$location} .= ">".$_;
			}
		}
}

if ($type =~ /dat/){
	$/ = "\@",<FNN>;
	while(<FNN>){
		$location = (split /\n/,$_)[0];
		$_ =~  s/\@//;
		#$location =~ s/^Chr/chromosome/;
		#$location =~ s/^p/plasmid/;
		if (not exists $hash{$location}){
			$hash{$location} = "\@".$_;
		}else{
			$hash{$location} .= "\@".$_;
		}
	}
}

if ($type =~ /struc/){
	open(OUT,">$outputdir"."/$sample"."$nanalysis"."$type") or die;
	$/ = "\n\n";
	while(<FNN>){
		my @arr = split /\./,$_,2;
		#$arr[0] =~ s/^Chromosome/Chr/;
	 	#$arr[0] =~ s/^Plasmid/p/;
	 	my $seq = join ".",@arr;
	 	#$seq =~  s/\n\n//;
		#if (not exists $hash{$arr[0]}){
		#	$hash{$arr[0]} = $_."\n\n";
		#}else{
		#	$hash{$arr[0]} .= $_."\n\n";
		#	}
		print OUT "$seq";
	}
}
close FNN;
if (%hash){
	for my $k1 (keys %hash){
		open(OUT,">$outputdir"."/$sample"."_$k1"."$nanalysis"."$type") or die;
		print OUT $hash{$k1};
		close OUT;
	}
}
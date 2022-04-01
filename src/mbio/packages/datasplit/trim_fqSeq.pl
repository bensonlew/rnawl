#!/usr/bin/perl -w
use strict;

my ($input,$output,$start,$length,$min,$max,$end);
my $opt;

while($opt = shift){
	if($opt eq "-i"){
		$input = shift;
	}elsif($opt eq "-o"){
		$output = shift;
	}elsif($opt eq "-s"){
		$start = shift;
	}elsif($opt eq "-l"){
		$length = shift;
	}elsif($opt eq "-m"){
		$min = shift;
	}elsif($opt eq "-x"){
		$max = shift;
	}elsif($opt eq "-e"){
		$end = shift;
	}elsif($opt eq "-h"){
		&usage;exit;
	}else{
		print "Unknown parameter: $opt\n\n";
		&usage;exit;
	}
}

unless($input and $output ){
	&usage;
	exit;
}

$start = $start?$start:1;
#$length = $length?$length:"";
$min = $min?$min:0;
$end = $end?$end:0;
if($length and $min){
	die "min length is longer than valid length: $min >> $length\n" if($min > $length);
}

if($input =~ /gz$/){
	open IN,"gzip -dc $input |" or die "$!\n";
}else{
	open IN,$input or die "$!\n";
}
open OUT,"> $output" or die "$!\n";
while(<IN>){
	chomp;
	my $head = $_;
	chomp(my $seq = <IN>);
	chomp(my $direction = <IN>);
	chomp(my $quality = <IN>);
	my $temp_seq = substr($seq,$start - 1,length($seq) - $end - $start + 1);
	next if(length($temp_seq) < $min);
	next if($max and $max < length($temp_seq));
	my $temp_quality = substr($quality,$start - 1,length($quality) - $end - $start + 1);
	if($length){
		print OUT "$head\n";
		print OUT substr($temp_seq,0,$length),"\n";
		print OUT "$direction\n";
		print OUT substr($temp_quality,0,$length),"\n";
	}else{
		print OUT "$head\n$temp_seq\n$direction\n$temp_quality\n";
	}
}

sub usage{
print <<EOD

	usage: perl $0 -i input.fq -o output -s start.pos -l valid.len -m min.len -x max.len -e trim.end.len
		-i	input fastq file,required
		-o	output fastq file,required
		-s	start position substr seqs,default 1
		-l	length to substr seqs,default length(seq)
		-m	trim seqs shorter than minlength default:0
		-x	trim seqs longer than maxlegth   default:length(seq)
		-e	trim length sequence on the end,default 0
EOD
}

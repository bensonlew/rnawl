#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use File::Basename qw(basename dirname);
use Data::Dumper;
my ($help,$input,$type,$output);

#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
GetOptions(
				"help|?" =>\&USAGE,
				"input:s"=>\$input,
				"type:s"=>\$type,
				"output:s"=>\$output,
				) or &USAGE;
&USAGE unless ($input and $output);
$input = ABSOLUTE_DIR($input);
$output = ABSOLUTE_DIR($output);
$type ||= "rc";
#######################################################################################
open REF,$input;
open OUT,">>$output/output.fasta";
my (%RAWfasta,%NEWfasta,$seq);
while (<REF>) {
	chomp;
	next if ($_ eq "" || /^$/);
	if($_ =~ /\>/){
		$RAWfasta{$_} = "";
		$seq = $_;
	}else{
		$RAWfasta{$seq} .= $_;
	}
}
close REF;

foreach  my $name (keys %RAWfasta) {
	$NEWfasta{$name} = $RAWfasta{$name};
	if($type eq "r"){
		$NEWfasta{$name} = reverse $RAWfasta{$name};
	}elsif($type eq "c"){
		$NEWfasta{$name} = $RAWfasta{$name};
		$NEWfasta{$name} =~ tr/ATCGatcg/TAGCtagc/;
	}elsif($type eq "rc"){
		$NEWfasta{$name} = reverse $RAWfasta{$name};
		$NEWfasta{$name} =~ tr/ATCGatcg/TAGCtagc/;
	}else{
		print "Please provide the right type!";
	}
	print OUT "$name\n$NEWfasta{$name}\n";
}

close OUT;

sub ABSOLUTE_DIR #$pavfile=&ABSOLUTE_DIR($pavfile);
{
	my $cur_dir=`pwd`;chomp($cur_dir);
	my ($in)=@_;
	my $return="";
	if(-f $in){
		my $dir=dirname($in);
		my $file=basename($in);
		chdir $dir;$dir=`pwd`;chomp $dir;
		$return="$dir/$file";
	}elsif(-d $in){
		chdir $in;$return=`pwd`;chomp $return;
	}else{
		warn "Warning just for file and dir \n$in";
		exit;
	}
	chdir $cur_dir;
	return $return;
}


sub USAGE{#
	my $usage =<<"USAGE";
Usage :
	-help	help information
	-input	<file> input fasta file
	-type	<string>	r: reverse, c:complement, rc:reverse and complement
	-output	<file>	output file
USAGE
	print $usage;
	exit;
}
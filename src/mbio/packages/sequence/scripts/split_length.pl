#!/usr/bin/perl -w 
use strict;
unless(@ARGV==3) {
	&usage;
	exit;
}
my($file,$filter,$seq,$len,@new);
open IN,$ARGV[0] or die "$ARGV[0] $!\n";
$filter = $ARGV[1];
open OT1,">$ARGV[2].more$ARGV[1].fa" or die "$ARGV[2].more$ARGV[1] $!\n";
open OT2,">$ARGV[2].less$ARGV[1].fa" or die "$ARGV[2].less$ARGV[1] $!\n";
$/ = ">";
<IN>;
while(<IN>) {
	chomp;
	@new = split /\n/;
	shift @new;
	$seq = join "",@new;
	$len = rindex $seq."\$","\$";
	if($len >= $filter) {
		print OT1 ">",$_;
	}else{
		print OT2 ">",$_; 
	}
}
$/ = "\n";
close IN;
close OT1;
close OT2;

sub usage{
	print STDERR<<USAGE;
	
	perl $0 [fa.file] [cut by length] [out name]
	
USAGE
}

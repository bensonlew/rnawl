#! /usr/bin/perl -w 
use strict;

my $usage = "usage: perl $0 [type] [file num] [old_file.fa] [new_file_prefix] [new_file_rout]
      type
        n\t\t number of file
        s\t\tfile size per file\n";
die $usage if (@ARGV != 5 or ($ARGV[0] ne 'n' and $ARGV[0] ne 's'));

my ($num, $len) = (0, 0);
open IN, $ARGV[2] or die "$ARGV[2]: $!\n";
while (<IN>) {
    chomp;
    if (/^>/) {
	$num ++;
    } else {
	$len += length $_;
    }
}
close IN;
`mkdir -p $ARGV[4]`;
print "num: $num\nlen: $len\n";
my ($file_num, $file_len) = (0, 0);
if ($ARGV[0] eq 'n') {
    $file_num = sprintf ("%1.0f", $num / $ARGV[1]);
    print "file_num: $file_num\n";
    open IN, $ARGV[2] or die "ag: $ARGV[2]: $!\n";
    my $i=1;
    #my $flag = sprintf("%05d",$i);
    my $flag = $i;
    my $count = -1;
    open OUT, ">$ARGV[4]/$ARGV[3]\_$flag.fa" or die "write: $ARGV[3]\_$flag.fa error: $!\n";
    while (<IN>) {
	chomp;
	if (/^>/) {
	    $count ++;
	    if ($count > $file_num and $flag < $ARGV[1] and $count % 2 == 0 ) {
		close OUT;
		$i ++;
		#$flag = sprintf("%05d",$i);
		$flag = $i;
		open OUT, ">$ARGV[4]/$ARGV[3]\_$flag.fa" or die "write: $ARGV[3]\_$flag.fa error: $!\n";
		$count = 0;
	    }
	}
	print OUT "$_\n";
    }
    close IN;
    close OUT;
} else {
    $file_len = sprintf ("%1.0f", $len / $ARGV[1]);
    open IN, $ARGV[2] or die "ag: $ARGV[2]: $!\n";
    my $i=1;
    #my $flag = sprintf("%05d",$i);
    my $flag = $i;
    my $size = 0;
    my $seq = '';
    my $slen = 0;
    open OUT, ">$ARGV[4]/$ARGV[3]\_$flag.fa" or die "write: $ARGV[3]\_$flag.fa error: $!\n";
    while (<IN>) {
	chomp;
	if (/^>/) {
	    if ($seq ne '') {
		print OUT "$seq";
		if ($slen > $file_len and $flag < $ARGV[1]) {
		    close OUT;
		    $i ++;
		    $flag = sprintf("%05d",$i);
		    open OUT, ">$ARGV[4]/$ARGV[3]\_$flag.fa" or die "write: $ARGV[3]\_$flag.fa error: $!\n";
		    $slen = 0;
		}
	    }
	    $seq = "$_\n";
	} else {
	    $seq .= "$_\n";
	    $slen += length $_;
	}
    }
    print OUT "$seq";
    close OUT;
    close IN;
}
print STDERR "finish\t","\n";

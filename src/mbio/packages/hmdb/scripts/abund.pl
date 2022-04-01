#!/usr/bin/perl -w
use strict;

my $script_name = $0;
# my ($arg, $pe, $se, $set, $out, $samp);
my ($arg, $pe, $se, $set, $out);

while ($arg = shift){
	if ($arg eq "-h") { print_usage(); }
	elsif ($arg eq "-pe") { $pe = shift;}
	elsif ($arg eq "-se") {$se = shift;}
	elsif ($arg eq "-set") { $set = shift;}
	elsif ($arg eq "-o") { $out = shift;}
	# elsif ($arg eq "-s") {$samp = shift;}
}

# (($pe or $se) and $set and $out and $samp) || print_usage();
(($pe or $se) and $set and $out) || print_usage();

sub print_usage {
	print <<EOD;
Usage: perl $script_name options
		-pe input pe sam file, required pe or se
		-se input se sam file, required pe or se
		-o output directory, required
		-set geneset, required
		# -s sample name, required
		-h print this help
EOD
exit;
}

my (%reads, @genename, %fq);
# my (%length);
my $start_time = time;

open INS, $set or die "Can not read $set: $!\n";
$/ = ">";
my $i = 0;
<INS>;
while(<INS>) {
	$i++;
	chomp;
	my @list = split /\n/;
	my $name = (split /\s+/, (shift @list))[0];
	# print "$name\n";
	my $seq = join "", @list;
	my $len = length($seq);
	# $lengths{$name} = $len;
	push @genename,$name;
}
close INS;
$/ = "\n";

my $start_read_time = time;
my $cost_time = $start_read_time - $start_time;
print "read geneset over,cost time: $cost_time\n";
if(defined $pe){
print "now start read $pe\n";
open IN, $pe or die "Can not read $pe: $!\n";
my $last_read_name = "name";
while(<IN>){
	chomp;
	my @tmp = split;
	if($last_read_name ne $tmp[0]){
		$reads{$tmp[2]} += 1;
	}
	$last_read_name = $tmp[0];
}
close IN;
}

if(defined $se){
	print "now start read $se\n";
	open IN, $se or die "Can not read $se: $!\n";
	while(<IN>){
		chomp;
		my @tmp = split;
		$reads{$tmp[2]} += 1;
	}
	close IN;
}

foreach my $gen(@genename) {
	$reads{$gen} = 0 unless(defined $reads{$gen})
}

my $start_output_time = time;
$cost_time = $start_output_time - $start_read_time;
print "read pe,se over,cost time: $cost_time\n";
print "now write output\n";
open OUT,">$out" or die "Can not write $out: $!\n";
#print OUT "\t$samp\n";
print OUT "gene\treads_num\n";
foreach my $gen(@genename){
	my $abund = $reads{$gen};
	if($abund != 0){
		print OUT "$gen\t$abund\n";
	}
}
close OUT;
$cost_time = time - $start_output_time;
my $para_time = time - $start_time;
print "output end,cost time: $cost_time\n";
print "All spend time: $para_time\n";

#!/usr/bin/perl -w
use strict;
use warnings;

my $script_name = $0;
my $arg;
my $tRNA;
my $rRNA;
my $location;
while ($arg=shift) {
        if    ($arg eq "-h" ) { print_usage(); }
	elsif ($arg eq "-t" ) { $tRNA = shift; }
	elsif ($arg eq "-r" ) { $rRNA = shift; }
    elsif ($arg eq "-l") { $location = shift; }
}
($tRNA and $tRNA and $location) || print_usage();

my %colors = ("tRNA"=>"#FF0000","5S_rRNA"=>"#008000","16S_rRNA"=>"#0000FF","23S_rRNA"=>"#00BFFF");
#my %colors = ("tRNA"=>"255,0,0","5S_rRNA"=>"0,128,0","16S_rRNA"=>"0,0,255","23S_rRNA"=>"0,191,255");
open SEN,">ncRNA.txt" || die $!;
if ( $tRNA ) {
open T,"<$tRNA" || die "don't open $tRNA!\n";
while(<T>){
	chomp;
	next if(/^Sequence/);
	next if(/^Name/);
	next if(/^---/);
	s/ //g;
	my @line = split /\t/,$_;
    if($line[1]=~/^$location/){
	print SEN "$location $line[9] $line[10] fill_color=(255,0,0)\n" if($line[9]<$line[10]);
	print SEN "$location $line[10] $line[9] fill_color=(255,0,0)\n" if($line[9]>$line[10]);
    }
}
close T;
}
if ( $rRNA ) {
open R,"<$rRNA" || die "don't open $rRNA!\n";
while(<R>){
	chomp;
	next if(/^##gff-version/);
	my @a = split /\t/,$_;
    if($a[1]=~/^$location/){
	if($a[5]=~/\+/){
	    print($a[5]);
		print SEN "$location $a[8] $a[9] fill_color=(0,128,0)\n" if($a[7]=~/5S_rRNA/i);
		print SEN "$location $a[8] $a[9] fill_color=(0,0,255)\n" if($a[7]=~/16S_rRNA/i);
		print SEN "$location $a[8] $a[9] fill_color=(0,191,255)\n" if($a[7]=~/23S_rRNA/i);
	}elsif($a[5]=~/\-/){
		print SEN "$location $a[8] $a[9] fill_color=(0,128,0)\n" if($a[7]=~/5S_rRNA/i);
		print SEN "$location $a[8] $a[9] fill_color=(0,0,255)\n" if($a[7]=~/16S_rRNA/i);
		print SEN "$location $a[8] $a[9] fill_color=(0,191,255)\n" if($a[7]=~/23S_rRNA/i);
	}
    }
}
close R;
}
close SEN;
system("sort -t ' ' -n -k 2 ncRNA.txt -o ncRNA.txt");

sub print_usage {
        print <<EOD;
Usage: perl $script_name options
    -t tRNA.gff
    -r rRNA.gff
    -l location
     -h print this help

EOD
exit;
}

#########################Draw ncRNA legened########################
#use SVG;
#
#my @type= ("tRNA", "5S_rRNA", "16S_rRNA", "23S_rRNA");
#my $legend_num = scalar(@type);
#my $height = 25*$legend_num+10*($legend_num-1)+14;
#my $svg = SVG->new(width=>250, height=>$height);
#for(my $i=0;$i<$legend_num;$i++){
#	$svg->rect(x => 10, y => 7+35*$i, width => 25, height => 25,style=>{stroke=>'black','stroke-width',1,fill=>"$colors{$type[$i]}"});
#	$svg->text(x => 40, y => 7+35*$i+25,'font-size'=>23, 'stroke', 'black', 'stroke-width',0.3, '-cdata', "$type[$i]");
#}
#my $out = $svg->xmlify;
#open SVGFILE,">ncRNA_legend.svg"||die $!;
#print SVGFILE $out;
#close SVGFILE;
#`/usr/bin/convert ncRNA_legend.svg ncRNA_legend.png`;

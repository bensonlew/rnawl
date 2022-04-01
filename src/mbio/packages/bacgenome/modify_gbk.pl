#! /usr/bin/perl -w
use Getopt::Long;
my %opts;
my $VERSION="1.0";
GetOptions( \%opts,"i=s", "des=s","source=s","organism=s","author=s","title=s","journal=s","o=s","h!");
my $usage = <<"USAGE";
       Program : $0
       Version : $VERSION
       Contact : gaohao\@majorbio.com
       Lastest modify:2019-01-08
       Discription:gbk file 
       Usage :perl $0 [options]
			-i*	inputfile	input file;gbk file
			-des		
			-source
			-o	outputdir	output file
			-organism
			-author		
			-title	      
            -journal
			-h			display help message			
        exmaple:perl $0 -i *.gbk -des aaa -source bbb -organism cccc -author dddd -title eeee -journal ffff -o *.gbk
USAGE
die $usage if((! $opts{i}) || $opts{h} );

open (IN,$opts{i}) || die $!;
open (OUT,">$opts{o}") || die $!;
while(<IN>){
  chomp;
  if (/^DEFINITION/){
  print OUT "DEFINITION  $opts{des}\n";
}elsif(/^SOURCE/){
  print OUT "SOURCE      $opts{source}\n";
}elsif(/ORGANISM/){
  print OUT "  ORGANISM  $opts{organism}\n";
}elsif(/AUTHORS/){
  print OUT "  AUTHORS   $opts{author}\n";
}elsif(/TITLE/){
  print	OUT "  TITLE     $opts{title}\n";
}elsif(/JOURNAL/){
  print	OUT "  JOURNAL   $opts{journal}\n";
}else{
  print OUT "$_\n";
}
}

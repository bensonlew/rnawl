#! /usr/bin/perl -w
use warnings;
use strict;

if(@ARGV < 4) {
    print STDERR "fasta-rename.pl  fastaFile prefix\tsort(Y/N)\toutput \n";
    exit;
}
my $file = shift;
my $prefix = shift;
my $sort = shift;
my $output =shift;
open(FILE,"<$file") or die;
open(OUT,">$output ") or die;
open(OUT2,">all.seqid.xls") or die;
print OUT2 "Init_name\tRename\n";
my $num =1;
my $fas;
my %seq;
my %len;
my %index;
$fas = <FILE>;
chomp $fas;
while(<FILE>){
	chomp $_;
	if($_ =~/\>/){
		$index{$_} = $num;
		$num++;
		$fas = $_;
	}else{
		$seq{$fas} .= $_;
		$len{$fas} = length($seq{$fas});
	}
}
close(FILE);
my $index = 0;
if($sort eq "Y"){
	foreach my $fas (sort {$len{$b}<=>$len{$a}} keys %len ){
			$index++;
			print OUT "\>$prefix$index\n",uc($seq{$fas}),"\n";
			my ($id) =$fas =~/^>(.*)$/;
			print OUT2 "$id\t$prefix$index\n";
	}
}else{
	foreach my $fas (sort {$len{$b}<=>$len{$a}} keys %len){
			$index++;
			print OUT "$fas\n",uc($seq{$fas}),"\n";
	}
}

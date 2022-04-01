#! /usr/bin/perl
use strict;
use warnings;

if(@ARGV!=5) {
    print STDERR "extract_seq_byass.pl seq_path seq_list sample_new seq_type outputdir\n";
    exit;
}
my ($file_path,$seq_list,$sample_new,$seq_type,$outputdir)=@ARGV;

my ($name);
my @gene = split /\,/,$seq_list;
if (@gene>1){
	$name = @gene."_".$seq_type;
	}else{
	$name =$seq_list;
	}
open(FNN,$file_path) or die;
open(FNNO,">$outputdir"."/$sample_new"."_$name".".fnn") or die;
$/ = ">",<FNN>;
while(<FNN>){
	my @arr = split /\s/,$_;
		my $gene = $arr[0];
		$_ =~ s/>//;
		if(grep /^$gene$/, @gene){
			print FNNO ">".$_;
		}
}
close FNN;
close FNNO;

#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
my %opts;
my $VERSION="2.0";
GetOptions( \%opts,"t=s","n=i","i=s","headi=s","headt=s","fill=s","h!");


my $usage = <<"USAGE";
        Program : $0
        Version : $VERSION
        Contact: liubinxu
        Lastest modify:2013-04-09

        Discription: select table info of index from
        Usage:perl $0 -i index_list -t file -n index_col
        -i      indexfile
        -t      tablefile
        -n      index line num in tablefile
        -headi	T or F
	-headt	T or F
	-fill	defaut "T",only used when -headi is "T"
        -h      Display this usage information

USAGE

die $usage if ( !( $opts{i} && $opts{t} && $opts{n} ) || $opts{h} );


my %indexchange;

my $index_file = $opts{i};
my $change_file = $opts{t};
my $col_index = $opts{n} - 1;
my $index_head = $opts{headi}?$opts{headi}:"T";
my $table_head = $opts{headt}?$opts{headt}:"T";
my $fill = $opts{fill}?$opts{fill}:"T";
if($index_head ne "T"){
	$fill = "F";
}
#print $fill;
my $element_num;

open (FIN, "<$index_file");
my $index_head_line = "";
if($index_head eq "T"){
        $index_head_line = <FIN>;
	chomp($index_head_line);
	#print $index_head_line;
	my @line = split(/\t/,$index_head_line);
	shift @line;
	$element_num = @line;
	
	$index_head_line = join("\t",@line);

	#print $index_head_line."***\n";
}

#print $index_head_line."***\n";
while(<FIN>){
	chomp;
	my $line = $_;
	if($_=~/^([^\t]*)\t([\s\S]*)$/){
		$indexchange{$1} = $2;
		if($fill eq "T"){
			my @elements = split(/\t/,$2);
			#my $num = @elements;
			while ($element_num != @elements){
				push (@elements,"_");
			}
			my $element_result=join("\t",@elements);
			$indexchange{$1} = $element_result; 
		}
		# print $1."\n";
	}else{
		die "the $index_file format is wrong";
	}		
}
close FIN;

open (FIN2, "<$change_file");
my $table_head_line = "";
if($table_head eq "T"){
        my $table_head_line = <FIN2>; 
	chomp($table_head_line);
	#print $index_head_line."***";
	print $table_head_line."\t".$index_head_line."\n";
}

#if($index_head eq "T" && $table)
#open (OUT1, ">$change_file.new");

while(<FIN2>){
	chomp;
	#print $_;
	my @line=split(/\t/,$_);
	#print $#line;
	if($#line >= $col_index  && $line[$col_index] ne ""){
		if (exists $indexchange{$line[$col_index]}){
			#print $line[$col_index]."\n";
			#print $indexchange{$line[$col_index]}."\n";
			#print $1;
			print $_."\t".$indexchange{$line[$col_index]}."\n";
		}else{
			if($fill eq "F"){
				print "$_"."\n";	
			}else{
				my @elements;
				while ($element_num != @elements){
					push (@elements,"_");
				}
				my $element_result=join("\t",@elements);
				print "$_\t".$element_result."\n";
			}
		}
	}else{
		if (exists $indexchange{$_}){
			#print $1;
			print "$indexchange{$_}"."\n";
		}else{
			if($fill eq "F"){
				print "$_"."\n";	
			}else{
				my @elements;
				while ($element_num != @elements){
					push (@elements,"_");
				}
				my $element_result=join("\t",@elements);
				print "$_\t".$element_result."\n";
			}
			
		}
	}
}	
close FIN2;
#close OUT1;

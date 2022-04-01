#!/usr/bin/perl -w
use strict;
use warnings;
  
die "perl $0 XX.scaf2contig  XX.agp  Final.fa cutoff(200bp)\n" if(@ARGV!=4);
my ($line,@inf,%contig,%scaffold,%scaffold_id);

open IN, "$ARGV[0]" or die "can not open file: $ARGV[0]\n";
while($line=<IN>){
	chomp $line;
	$line=~s/^>//;
	my $id=$line;
	$line=<IN>;
	chomp $line;
	$contig{$id}=$line;
}
close IN;

open IN, "$ARGV[1]" or die "can not open file: $ARGV[1]\n";
my $i=0;
my $pre_id="";
while($line=<IN>){
	chomp $line;@inf=split /\t/,$line;
	if($pre_id ne $inf[0]){
		$i++;
		$scaffold_id{$i}=$inf[0];
		$pre_id=$inf[0];
	}
	if($inf[4] eq "W" && $inf[7] >=200){
		$scaffold{$inf[0]}.=$line."\n";
	}elsif($inf[4] eq "N"){
		$scaffold{$inf[0]}.=$line."\n";
	}else{
		$line=join("\t",@inf[0..3]);
		$line.="\tN\t".$inf[6]."\tfragment\tyes";
		$scaffold{$inf[0]}.=$line."\n";
	}
}
close IN;

open OA, ">$ARGV[2]" or die "can not open file: $ARGV[2]\n";
foreach my $i (sort {$a<=>$b} keys %scaffold_id){
	@inf=split /\n/,$scaffold{$scaffold_id{$i}};
	my $seq="";
	for(my $j=0;$j<=$#inf;$j++){
		my @ele=split /\t/,$inf[$j];
		if($ele[4] eq "W"){
			last;
		}
		if($ele[4] eq "N"){
			$inf[$j]="";
		}
	}
	for(my $j=$#inf;$j>=0;$j--){
		next if($inf[$j] eq "");
		my @ele=split /\t/,$inf[$j];
		if($ele[4] eq "W"){
			last;
		}
		if($ele[4] eq "N"){
			$inf[$j]="";		
		}
	}
	for(my $j=0;$j<=$#inf;$j++){
		next if($inf[$j] eq "");	
		my @ele=split /\t/,$inf[$j];
		if($ele[4] eq "W"){
			$seq.=$contig{$ele[5]};
		}elsif($ele[4] eq "N"){
			$seq.="N"x$ele[5];		
		}
	}
	if($seq ne ""){
		print OA ">$scaffold_id{$i}\n$seq\n";
	}
}
close OA;




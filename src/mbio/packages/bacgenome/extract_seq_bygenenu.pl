#! /usr/bin/perl
use strict;
use warnings;

if(@ARGV!=6) {
    print STDERR "extract_seq_bygenenu.pl seq_prefix_path type gene_list sample gene_new outputdir\n";
    exit;
}
my ($file_path,$type,$gene_list,$sample,$gene_new,$outputdir)=@ARGV;
my %hash;
if ($gene_new){
	my @arr = split /,/,$gene_new;
	foreach(@arr){
		my @a = split /:/,$_;
		$hash{$a[0]} = $a[1];
	}
}
my ($name,@nu,$loc);
my @gene = split /\,/,$gene_list;
if (@gene>1){
	$name = @gene."_gene";
	}else{
	$name =$gene_list;
	}
open(FNN,"<$file_path.fnn") or die;
open(FNNO,">$outputdir"."/$sample"."_$name".".fnn") or die;
$/ = ">",<FNN>;
while(<FNN>){
	#$nu++;
	my @arr = split /\s/,$_;
	my $gene = $arr[0];
    $loc = $arr[3];
    $loc =~ s/_tRNA[0]*/\.trna/;
	$_ =~ s/>//;
	if ($gene =~ /^(\S+)[0-9]{4,}/){
            my $prefix = $1;
		if (exists($hash{$1})){
			$_ =~ s/^$prefix/$hash{$prefix}/;
		}
	}
	if(grep /^$gene$/, @gene){
		print FNNO ">".$_;
		push(@nu,$loc);
        print $loc;
        #print $nu."\n";
	}
}
close FNN;
close FNNO;

if ($type =~ /^gene$/){
	open(FAA,"<$file_path.faa") or die;
	open(FAAO,">$outputdir"."/$sample"."_$name".".faa") or die;
	$/ = ">",<FAA>;
	while(<FAA>){
		my @arr = split /\s/,$_;
		my $gene = $arr[0];
		$_ =~ s/>//;
		if ($gene =~ /^(\S+)[0-9]{4,}/){
            my $prefix = $1;
			if (exists($hash{$1})){
				$_ =~ s/^$prefix/$hash{$prefix}/;
				}
			}
		if(grep /^$gene$/, @gene){
			print FAAO ">".$_;
		}
	}
	close FAA;
	close FAAO;
}

if ($type =~ /^trna$/){
	open(STRUC,"<$file_path.struc") or die;
	open(STRUCO,">$outputdir"."/$sample"."_$name".".struc") or die;
	$/ = "\n\n";
	while(<STRUC>){
		my @arr = split /\s/,$_;
		my $gene = $arr[0];
		$_ =~ tr/>//;
		if(grep /^$gene$/, @nu){
			print STRUCO $_;
		}
	}
	close STRUC;
	close STRUCO;
}




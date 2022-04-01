#!/usr/bin/perl
use strict;
use warnings;
 
if(@ARGV!=3){
	print STDERR "perl orthomcl_stat.pl mclFile faa.list outDir\n";
	exit;
}

my $mclFile = shift;
my $faalist= shift;
my $outDir= shift;

if (!-d $outDir){
    system ("mkdir -p $outDir");
}

(-s $faalist) || die "$faalist is not exists!\n";
open LI,"<$faalist" || die $!;
my (@spec_names,%choose);
while(<LI>){
	chomp;
	my @l=split /\t/;
	push(@spec_names,$l[0]);
	(-s $l[1]) || die "$l[1] is not exists!\n";
	open FA,"<$l[1]" || die $!;
	my $id;
	while(<FA>){
		chomp;
		if(/^>/){
			$id = (split /\s/,$_)[0];
			$id =~ s/>//g;
			$choose{$l[0]}{$id}=0;
		}
	}
	close FA;
}
close LI;

my $sample_len = @spec_names;
print "$sample_len\n";

my (%family_number,%gene_number);
open OUT1,"> $outDir/orthomcl_stat.xls" || die $!;
print OUT1 "#ClusterID\t",join("\t",@spec_names),"\tSample_exists\n";
open(ORTHOMCL,"<$mclFile") or die;
while(<ORTHOMCL>){
	my %species_orf;
	for(my $j=0;$j<@spec_names;$j++){
		$species_orf{$spec_names[$j]} = "";
	}
	chomp $_;
	my @line = split /\t/, $_;
	my @clsuter = split(/\s/,$line[0]);
	my $newcluster = $clsuter[0].")";
	print OUT1 "$newcluster";
	my @line2 = split /\s+/,$line[1];
	for(my $i=0;$i<@line2;$i++){	
			if($line2[$i] =~ /(\S+)\((\S+)\)/){
				my @gene_name = split(/\|/,$1) ;
				#print $gene_name[1] ;
				#$choose{$2}{$1} = 1;
				$choose{$2}{$gene_name[1]} = 1;
				$species_orf{$2} .= $1.";";	
			}	
	}
	my (@sample_exists, @tmp_gene_num,$tmp_len);
	$tmp_len =0;
	for (my $m=0;$m<@spec_names;$m++){
		if ($species_orf{$spec_names[$m]} ne ""){
			my $sample_genes = $species_orf{$spec_names[$m]};
			$sample_genes =~ s/;$//;
			print OUT1 "\t$sample_genes";
			@tmp_gene_num = split(/;/,$sample_genes);
			$tmp_len += @tmp_gene_num;
			push (@sample_exists,$spec_names[$m]);
		}else{
			print OUT1 "\t-";
		}
	}
	my $sample_exists = join(";",@sample_exists);
	print OUT1 "\t$sample_exists\n";
	if (!exists $family_number{$sample_exists}){
		$gene_number{$sample_exists} = $tmp_len;
		$family_number{$sample_exists} = 1;
	}else{
		$gene_number{$sample_exists} += $tmp_len;
		$family_number{$sample_exists} += 1;
	}
}


for(my $i=0;$i<@spec_names;$i++){
	my $single_num = 1;
	for my $j(keys %{$choose{$spec_names[$i]}}){
		if($choose{$spec_names[$i]}{$j} == 0){
			#print OUT1 "SINGLE"."-"x($i+1)."\t$j\n";
			print OUT1 "SINGLE".$single_num."\t-"x$i."\t"."$spec_names[$i]|".$j."\t-"x($sample_len - $i-1)."\t$spec_names[$i]\n";
			$single_num += 1;
			if (!exists $family_number{$spec_names[$i]}){
				$family_number{$spec_names[$i]} = 1;
				$gene_number{$spec_names[$i]} = 1;
			}
			else{
				$family_number{$spec_names[$i]} += 1;
				$gene_number{$spec_names[$i]} += 1;
			}
		}
	}
}
close(ORTHOMCL);
close OUT1;

##venn##

open OUT2,">$outDir/orthogenes_venn.xls" || die $!;
print OUT2 "#Lables\tFamily_number\tGene_number\n";

my (@all_condition,$condition);
for(my $i=0;$i<@spec_names;$i++){
	my $result = choose($i+1,@spec_names);
	if($i eq 0){
		for (@$result){
			$condition = join("",@$_)." only";
			push (@all_condition,$condition);
			}
	}else{
		for (@$result){
			$condition = join(" & ",@$_);
			push (@all_condition,$condition);
		}
	}
}

for my $each_conditin(@all_condition){
	my $origin_condition = $each_conditin;
	#print $each_conditin,"\n";
	if ($each_conditin =~ /&/){
		$each_conditin =~ s/ & /;/g;
	}elsif($each_conditin =~ /only/){
		$each_conditin =~ s/ only//g;
	}
	my ($condtion_fa_num,$condtion_gene_num);
	if (exists $family_number{$each_conditin}){
		$condtion_fa_num = $family_number{$each_conditin};
		$condtion_gene_num = $gene_number{$each_conditin};
	}else{
		$condtion_fa_num = 0;
		$condtion_gene_num = 0;
	}
	print OUT2 "$origin_condition\t$condtion_fa_num\t$condtion_gene_num\n";
}

close OUT2;

sub choose {
    my($n, @data) = @_;    # 需要从 @data 中取出 $n 项
    my @result;
    #print $n,">>>",join(";",@data),"\n";
    return [map {[$_]} @data] if $n == 1;  # 只取一个时用
    while (1) {
        last if @data < $n; # 退出条件
        my $item = shift @data;
        my $ret = choose($n-1, @data);
        for (@$ret) {
            unshift @$_, $item;
            push @result, $_;
            }
        }
    return \@result;
}



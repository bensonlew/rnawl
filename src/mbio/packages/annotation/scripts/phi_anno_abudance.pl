#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
my %opts;
GetOptions(\%opts,"q=s","p=s","o=s");

my $usage = <<"USAGE";
        Program : $0
        Contact : haidong.gu\@majorbio.com
        Discription: Calculate phi abundance.
        Usage:perl $0 [options]
                -q     input gene_phi_anno detail file
                -p     input gene profile file
                -o     output Dir

USAGE

die $usage  if (!($opts{q}&&$opts{o}&&$opts{p}));

`mkdir -p $opts{o}` unless -e ($opts{o});

my (%gene_profile,@sams,%profile_genes);
open INP,$opts{p} or die "can not open $opts{p}!\n";
while(<INP>){
    chomp;
    if (/^GeneID/){
    @sams = split /\t/;
    #### 筛除head的GeneID列和最后的total列
    shift @sams;
    #pop @sams;
    }else{
    my @tmp = split /\t/;
    my $gene = $tmp[0];
    $profile_genes{$gene} = 1;
    shift @tmp;
    my $tmp_profile = join("\t",@tmp);
    $gene_profile{$gene} = $tmp_profile;
    }
}
close INP;

my (%host, %host_des, %pathogen, %phenotype, %protein);
open INF, $opts{q} or die 'can not open $opts{q}!\n';
# open OUT, ">$opts{o}/gene_phi_anno.xls" or die "can not open $opts{o}/gene_phi_anno.xls\n";
while(<INF>){
	chomp;
	next if(/^#/);
	my @tmp = split /\t/;
	my $query = $tmp[0];
	# $query =~ s/_1//;
	# print OUT "$query\t", join("\t", @tmp[1:]), "\n";
	if(exists $profile_genes{$query}){
	    my $protein = $tmp[2];
		my $host = $tmp[9];
		my $host_description = $tmp[7];
		my $pathogen = $tmp[5];
		my $phenotype = $tmp[6];
		$protein{$protein} .= ",".$query;
		my @host_list = split /;/, $host;
		my $tmp_length = @host_list;
		if ($tmp_length > 1){
			my @host_des_list = split /;/, $host_description;
			for(my $i = 0;$i < @host_list; $i ++){
				next if (exists $host{$host_list[$i]} and $host{$host_list[$i]} =~ /$query$/);
				$host{$host_list[$i]} .= ",".$query;
				if (@host_des_list > 1){
					$host_des{$host_list[$i]} = $host_des_list[$i];
				}else{
					$host_des{$host_list[$i]} = $host_description;
				}
			}
		}else{
			$host{$host} .= ",".$query;
			$host_des{$host} = $host_description;
		}
		my @phenotype_list = split /;/, $phenotype;
		$tmp_length = @phenotype_list;
		if ($tmp_length > 1){
			for (my $i = 0; $i < @phenotype_list; $i ++){
				next if (exists $phenotype{$phenotype_list[$i]} and $phenotype{$phenotype_list[$i]} =~ /$query$/);
				$phenotype{$phenotype_list[$i]} .= ",".$query;
			}
		}else{
			$phenotype{$phenotype} .= ",".$query;
		}
		my @pathogen_list = split /;/, $pathogen;
		$tmp_length = @pathogen_list;
		if ($tmp_length > 1){
			for (my $i = 0;$i < @pathogen_list;$i ++){
				next if (exists $pathogen{$pathogen_list[$i]} and $pathogen{$pathogen_list[$i]} =~ /$query$/);
				$pathogen{$pathogen_list[$i]} .= ",".$query;
			}
		}else{
			$pathogen{$pathogen} .= ",".$query;
		}
		#$host{$host} .= ",".$query;
		#$host_des{$host} = $host_description;
		#$phenotype{$phenotype} .= ",".$query;
		# $pathogen{$pathogen} .= ",".$query;
	}
}
# close OUT;
close INF;

open OUTPR, ">$opts{o}/phi_protein_profile.xls" or die "can not open $opts{o}/phi_protein_profile.xls!\n";
open OUTH, ">$opts{o}/phi_host_profile.xls" or die "can not open $opts{o}/phi_host_profile.xls!\n";
open OUTPA, ">$opts{o}/phi_pathogen_profile.xls" or die "can not open $opts{o}/phi_pathogen_profile.xls!\n";
open OUTPH, ">$opts{o}/phi_phenotype_profile.xls" or die "can not open $opts{o}/phi_phenotype_profile.xls!\n";
my $samples = join("\t", @sams);
print OUTPR "#Protein\t$samples\n";
#print OUTH "#Host\t$samples\tDescription\n";
print OUTH "#Host\t$samples\tHost_classification\n";
print OUTPA "#Pathogen\t$samples\n";
print OUTPH "#Phenotype\t$samples\n";

foreach my $eachpr (sort keys %protein){
	print OUTPR "$eachpr";
	my $pr_gene = $protein{$eachpr};
	my @genes = split(',', $pr_gene);
	shift @genes;
	for (my $i = 0; $i < @sams; $i++){
		my $pr_abu;
		foreach my $gene (@genes){
			my @sams_profile = split(/\t/, $gene_profile{$gene});
			my $sam_profile = $sams_profile[$i];
			$pr_abu += $sam_profile;
		}
		print OUTPR "\t$pr_abu";
	}
	print OUTPR "\n";
}
close OUTPR;

foreach my $eachhost (sort keys %host){
    print OUTH "$eachhost";
	my $host_gene = $host{$eachhost};
	my @genes = split(',', $host_gene);
	shift @genes;
	for (my $i = 0; $i < @sams; $i++){
		my $host_abu;
		foreach my $gene (@genes){
			my @sams_profile = split(/\t/, $gene_profile{$gene});
			my $sam_profile = $sams_profile[$i];
			$host_abu += $sam_profile;
		}
		print OUTH "\t", $host_abu;
	}
	print OUTH "\t$host_des{$eachhost}\n";
}
close OUTH;

foreach my $eachpathogen (sort keys %pathogen){
	print OUTPA "$eachpathogen";
	my $pathogen_gene = $pathogen{$eachpathogen};
	my @genes = split(/,/, $pathogen_gene);
	shift @genes;
	for (my $i = 0; $i < @sams; $i++){
		my $pathogen_abu;
		foreach my $gene (@genes){
			my @sams_profile = split(/\t/, $gene_profile{$gene});
			my $sam_profile = $sams_profile[$i];
			$pathogen_abu += $sam_profile;
		}
		print OUTPA "\t", $pathogen_abu;
	}
	print OUTPA "\n";
}
close OUTPA;

foreach my $eachphenotype (sort keys %phenotype) {
	print OUTPH "$eachphenotype";
	my $phenotype_gene = $phenotype{$eachphenotype};
	my @genes = split(/,/, $phenotype_gene);
	shift @genes;
	for (my $i = 0; $i < @sams; $i++){
		my $phenotype_abu;
		foreach my $gene (@genes){
			my @sams_profile = split(/\t/, $gene_profile{$gene});
			my $sam_profile = $sams_profile[$i];
			$phenotype_abu += $sam_profile;
		}
		print OUTPH "\t", $phenotype_abu;
	}
	print OUTPH "\n";
}
close OUTPH;

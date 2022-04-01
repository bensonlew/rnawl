#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
my %opts;
GetOptions(\%opts,"i=s","p=s","o=s");

my $usage = <<"USAGE";
        Program : $0
        Contact : shaohua.yuan\@majorbio.com
        Discription: Calculate probio abundance.
        Usage:perl $0 [options]
                -i     input gene_probio_anno detail file
                -p     input gene profile file
                -o     output dir

USAGE

die $usage  if (!($opts{i}&&$opts{o}&&$opts{p}));

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

my (%probio_detail,%probio_genus,%probio_genes);
open INF,$opts{i} or die "can not open $opts{i}!\n";
while(<INF>){
    chomp;
    next if(/^#/);
    my @tmp = split /\t/;
    my $query = $tmp[0];
    if(exists $gene_profile{$query}){
        my $probio = $tmp[1];
        print "$probio\n";
        my $strain = $tmp[3];
        my $genus = $tmp[4];
        $genus =~ s/g__//;
        my $name = $probio.";".$strain;
        shift @tmp;
        $probio_detail{$name} = join("\t",@tmp);
        $probio_genus{$name} = $genus;
        $probio_genes{$name} .= ",".$query;
    }
}
close INF;


open OUTC, "> $opts{o}/probio_abun.xls" or die "can not open $opts{o}/probio_abun.xls!\n";
open OUTA, "> $opts{o}/probio_anno.xls" or die "can not open $opts{o}/probio_anno.xls!\n";
print OUTA "#Probiotic_name\tBrand\tStrain\tGenus\tCommercial Development Stage\tUse in\tProbiotic Effect\tDisease_class\tICD 10 Disease Code\tLineage\tProbio_ID\n";
my $samples = join("\t",@sams);
print OUTC "#Probiotics\tProbio genus\t","$samples\n";
foreach my $eachpro (sort keys %probio_detail){
    my @species_strain = split(/;/, $eachpro);
    my $species = $species_strain[0];
    my $strain = $species_strain[1];
    my $probio_name;
    if($strain =~ /-/){
        $probio_name = $species;
    }else{
        $probio_name = $eachpro;
        $probio_name =~ s/;/_/g;
    }
    print OUTA "$probio_detail{$eachpro}\n";
    my $genus = $probio_genus{$eachpro};
    print OUTC "$probio_name\t$genus";
    my $probio_gene = $probio_genes{$eachpro};
    my @genes = split(/,/,$probio_gene);
    shift @genes;
    $probio_gene =~ s/,//;
    my $gene_count = @genes;
    for(my $i = 0;$i < @sams ;$i++){
        my $pro_abu;
        foreach my $gene (@genes){
            my @sams_profile = split(/\t/,$gene_profile{$gene});
            my $sam_profile = $sams_profile[$i];
            $pro_abu += $sam_profile;
       }
    print OUTC "\t",$pro_abu;
    }
    print OUTC "\n";
}
close OUTC;
close OUTA;



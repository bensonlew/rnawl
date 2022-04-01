#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
my %opts;
GetOptions(\%opts,"q=s","p=s","o=s");

my $usage = <<"USAGE";
        Program : $0
        Contact : shaohua.yuan\@majorbio.com
        Discription: Calculate OPU unique protein abundance.
        Usage:perl $0 [options]
                -q     input gene_ardb_anno detail file
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

my (%type,%class,%ARG,%class_des,%antitype,%anticlass);
open INF,$opts{q} or die "can not open $opts{q}!\n";
while(<INF>){
    chomp;
    next if(/^#/);
    my @tmp = split /\t/;
    my $query = $tmp[0];
    if(exists $profile_genes{$query}){
        my $type = $tmp[2];
        my $class = $tmp[7];
        my $class_des = $tmp[8];
        my $antitype = $tmp[5];
        $class_des{$class} = $class_des;
        $antitype{$type} = $antitype;
        my $ARG = $tmp[1];
        $ARG{$ARG} .= ",".$query;
        $type{$type} .= ",".$query;
        $class{$class} .=",".$query;
        my $anti_class = $tmp[9];
        $anticlass{$anti_class} .= ",".$query;
    }
}
close INF;


open OUTC, "> $opts{o}/ardb_class_profile.xls" or die "can not open $opts{o}/ardb_class_profile.xls!\n";
open OUTT, "> $opts{o}/ardb_type_profile.xls" or die "can not open $opts{o}/ardb_type_profile.xls!\n";
open OUTG, "> $opts{o}/ardb_ARG_profile.xls" or die "can not open $opts{o}/ardb_ARG_profile.xls!\n";
open OUTS,"> $opts{o}/gene_ardb_class_stat.xls" or die "can not open $opts{o}/gene_ardb_class_stat.xls!\n";
open OUTM,"> $opts{o}/ardb_Antibiotic_class_profile.xls" or die "can not open $opts{o}/ardb_Antibiotic_class_profile.xls!\n";

my $samples = join("\t",@sams);
print OUTC "#Class\t","$samples\tDescription\n";
print OUTT "#Type\t","$samples\tAntibiotic_type\n";
print OUTG "#ARG\t","$samples\n";
print OUTS "#Class\tCount\tGene_list\n";
print OUTM "Antibiotic_class\t","$samples\n";

foreach my $eachclass (sort keys %class){
    print OUTC "$eachclass\t";
    print OUTS "$eachclass\t";
    my $class_gene = $class{$eachclass};
    my @genes = split(/,/,$class_gene);
    shift @genes;
    $class_gene =~ s/,//;
    my $gene_count = @genes;
    print OUTS "$gene_count\t$class_gene\n";
    for(my $i = 0;$i < @sams ;$i++){
        my $class_abu;
        foreach my $gene (@genes){
            my @sams_profile = split(/\t/,$gene_profile{$gene});
            my $sam_profile = $sams_profile[$i];
            $class_abu += $sam_profile;
       }
    print OUTC $class_abu,"\t";
    }
    print OUTC "$class_des{$eachclass}\n";
}
close OUTC;
close OUTS;

foreach my $eachtype (sort keys %type){
    print OUTT "$eachtype\t";
    my $type_gene = $type{$eachtype};
    my @genes = split(/,/,$type_gene);
    shift @genes;
    for(my $i = 0;$i < @sams ;$i++){
        my $type_abu;
        foreach my $gene (@genes){
            my @sams_profile = split(/\t/,$gene_profile{$gene});
            my $sam_profile = $sams_profile[$i];
            $type_abu += $sam_profile;
       }
    print OUTT $type_abu,"\t";
    }
    print OUTT "$antitype{$eachtype}\n";
}
close OUTT;

foreach my $eachARG (sort keys %ARG){
    print OUTG "$eachARG";
    my $ARG_gene = $ARG{$eachARG};
    my @genes = split(/,/,$ARG_gene);
    shift @genes;
    for(my $i = 0;$i < @sams ;$i++){
        my $ARG_abu;
        foreach my $gene (@genes){
            my @sams_profile = split(/\t/,$gene_profile{$gene});
            my $sam_profile = $sams_profile[$i];
            $ARG_abu += $sam_profile;
       }
    print OUTG "\t",$ARG_abu;
    }
    print OUTG "\n";
}
close OUTG;

foreach my $eachanticlass (sort keys %anticlass){
    print OUTM "$eachanticlass";
    my $class_gene = $anticlass{$eachanticlass};
    my @genes = split(/,/,$class_gene);
    shift @genes;
    for(my $i = 0;$i < @sams ;$i++){
        my $class_abu;
        foreach my $gene (@genes){
            my @sams_profile = split(/\t/,$gene_profile{$gene});
            my $sam_profile = $sams_profile[$i];
            $class_abu += $sam_profile;
       }
    print OUTM "\t",$class_abu;
    }
    print OUTM "\n";
}
close OUTM;

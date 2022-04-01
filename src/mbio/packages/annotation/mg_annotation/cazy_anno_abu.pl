#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
my %opts;
GetOptions(\%opts,"q=s","p=s","des=s","o=s");

my $usage = <<"USAGE";
        Program : $0
        Contact : shaohua.yuan\@majorbio.com
        Discription: Calculate cazy abundance.
        Usage:perl $0 [options]
                -q     input gene_ardb_anno detail file
                -p     input gene profile file
                -o     output Dir
                -des   FamInfo.txt path

USAGE

die $usage  if (!($opts{q}&&$opts{o}&&$opts{p}&&$opts{des}));

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

my (%fam_des);
open IND,$opts{des} or die "can not open $opts{des}!\n";
while(<IND>){
    chomp;
    my @tmp = split /\t/;
    my $fa = $tmp[0];
    my $fa_des = $tmp[4];
    $fam_des{$fa} = $fa_des;
}
close IND;

my (%family,%class,%class_des);
open INF,$opts{q} or die "can not open $opts{q}!\n";
while(<INF>){
    chomp;
    next if(/^#/);
    my @tmp = split /\t/;
    my $query = $tmp[0];
    if(exists $profile_genes{$query}){
        my $family = $tmp[1];
        # 暂时测试用
        my $class = $tmp[2];
        my $class_des = $tmp[3];
        $class_des{$class} = $class_des;
        $family{$family} .= ",".$query;
        $class{$class} .=",".$query;
    }
}
close INF;

open OUTC, "> $opts{o}/cazy_class_profile.xls" or die "can not open $opts{o}/cazy_class_profile.xls!\n";
open OUTT, "> $opts{o}/cazy_family_profile.xls" or die "can not open $opts{o}/cazy_family_profile.xls!\n";
open OUTS, "> $opts{o}/gene_cazy_class_stat.xls" or die "can not open $opts{o}/gene_cazy_class_stat.xls!\n";
open OUTP, "> $opts{o}/gene_cazy_family_stat.xls" or die "can not open $opts{o}/gene_cazy_family_stat.xls!\n";

my $samples = join("\t",@sams);
print OUTC "#Class\t","$samples\tDescription\n";
print OUTT "#Family\t","$samples\tDescription\n";
print OUTS "Class\tGene_counts\tGene_list\tDescription\n";
print OUTP "Family\tGene_counts\tGene_list\tDescription\n";
foreach my $eachclass (sort keys %class){
    print OUTC "$eachclass\t";
    my $class_gene = $class{$eachclass};
    my @genes = split(/,/,$class_gene);
    shift @genes;
    my $cl_counts = @genes;
    $class_gene =~ s/,//;
    print OUTS "$eachclass\t$cl_counts\t$class_gene\t$class_des{$eachclass}\n";
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

foreach my $eachfamily (sort keys %family){
    print OUTT "$eachfamily\t";
    my $family_gene = $family{$eachfamily};
    my @genes = split(/,/,$family_gene);
    shift @genes;
    my $fa_counts = @genes;
    $family_gene =~ s/,//;
    print OUTP "$eachfamily\t$fa_counts\t$family_gene\t$fam_des{$eachfamily}\n";
    for(my $i = 0;$i < @sams ;$i++){
        my $family_abu;
        foreach my $gene (@genes){
            my @sams_profile = split(/\t/,$gene_profile{$gene});
            my $sam_profile = $sams_profile[$i];
            $family_abu += $sam_profile;
       }
    print OUTT $family_abu,"\t";
    }
    print OUTT "$fam_des{$eachfamily}\n";
}
close OUTT;

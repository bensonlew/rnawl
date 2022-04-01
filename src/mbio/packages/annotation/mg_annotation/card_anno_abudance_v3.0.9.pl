#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
my %opts;
GetOptions(\%opts,"q=s","p=s","o=s");

my $usage = <<"USAGE";
        Program : $0
        Contact : shaohua.yuan\@majorbio.com
        Discription: Calculate card abundance.
        Usage:perl $0 [options]
                -q     input gene_card_anno detail file
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

my (%class,%ARO,%ARO_name,%ARO_des, %class_des, %amr_family, %drug_class, %type,%type_aro, %resistance, %resistance_aro );
open INF,$opts{q} or die "can not open $opts{q}!\n";
while(<INF>){
    chomp;
    next if(/^#/);
    my @tmp = split /\t/;
    my $query = $tmp[0];
    if(exists $profile_genes{$query}){
        my $class = $tmp[5];
#        my $class_des = $tmp[6];
        my $ARO = $tmp[1];
        my $ARO_name = $tmp[2];
        my $des = $tmp[3];
        my $amr = $tmp[4];
        my $drug = $tmp[5];
        my $type_des = $tmp[6];
        my $resis = $tmp[7];
        $ARO_name{$ARO} =  $ARO_name;
        $ARO_des{$ARO} = $des;
        $amr_family{$ARO} = $amr;
        $drug_class{$ARO} = $drug;
        $type{$type_des} .= ",".$query;
        $type_aro{$ARO} = $type_des;
        $resistance_aro{$ARO} = $resis;
        $ARO{$ARO} .= ",".$query;
        my @classes = split /;/,$class;
        my $tmp_length = @classes;
        if ($tmp_length > 1){
            for(my $i = 0;$i < @classes;$i++){
            $class{$classes[$i]} .=",".$query;
            }
        }else{
            $class{$class} .=",".$query;
        }
        my @resis_classes = split /;/,$resis;
        my $tmp_length2 = @resis_classes;
        if ($tmp_length2 > 1){
            for(my $i = 0;$i < @resis_classes;$i++){
            $resistance{$resis_classes[$i]} .=",".$query;
            }
        }else{
            $resistance{$resis} .=",".$query;
        }
    }
}
close INF;

open OUTC, "> $opts{o}/card_Drug_class_profile.xls" or die "can not open $opts{o}/card_class_profile.xls!\n";
open OUTG, "> $opts{o}/card_ARO_profile.xls" or die "can not open $opts{o}/card_ARO_profile.xls!\n";
open OUTN, "> $opts{o}/card_ARO_gene_number.xls" or die "can not open $opts{o}/card_ARO_gene_number.xls!\n";
open OUTA, "> $opts{o}/card_ARO_profile_all.xls" or die "can not open $opts{o}/card_ARO_profile_all.xls!\n";
open OUTR, "> $opts{o}/card_Resistance_Mechanism_profile.xls" or die "can not open $opts{o}/card_ARO_profile.xls!\n";
open OUTT, "> $opts{o}/card_Antibiotic_class_profile.xls" or die "can not open $opts{o}/card_Antibiotic_type_profile.xls!\n";
my $samples = join("\t",@sams);
print OUTC "#Drug_class\t","$samples\n";
print OUTG "#ARO\tAro_name\t$samples\tDescription\n";
print OUTN "#ARO\tGene_number\tGenes\n";
print OUTA "#ARO\tAro_name\tDescription\tAMR_Gene_Family\tDrug_Class\tAntibiotic_class\tResistance_Mechanism\t$samples\n";
print OUTR "Resistance_Mechanism\t$samples\n";
print OUTT "Antibiotic_class\t$samples\n";

foreach my $eachclass (sort keys %class){
    print OUTC "$eachclass\t";
    my $class_gene = $class{$eachclass};
    my @genes = split(/,/,$class_gene);
    shift @genes;
    for(my $i = 0;$i < @sams ;$i++){
        my $class_abu;
        foreach my $gene (@genes){
            my @sams_profile = split(/\t/,$gene_profile{$gene});
            my $sam_profile = $sams_profile[$i];
            $class_abu += $sam_profile;
       }
    print OUTC $class_abu,"\t";
    }
#    print OUTC "$class_des{$eachclass}\n";
    print OUTC "\n";
}
close OUTC;

foreach my $eachARO (sort keys %ARO){
    print OUTG "$eachARO\t$ARO_name{$eachARO}\t";
    print OUTA "$eachARO\t$ARO_name{$eachARO}\t$ARO_des{$eachARO}\t$amr_family{$eachARO}\t$drug_class{$eachARO}\t$type_aro{$eachARO}\t$resistance_aro{$eachARO}\t";
    print OUTN "$eachARO\t";
    my $ARO_gene = $ARO{$eachARO};
    my @genes = split(/,/,$ARO_gene);
    shift @genes;
    my $gene_names = join(",",@genes);
    my $number = @genes;
    print OUTN "$number\t$gene_names\n";

    for(my $i = 0;$i < @sams ;$i++){
        my $ARO_abu;
        foreach my $gene (@genes){
            my @sams_profile = split(/\t/,$gene_profile{$gene});
            my $sam_profile = $sams_profile[$i];
            $ARO_abu += $sam_profile;
       }
    print OUTG $ARO_abu,"\t";
    print OUTA $ARO_abu,"\t";
    }
    print OUTG "$ARO_des{$eachARO}\n";
    print OUTA "\n";
}
close OUTG;

foreach my $eachtype (sort keys %type){
    print OUTT "$eachtype\t";
    my $type_gene = $type{$eachtype};
    my @genes = split(/,/,$type_gene);
    shift @genes;
    for(my $j = 0;$j < @sams ;$j++){
        my $type_abu;
        foreach my $gene (@genes){
            my @sams_profile = split(/\t/,$gene_profile{$gene});
            my $sam_profile = $sams_profile[$j];
            $type_abu += $sam_profile;
       }
    print OUTT $type_abu,"\t";
    }
    print OUTT "\n";
}
close OUTT;

foreach my $eachresistance (sort keys %resistance){
    print OUTR "$eachresistance\t";
    my $resistance_gene = $resistance{$eachresistance};
    my @genes = split(/,/,$resistance_gene);
    shift @genes;
    for(my $m = 0;$m < @sams ;$m++){
        my $resistance_abu;
        foreach my $gene (@genes){
            my @sams_profile = split(/\t/,$gene_profile{$gene});
            my $sam_profile = $sams_profile[$m];
            $resistance_abu += $sam_profile;
       }
    print OUTR $resistance_abu,"\t";
    }
    print OUTR "\n";
}
close OUTR;
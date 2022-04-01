#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
my %opts;
GetOptions(\%opts,"q=s","p=s","o=s","e=s","m=s","path=s");

my $usage = <<"USAGE";
        Program : $0
        Contact : shaohua.yuan\@majorbio.com
        Discription: Calculate nr taxon abundance.
        Usage:perl $0 [options]
                -q     input gene_kegg_anno detail file
                -p     input gene profile file
                -e     input enzyme profile file
                -m     input module profile file
                -path   input pathway profile file
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

my (%enzyme_des,%module_des);
open INFE,$opts{e} or die "can not open $opts{e}!\n";
while(<INFE>){
    chomp;
    next if(/^#/);
    my @tmp = split /\t/;
    my $enzyme = $tmp[0];
    my $e_des = $tmp[-1];
    if (!exists $enzyme_des{$enzyme}){
        $enzyme_des{$enzyme} = $e_des;
    }
}
close INFE;

open INFM,$opts{m} or die "can not open $opts{m}!\n";
while(<INFM>){
    chomp;
    next if(/^#/);
    my @tmp = split /\t/;
    my $module = $tmp[0];
    my $m_des = $tmp[-1];
    if (!exists $module_des{$module}){
        $module_des{$module} = $m_des;
    }
}
close INFM;

my (%Gene,%gene_KO,%query_KO,%KO,%KO_des,%KO_link,%pathway,%path_KO,%path_des,%enzyme,%module,%level1,%level2,%level3,%l2_l1,%l3_l1,%l3_l2);
open INF,$opts{q} or die "can not open $opts{q}!\n";
while(<INF>){
    chomp;
    next if(/^#/);
    my @tmp = split /\t/;
    my $query = $tmp[0];
    if(exists $profile_genes{$query}){
        my $Gene = $tmp[1];
        my $KO = $tmp[2];
        $gene_KO{$Gene} = $KO;
        $query_KO{$query} = $KO;
        my $KO_des = $tmp[3];
        my $pathway = $tmp[4];
        my $enzyme = $tmp[5];
        my $module = $tmp[6];
        my $link = $tmp[7];
        my $level1 = $tmp[-3];
        my $level2 = $tmp[-2];
        my $level3 = $tmp[-1];
        $KO_des{$KO} = $KO_des;
        $KO_link{$KO} = $link;
        my @pathways = split(/;/,$pathway);
        my @l1 = split(/;/,$level1);
        my @l2 = split(/;/,$level2);
        my @l3 = split(/;/,$level3);
        for(my $i = 0;$i < @pathways;$i++){
            my $path = $pathways[$i];
            $pathway{$path} .=",".$query;
            $path_KO{$path} .=",".$KO;
            $path_des{$path} = $l3[$i];
            $level1{$l1[$i]} .=",".$query;
            $level2{$l2[$i]} .=",".$query;
            $level3{$l3[$i]} .=",".$query;
            $l2_l1{$l2[$i]} = $l1[$i];
            $l3_l1{$l3[$i]} = $l1[$i];
            $l3_l2{$l3[$i]} = $l2[$i];
        }
        $Gene{$Gene} .= ",".$query;
        $KO{$KO} .= ",".$query;
        my @enzymes = split(/;/,$enzyme);
        foreach my $enzy (@enzymes){
            $enzyme{$enzy} .= ",".$query;
        }
        my @module = split(/;/,$module);
        foreach my $modu (@module){
            $module{$modu} .= ",".$query;
        }
    }
}
close INF;


my $samples = join("\t",@sams);
open OUTG, "> $opts{o}/kegg_gene_profile.xls" or die "can not open $opts{o}/kegg_gene_profile.xls!\n";
print OUTG "#Gene\t","$samples\tKO\n";
foreach my $eachgene (sort keys %Gene){
        print OUTG "$eachgene\t";
        my $Gene_gene = $Gene{$eachgene};
        my @genes = split(/,/,$Gene_gene);
        shift @genes;
        for(my $i = 0;$i < @sams ;$i++){
            my $Gene_abu;
            foreach my $gene (@genes){
                my $sam_profile;
                my @sams_profile = split(/\t/,$gene_profile{$gene});
                $sam_profile = $sams_profile[$i];
                $Gene_abu += $sam_profile;
           }
        print OUTG $Gene_abu,"\t";
        }
        print OUTG "$gene_KO{$eachgene}\n";
}
close OUTG;

my %KO_abu;
open OUTK, "> $opts{o}/kegg_KO_profile.xls" or die "can not open $opts{o}/kegg_KO_profile.xls!\n";
print OUTK "#KO\t","$samples\tDescription\tHyperlink\n";
foreach my $eachKO (sort keys %KO){
    next if($eachKO ~~ "-");
    print OUTK "$eachKO\t";
    my $KO_gene = $KO{$eachKO};
    my @genes = split(/,/,$KO_gene);
    shift @genes;
    for(my $i = 0;$i < @sams ;$i++){
        my $KO_abu;
        foreach my $gene (@genes){
            my $sam_profile;
            if(exists $gene_profile{$gene}){
                my @sams_profile = split(/\t/,$gene_profile{$gene});
                $sam_profile = $sams_profile[$i];
            }else{
                $sam_profile = 0;
            }
            $KO_abu += $sam_profile;
       }
    print OUTK $KO_abu,"\t";
    $KO_abu{$eachKO} += $KO_abu/2;   ### 包含Total列，除以2
    }
    print OUTK "$KO_des{$eachKO}\t$KO_link{$eachKO}\n";
}
close OUTK;

open OUTE, "> $opts{o}/kegg_enzyme_profile.xls" or die "can not open $opts{o}/kegg_enzyme_profile.xls!\n";
print OUTE "#Enzyme\t","$samples\tDescription\n";
foreach my $eachenzyme (sort keys %enzyme){
    next if($eachenzyme ~~ "-");
    print OUTE "$eachenzyme\t";
    my $enzyme_gene = $enzyme{$eachenzyme};
    my @genes = split(/,/,$enzyme_gene);
    shift @genes;
    for(my $i = 0;$i < @sams ;$i++){
        my $enzyme_abu;
        foreach my $gene (@genes){
            my $sam_profile;
            if(exists $gene_profile{$gene}){
                my @sams_profile = split(/\t/,$gene_profile{$gene});
                $sam_profile = $sams_profile[$i];
            }else{
                $sam_profile = 0;
            }
            $enzyme_abu += $sam_profile;
       }
    print OUTE $enzyme_abu,"\t";
    }
    print OUTE "$enzyme_des{$eachenzyme}\n";
}
close OUTE;

open OUTM, "> $opts{o}/kegg_module_profile.xls" or die "can not open $opts{o}/kegg_module_profile.xls!\n";
print OUTM "#Module\t","$samples\tDescription\n";
foreach my $eachmodule (sort keys %module){
    next if($eachmodule ~~ "-");
    print OUTM "$eachmodule\t";
    my $module_gene = $module{$eachmodule};
    my @genes = split(/,/,$module_gene);
    shift @genes;
    for(my $i = 0;$i < @sams ;$i++){
        my $module_abu;
        foreach my $gene (@genes){
            my $sam_profile;
            if(exists $gene_profile{$gene}){
                my @sams_profile = split(/\t/,$gene_profile{$gene});
                $sam_profile = $sams_profile[$i];
            }else{
                $sam_profile = 0;
            }
            $module_abu += $sam_profile;
       }
    print OUTM $module_abu,"\t";
    }
    #print "$eachmodule\t$module_des{$eachmodule}\n";
    print OUTM "$module_des{$eachmodule}\n";
}
close OUTM;

open OUTP, "> $opts{o}/kegg_pathway_profile.xls" or die "can not open $opts{o}/kegg_pathway_profile.xls!\n";
open OUTPE, "> $opts{o}/kegg_pathway_eachmap.xls" or die "can not open $opts{o}/kegg_pathway_eachmap.xls!\n";
print OUTP "#Pathway\t","$samples\tDescription\tPathwaymap\n";
my @sams1 = @sams;
pop @sams1;
my $sam_map = join("_map\t",@sams1)."_map";
print OUTPE "Pathway\tDescription\tKO_list\tTotal_abundance\tQuery_count\t$sam_map\n";
foreach my $eachpathway (sort keys %pathway){
    my @pathway_map;
    next if($eachpathway ~~ "-");
    my $total_abun = 0;
    print OUTP "$eachpathway\t";
    print OUTPE "$eachpathway\t$path_des{$eachpathway}\t";
    my $KOs = $path_KO{$eachpathway};
    my @KOs = split(/,/,$KOs);
    shift @KOs;
    my %KO_count;
    @KOs = grep { ++$KO_count{$_} < 2 } @KOs;
    foreach my $each (@KOs){
    $total_abun += $KO_abu{$each};
    }
    #my $newKOs = join("+",@KOs);
    my $KO_list = join(";",@KOs);
    #my $map = "http://www.genome.jp/kegg-bin/show_pathway?".$eachpathway."+".$newKOs ;
    my $pathway_gene = $pathway{$eachpathway};
    my @genes = split(/,/,$pathway_gene);
    shift @genes;
    my $query_count = @genes;
    print OUTPE "$KO_list\t$total_abun\t$query_count";
    #print OUTPE "$query_count\t$total_abun\t$KO_list";
    for(my $i = 0;$i < @sams ;$i++){
        my $pathway_abu;
        my @sam_KO;
        my $eachmap;
        foreach my $gene (@genes){
            my $sam_profile;
            if(exists $gene_profile{$gene}){
                my @sams_profile = split(/\t/,$gene_profile{$gene});
                $sam_profile = $sams_profile[$i];
            #print "$sam_profile\n";
                if ($sam_profile != 0){
                #print "$gene\t$query_KO{$gene}\n";
                    push (@sam_KO,$query_KO{$gene});
                    push (@pathway_map,$query_KO{$gene});
                #print "sam_profile!=0\t$sam_profile\t@sam_KO\n";
                }
            }else{
                $sam_profile = 0;
            }
            $pathway_abu += $sam_profile;
       }
    if (@sam_KO){
        my %count;
        @sam_KO = grep { ++$count{$_} < 2 } @sam_KO;
        $eachmap = join("+",@sam_KO);
        #print "$sams[$i]\t$eachpathway\teachmap:$eachmap","\n";
        print OUTPE "\thttp://www.genome.jp/kegg-bin/show_pathway?".$eachpathway."+".$eachmap;
    }else{
        print OUTPE "\tNA";
    }
    print OUTP $pathway_abu,"\t";
    }
    if (@pathway_map){
        my %count;
        @pathway_map = grep { ++$count{$_} < 2 } @pathway_map;
        my $pathway_map = join("+",@pathway_map);
        my $map = "http://www.genome.jp/kegg-bin/show_pathway?".$eachpathway."+".$pathway_map ;
        print OUTP "$path_des{$eachpathway}\t$map\n";
    }else{
    print OUTP "$path_des{$eachpathway}\t-\n";
    }
    print OUTPE "\n";
}
close OUTP;
#close OUTI;
close OUTPE;

open OUTLO, "> $opts{o}/kegg_level1_profile.xls" or die "can not open $opts{o}/kegg_level1_profile.xls!\n";
print OUTLO "#Level1\t","$samples\n";
foreach my $eachlevel1 (sort keys %level1){
    next if($eachlevel1 ~~ "-");
    print OUTLO "$eachlevel1";
    my $level1_gene = $level1{$eachlevel1};
    my @genes = split(/,/,$level1_gene);
    my %count;
    my @uniq_genes = grep { ++$count{ $_ } < 2; } @genes;
    shift @genes;
    for(my $i = 0;$i < @sams ;$i++){
        my $level1_abu;
        foreach my $gene (@uniq_genes){
            my $sam_profile;
            if(exists $gene_profile{$gene}){
                my @sams_profile = split(/\t/,$gene_profile{$gene});
                $sam_profile = $sams_profile[$i];
            }else{
                $sam_profile = 0;
            }
            $level1_abu += $sam_profile;
       }
    print OUTLO "\t$level1_abu";
    }
    print OUTLO "\n";
}
close OUTLO;

open OUTLW, "> $opts{o}/kegg_level2_profile.xls" or die "can not open $opts{o}/kegg_level2_profile.xls!\n";
print OUTLW "#Level1\tLevel2\t","$samples\n";
foreach my $eachlevel2 (sort keys %level2){
    next if($eachlevel2 ~~ "-");
    print OUTLW "$l2_l1{$eachlevel2}\t$eachlevel2";
    my $level2_gene = $level2{$eachlevel2};
    my @genes = split(/,/,$level2_gene);
    shift @genes;
    my %count;
    my @uniq_genes = grep { ++$count{ $_ } < 2; } @genes;
    for(my $i = 0;$i < @sams ;$i++){
        my $level2_abu;
        foreach my $gene (@uniq_genes){
            my $sam_profile;
            if(exists $gene_profile{$gene}){
                my @sams_profile = split(/\t/,$gene_profile{$gene});
                $sam_profile = $sams_profile[$i];
            }else{
                $sam_profile = 0;
            }
            $level2_abu += $sam_profile;
       }
    print OUTLW "\t$level2_abu";
    }
    print OUTLW "\n";
}
close OUTLW;

open OUTLT, "> $opts{o}/kegg_level3_profile.xls" or die "can not open $opts{o}/kegg_level3_profile.xls!\n";
print OUTLT "#Level1\tLevel2\tLevel3\t","$samples\n";
foreach my $eachlevel3 (sort keys %level3){
    next if($eachlevel3 ~~ "-");
    print OUTLT "$l3_l1{$eachlevel3}\t$l3_l2{$eachlevel3}\t$eachlevel3";
    my $level3_gene = $level3{$eachlevel3};
    my @genes = split(/,/,$level3_gene);
    shift @genes;
    my %count;
    my @uniq_genes = grep { ++$count{ $_ } < 2; } @genes;
    for(my $i = 0;$i < @sams ;$i++){
        my $level3_abu;
        foreach my $gene (@uniq_genes){
            my $sam_profile;
            if(exists $gene_profile{$gene}){
                my @sams_profile = split(/\t/,$gene_profile{$gene});
                $sam_profile = $sams_profile[$i];
            }else{
                $sam_profile = 0;
            }
            $level3_abu += $sam_profile;
       }
    print OUTLT "\t$level3_abu";
    }
    print OUTLT "\n";
}
close OUTLT;

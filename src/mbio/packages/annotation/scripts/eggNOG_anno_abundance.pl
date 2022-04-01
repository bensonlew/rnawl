#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
my %opts;
GetOptions(\%opts,"q=s","p=s","o=s");

my $usage = <<"USAGE";
        Program : $0
        Contact : shaohua.yuan\@majorbio.com
        Discription: Calculate COG abundance.
        Usage:perl $0 [options]
                -q     input gene_cog_anno detail file
                -p     input gene profile file
                -o     output Dir

USAGE

die $usage  if (!($opts{q}&&$opts{o}&&$opts{p}));

`mkdir -p $opts{o}` unless -e ($opts{o});

#my (%gene_profile,@sams,@profile_genes);
my (%gene_profile,@sams,%profile_genes);
open INP,$opts{p} or die "can not open $opts{p}!\n";
while(<INP>){
    chomp;
    if (/^GeneID/){
    @sams = split /\t/;
    #### ɸ��head��GeneID�к�����total��
    shift @sams;
    #pop @sams;
    }else{
    my @tmp = split /\t/;
    my $gene = $tmp[0];
    #push (@profile_genes,$gene);
    $profile_genes{$gene} = 1;
    shift @tmp;
    my $tmp_profile = join("\t",@tmp);
    $gene_profile{$gene} = $tmp_profile;
    }
}
close INP;

my (%cog,%fun,%cate,%cog_des,%fun_des);
open INF,$opts{q} or die "can not open $opts{q}!\n";
while(<INF>){
    chomp;
    next if(/^#/);
    my (@funs,@fun_deses,@cates);
    # 每行表头
    # #Query,NOG,NOG_description,Function,Fun_description,Category,Identity(%),Align_len
    my @tmp = split /\t/,$_;
    my $query = $tmp[0];
    #if(grep /^$query$/, @profile_genes){
    if(exists $profile_genes{$query}){
        my $cog = $tmp[1];
        my $cog_des = $tmp[2];
        $cog_des{$cog} = "$cog_des\t$tmp[3]\t$tmp[4]\t$tmp[5]";
        my $fun = $tmp[3];
        my $fun_des = $tmp[4];
        $fun_des{$fun} = $fun_des;
        my $cate = $tmp[5];
        $cog{$cog} .= ",".$query;
        if (length($fun) > 1){
            my @funs = split /;/, $fun;
            my @fun_deses = split /;/, $fun_des;
            for(my $i = 0;$i < @funs;$i++){
            $fun{$funs[$i]} .= ",".$query;
            $fun_des{$funs[$i]} = $fun_deses[$i];
            }
            my @cates = split /;/, $cate;
            my %count;
            my @uniq_cates = grep { ++$count{ $_ } < 2; } @cates;
            for my $j (@uniq_cates){
                $cate{$j} .=",".$query;
            }
        }else{
            $fun{$fun}  .= ",".$query;
            $fun_des{$fun} = $fun_des;
            $cate{$cate} .=",".$query;
        }
    }
}
close INF;

open OUTC, "> $opts{o}/cog_nog_profile.xls" or die "can not open $opts{o}/cog_nog_profile.xls!\n";
open OUTT, "> $opts{o}/cog_function_profile.xls" or die "can not open $opts{o}/cog_function_profile.xls!\n";
open OUTG, "> $opts{o}/cog_category_profile.xls" or die "can not open $opts{o}/cog_category_profile.xls!\n";
my $samples = join("\t",@sams) ;
print OUTC "#NOG\t","NOG_description\tFunction\tFun_description\tCategory\t$samples\n";
print OUTT "#Function\t","$samples\t","Description\n";
print OUTG "#Category\t","$samples\n";
foreach my $eachnog (sort keys %cog){
        my $cog_gene = $cog{$eachnog};
        my @genes = split(/,/,$cog_gene);
        shift @genes;
        print OUTC "$eachnog\t$cog_des{$eachnog}";
        for(my $i = 0;$i < @sams ;$i++){
            #print $sams[$i],"\n";
            my $cog_abu;
            foreach my $gene (@genes){
                my @sams_profile = split(/\t/,$gene_profile{$gene});
                my $sam_profile = $sams_profile[$i];
                $cog_abu += $sam_profile;
           }
        print OUTC "\t$cog_abu";
        #print $cog_abu,"\n";
        }
        print OUTC "\n";
}
close OUTC;

foreach my $eachfun (sort keys %fun){
    # print $eachfun,"\n";
    # print $fun_des{$eachfun},"\n";
        print OUTT "$eachfun\t";
        my $fun_gene = $fun{$eachfun};
        my @genes = split(/,/,$fun_gene);
        shift @genes;
        for(my $i = 0;$i < @sams ;$i++){
            my $fun_abu;
            foreach my $gene (@genes){
                my @sams_profile = split(/\t/,$gene_profile{$gene});
                my $sam_profile = $sams_profile[$i];
                $fun_abu += $sam_profile;
           }
        print OUTT $fun_abu,"\t";
        }
        print OUTT "$fun_des{$eachfun}\n";
}
close OUTT;

foreach my $eachcate (sort keys %cate){
        print OUTG "$eachcate\t";
        my $cate_gene = $cate{$eachcate};
        my @genes = split(/,/,$cate_gene);
        shift @genes;
        for(my $i = 0;$i < @sams ;$i++){
            my $cate_abu;
            foreach my $gene (@genes){
                print "$sams[$i]\t$gene","\n";
                my @sams_profile = split(/\t/,$gene_profile{$gene});
                my $sam_profile = $sams_profile[$i];
                $cate_abu += $sam_profile;
           }
        print OUTG $cate_abu,"\t";
        }
        print OUTG "\n";
}
close OUTG;

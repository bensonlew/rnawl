#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my %opts;
GetOptions(
	\%opts,"tax=s","fun=s","tl=s","fl=s","p=s","Tt=i","Tf=i","o=s","help!");

my $usage = <<"USAGE";
              Version : v1.1_20171103
              Writer : shaohua.yuan\@majorbio.com

Usage: perl $0 [options]
                -tax 	 *Gene annotation classification (eg:nr.tax.xls)
                -fun 	 *function annotation  (eg:cog.gene.annotate.xls)
                -tl      *Taxon levels
                -fl      *function level(eg: Function,Pathway)
                -p	     *gene abundance (eg:gene_profile.reads_number.txt)  ## profileΪ��Total�ķ�ȱ�
                -Tt	     TOP of the first column (eg:10).Input "all" when calculate all.
                -Tf      TOP of the second column (eg:10).Input "all" when calculate all.
                -o       *output DIR
                -help	  Display this usage information
         eg:perl $0  -tl s -tax nr.tax.xls -fun eggNOG.anno.xls -fl Function -p reads_number.txt -Tt 10 -Tf 10 -o ourDir

USAGE

die $usage if ( !( $opts{tl} ) || !( $opts{tax} ) || !( $opts{fun} ) || !( $opts{fl} ) || !( $opts{p} )  || !( $opts{o} ));

`mkdir -p $opts{o}` unless -e ($opts{o});
#define defaults
$opts{Tt}=$opts{Tt}?$opts{Tt}:10;
$opts{Tf}=$opts{Tf}?$opts{Tf}:10;

my (%tf_genes,@top_tax,@top_fun,@all_tax,@all_fun,%gene_profile,@sams,%profile_genes);

my ($tax_list,$fun_list,$tax_abu,$fun_abu);
my $tax_list_file = $opts{o}."/Taxon_$opts{tl}.xls";
my $fun_list_file = $opts{o}."/Fun_$opts{fl}.xls";
$tax_abu = $opts{o}."/Tax_$opts{tl}_abu.xls";
$fun_abu = $opts{o}."/Fun_$opts{fl}_abu.xls";
&Tax_fun_gene($opts{tl},$opts{tax},$tax_list_file);
&Tax_fun_gene($opts{fl},$opts{fun},$fun_list_file);

my ($profile,$sams) = &function_abundance($opts{tl},$opts{tax},$opts{p},$tax_abu);
&function_abundance($opts{fl},$opts{fun},$opts{p},$fun_abu);

&get_top($opts{Tt},$opts{Tf},$tax_abu,$fun_abu);

my $final_output1= $opts{o}."/Taxon_function_abundance.xls";
my $final_output2= $opts{o}."/Function_taxon_abundance.xls";
print "start taxon_function abudance\n";
#&Taxon_fuction(\%top_tax,\%top_fun,$opts{tl},$opts{fl},$profile,$sams,$tax_list,$fun_list,$final_output1);
#&Taxon_fuction($top_fun,$top_tax,$opts{fl},$opts{tl},$profile,$sams,$fun_list,$tax_list,$final_output2);

&Taxon_fuction(1,$opts{tl},$opts{fl},$tax_list,$fun_list,$final_output1);
&Taxon_fuction(2,$opts{fl},$opts{tl},$fun_list,$tax_list,$final_output2);

#####���ݷ����ȡ����#####
sub Tax_fun_gene{
    my ($level,$input,$output)=@_;
    open (IN,$input) || die "can not open $input!\n";
    open (OUT,">$output") || die "can not open $output!\n";
    print OUT "#$level\tGene_ID\n";
    my ($leve_index, %tmp_tf_genes);
    while(<IN>){
        chomp;
        if (/^#/){
            my @head = split /\t/;
            ($leve_index)=grep{$head[$_] =~ $level} 0..$#head;
            print $level,"\n";
            print join(",",@head),"\n";
            print $leve_index,"\n";
            }else{
        my @temp=split /\t/;
        my $query = $temp[0];
        my $tax_fun=$temp[$leve_index];
        if ($tax_fun =~ /;/){
            $tax_fun =~ s/;\s+/;/g;
            $tax_fun =~ s/^\s//;
            my @tax_funs = split(/;/,$tax_fun);
            foreach my $each_tf(@tax_funs){
            $tf_genes{$each_tf} .= $query.";";
            $tmp_tf_genes{$each_tf} .= $query.";";
            }
        }else{
        $tf_genes{$tax_fun} .= $query.";";
        $tmp_tf_genes{$tax_fun} .= $query.";";
        }
        }
    }
    close IN;
    foreach (sort  keys %tmp_tf_genes){
        if ( $_ ne "-"){
            print OUT "$_\t$tf_genes{$_}\n";
        }
        }
    close OUT;
}

###��ȼ���
sub function_abundance{
    my ($level,$input,$abundance,$output)= @_;
    open (IN,$abundance) || die "can not open $abundance!\n";
    while(<IN>){
        chomp;
        if (/^GeneID/){
            @sams = split /\t/;
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
    close IN;
    my $samples = join("\t",@sams);
    open (INF,$input) || die $!;
    open (OUT,">$output") || die $!;
    print OUT "$level\t$samples\n";
    my (@arry2,$level_index);
    my %level_genes;
    while(<INF>){
        chomp;
        if (/^#/){
            my @head = split /\t/;
            my $test = join(",",@head);
            print "$test\n";
            ($level_index)=grep{$head[$_] =~ $level} 0..$#head;
        }else{
             my @temp=split /\t/;
             my $level_name = $temp[$level_index];
             if(exists $profile_genes{$temp[0]}){
             my @level_names = split(/;/,$level_name);
             for my $each (@level_names){
                 $each =~ s/^\s+|\s+$//;
                if ( $each ne "-"){
                    $level_genes{$each} .=",".$temp[0];
              }
              }
         }
        }
    }
    close INF;
    foreach my $eachlevel(sort keys %level_genes){
    print OUT "$eachlevel\t";
    my $level_gene = $level_genes{$eachlevel};
    my @genes = split(/,/,$level_gene);
    shift @genes;
    my %count;
    my @uniq_genes = grep { ++$count{ $_ } < 2; } @genes;
    for(my $i = 0;$i < @sams ;$i++){
        my $level_abu;
        foreach my $gene (@uniq_genes){
            my @sams_profile = split(/\t/,$gene_profile{$gene});
            my $sam_profile = $sams_profile[$i];
            $level_abu += $sam_profile;
       }
    print OUT $level_abu,"\t";
    }
    print OUT "\n";
    }
    close OUT;
}



sub get_top{
    my ($top1,$top2,$abu1,$abu2)=@_;
    open (IN,$abu1) ||die "can not open $abu1!\n";
    open (IN2,$abu2) ||die die "can not open $abu2!\n";
    my (%abundance,%abundance2);
    my (@arry,@arry2);
    my $name2=<IN>;
    while(<IN>){
    chomp;
    my @temp=split /\t/;
    $abundance{$temp[0]}=$temp[-1];
    }
    close IN;
    ## ���ܷ������
    foreach (sort {$abundance{$b} <=> $abundance{$a}} keys %abundance){
        push @arry,$_;
    }
    my $tmp_len = @arry;
    if(($top1 eq "all") or ($tmp_len < $top1)){
        $top1 =  $tmp_len;
    }
    for my $i (1..$tmp_len){
        #print "$top1\n";
        #print "$arry[$i-1]\n";
        if ($i <= $top1){
        push @top_tax,$arry[$i-1];
        push @all_tax,$arry[$i-1];
        #print "$arry[$i-1]\n";
        }else{
        push @all_tax,$arry[$i-1];
        }
    }
    my $name=<IN2>;
    my @name=split /\t/,$name;
    chomp($name);
    my $name_len = @name;
    while(<IN2>){
        chomp;
        my @temp=split /\t/;
        $abundance2{$temp[0]}=$temp[-1];
    }
    close IN2;

    foreach (sort {$abundance2{$b} <=> $abundance2{$a}} keys %abundance2){
        push @arry2,$_;
    }
    my $tmp3_len = @arry2;
    if(($top2 eq "all") or ($tmp3_len  < $top2)){
        $top2 =  $tmp3_len;
    }
    for my $i (1..$tmp3_len){
        if ($i <= $top2){
        push @top_fun,$arry2[$i-1];
        push @all_fun,$arry2[$i-1];
        }else{
        push @all_fun,$arry2[$i-1];
        }
    }
}


######����ÿ�����ֵĹ��ܵķ�ȱ���ÿ�����ܵ����ַ�ȱ�#####
sub Taxon_fuction{
    my ($type,$level1,$level2,$list1,$list2,$output)=@_;
    open (OUT1,">$output") ||die die "can not open $output!\n";
    print "$output\n";
    my (@top_1,@top_2,@all);
    if ($type == 1){ @top_1 = @top_tax;@top_2 = @top_fun;@all = @all_fun; }
    elsif($type == 2){ @top_1 = @top_fun;@top_2 = @top_tax; @all = @all_tax; }
    my $top_len = @top_2;
    my $all_fun_len = @all;
    my $sam_names = join("\t",@sams);
    print OUT1 "$level1\t$level2\t$sam_names\n";
    for my $i (1..@top_1){
        my %others_sam_abu;
        my $each_top = $top_1[$i-1];
        #print "$each_top\n";
        my @genes1 = split(/;/,$tf_genes{$each_top});
        my %hash_genes1 = map{$_=>1} @genes1;
        for my $j (1..@all){
            my $each_top2 = $all[$j-1];
            my @common_genes;
            my @genes2 = split(/;/,$tf_genes{$each_top2});
            @common_genes = grep {$hash_genes1{$_}} @genes2;
            my $common_len = @common_genes;
            if ($common_len>0){
               if($j <= $top_len){
                   print OUT1 "$each_top\t$each_top2";
                   for(my $i = 0;$i < @sams ;$i++){
                   my $sam_abu;
                   for my $eachgene (@common_genes){
                        my @sams_profile;
                        my $sam_profile;
                        if(exists $profile_genes{$eachgene}){
                            @sams_profile = split(/\t/,$gene_profile{$eachgene});
                            $sam_profile = $sams_profile[$i];
                            }else{   $sam_profile = 0;  }
                        $sam_abu += $sam_profile;
                        }
                            print OUT1 "\t$sam_abu";
                        }
                   print OUT1 "\n";
               }else{
                   for(my $i = 0;$i < @sams ;$i++){
                   for my $eachgene (@common_genes){
                        my @sams_profile;
                        my $sam_profile;
                        if(exists $profile_genes{$eachgene}){
                            @sams_profile = split(/\t/,$gene_profile{$eachgene});
                            $sam_profile = $sams_profile[$i];
                            }else{   $sam_profile = 0;  }
                        $others_sam_abu{$sams[$i]} += $sam_profile;
                        #print "-----------------------------------------$others_sam_abu{$sams[$i]}\n";
                        }
                        }
               }
            }else{
                if($j <= $top_len){
                    print OUT1 "$each_top\t$each_top2";
                    for(my $i = 0;$i < @sams ;$i++){
                        my $sam_abu = 0;
                        print OUT1 "\t$sam_abu";
                        }
                       print OUT1 "\n";
                }else{
                    for(my $i = 0;$i < @sams ;$i++){
                        $others_sam_abu{$sams[$i]} += 0;
                        }
                }
            }
        }
        if(@all>$top_len ) {
            print OUT1 "$each_top\tothers";
            for(my $i = 0;$i < @sams ;$i++){
                my $final_other = $others_sam_abu{$sams[$i]};
                print OUT1 "\t$final_other";
            }
            print OUT1 "\n";
        }
    }
}

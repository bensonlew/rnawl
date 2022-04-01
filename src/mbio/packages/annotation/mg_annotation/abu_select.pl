#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
my %opts;
GetOptions(\%opts,"i=s","num=s","pro=s","o=s");

my $usage = <<"USAGE";
        Program : $0
        Contact : shaohua.yuan\@majorbio.com
        Discription: Annotation abudance select.
        Usage:perl $0 [options]
                -i     input lowest_level profile
                -num   input keep,sample num ,abu cutoff.("1,3,500")
                -pro   input keep,abu_proportion cutoff.("1,0.5")
                -o     output file
                -database   others or nr. default:"others"

USAGE

die $usage  if (!($opts{i}&&$opts{o}));

$opts{num}=defined $opts{num}?$opts{num}:"1,0,0";
$opts{pro}=defined $opts{pro}?$opts{pro}:"1,0";
$opts{database}=defined $opts{database}?$opts{database}:"others";
if ($opts{num} =~ /all/){
    $opts{num} ="1,0,0";
}
if ($opts{pro} =~ /all/){
    $opts{pro} ="1,0";
}
#`mkdir -p $opts{o}` unless -e ($opts{o});

my @abu_num = split (/,/,$opts{num});
my $abu_keep_remove = int($abu_num[0]);
my $sam_cutoff = $abu_num[1];
my $abu_cutoff = $abu_num[2];
my (@levels,%level,%sam_level,@head,@sams,$total_index);
open INP,"$opts{i}" or die "can not open $opts{i}!\n";
while(<INP>){
    chomp;
    my @tmp = split /\t/;
    my $level = $tmp[0];
    if($opts{database} =~ /nr/){
        my @tmps = split(/;/,$level);
        $level = $tmps[-1];
        }
    my $num = 0;
    if ( $level =~ /^#/){
        @head = @tmp;
        ($total_index)=grep{$head[$_] eq "Total"} 0..$#head;
    }else{
            my @line = @tmp;
            shift @tmp;
            splice(@tmp,$total_index-1,1);    ## 删除基因名和Total
            #my @digit = grep /^[1-9]\d+\d$|[1-9]\d+\.\d+\d$|0\.\d+[1-9]\d+\d$/,@tmp;
            my @digit = grep /^\d+$|^\d+\.\d+$|0\.\d+$/,@tmp;
            my ($digit_index) = grep{$line[$_] ~~ $digit[0]} 0..$#line ;
            #print $level,"\n";
            #print join(",",@digit),"\n";
            my $sam_num = @digit;
            @sams = @head[$digit_index..($digit_index + $sam_num-1)] ;
            #print join(",",@sams),"\n";
            #print join(",",@tmp[0..($sam_num-1)]),"\n";
            my @match = grep{ $_ >= $abu_cutoff } @line[($digit_index)..($digit_index + $sam_num-1)];
            $num = @match;
            if ($abu_keep_remove ==1){
                if ($num >= $sam_cutoff ){
                push  (@levels,$level);
                for(my $i = 0;$i < @sams ;$i++){
                next if $level{$level}{$sams[$i]};  # 跳过重复行的加和 by xieshichang
                $level{$level}{$sams[$i]} = $tmp[$digit_index + $i - 1 ];
                #print "level\t",$level,"\t",$sams[$i],"\n";
                $sam_level{$sams[$i]} += $tmp[$digit_index + $i -1];
                    }
                }
            }elsif($abu_keep_remove == 0){
                if ($num <$sam_cutoff ){
                    push  (@levels,$level);
                    for(my $i = 0;$i < @sams ;$i++){
                        next if $level{$level}{$sams[$i]};  # 跳过重复行的加和 by xieshichang
                        $level{$level}{$sams[$i]} = $tmp[$digit_index + $i -1];
                        $sam_level{$sams[$i]} += $tmp[$digit_index + $i -1];
                    }
                }
            }else{ print "wrong -num fomat!\n"   }
        }
}
close INP;


my @proportion = split (/,/,$opts{pro});
my $pro_keep_remove = $proportion[0];
my $gene_proportion = $proportion[1];
open OUP,"> $opts{o}  " or die "can not open $opts{o}!\n";
print OUP "#name\n";
my $pro;
foreach my $each_level(@levels){
    for(my $i = 0;$i < @sams ;$i++){
    if($sam_level{$sams[$i]}!= 0 ){
    #print $each_level,"\t",$sams[$i],"\n";
    #print $level{$each_level}{$sams[$i]},"\n";
    $pro = $level{$each_level}{$sams[$i]}/$sam_level{$sams[$i]};
    }
    }
    if ($pro_keep_remove == 1){
        if ($pro >= $gene_proportion){
        print OUP $each_level,"\n";
        #print $pro,"\n"
        }
    }elsif($pro_keep_remove == 0){
        if ($pro < $gene_proportion){
            print OUP $each_level,"\n";
        }
    }else{ print "wrong -pro fomat!\n"  }
}



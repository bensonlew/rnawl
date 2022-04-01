#!/usr/bin/perl -w
use strict;

my $script_name = $0;
my $arg;
my $in1;
my $insert;
my $soap;
my $samp;
my $out;
my $ppm_arg; ### 新增ppm 计算方法

while ($arg=shift) {
    if    ($arg eq "-h" ) { print_usage(); }
	elsif ($arg eq "-i1" ) { $in1 = shift; }
	elsif ($arg eq "-i2" ) { $insert = shift; }
	elsif ($arg eq "-i3" ) { $soap = shift; }
	elsif ($arg eq "-n" ) { $samp = shift; }
	elsif ($arg eq "-o" ) { $out = shift; }
	elsif ($arg eq "-ppm" ) { $ppm_arg = shift; }
}

($in1 and $insert and $out and $soap and $samp) || print_usage();
if($ppm_arg){
    $ppm_arg = $ppm_arg;
}else{
   $ppm_arg = "F";
}


my (%lengths,%reads,$reads_sum,$reads_sample,%rpkm,$rpkm_sum,%reads_div,$reads_div_sum,@genename,%fq,%pefq,%sefq);

if($in1=~/\.gz$/){
    open I1,"gzip -dc $in1 |" or die "can't read $in1: $!\n";
}else{
    open I1,$in1 or die "can't read $in1: $!\n";
}

my $in1name=(split /\//,$in1)[-1];
$in1name=~s/\.gz$//;

$/=">";
<I1>;
my $i=0;
while(<I1>){
    $i++;
    chomp;
    my @list=split /\n/;
    my $name=(split /\s+/,(shift @list))[0];
    my $seq=join "",@list;
    my $len=length($seq);
    $lengths{$name}=$len;
    push @genename,$name;
}
close I1;
$/="\n";

my @soaps=split /,/,$soap;
foreach my $pese(@soaps){
if($pese =~ /\.pe$/ or $pese =~ /\.pe\.gz$/){
    if($pese =~ /\.pe$/){
    open PE,$pese or die "can't read $pese: $!\n";
    }elsif($pese =~ /\.pe\.gz$/){
    open PE,"gzip -dc $pese |" or die "can't read $pese: $!\n";
    }
    my @tmp;
    while(<PE>){
    chomp;
    @tmp=split;
    if (not exists($fq{$tmp[0]})){
    if($tmp[3] == 1){
        $reads{$tmp[7]} +=0.5;
        $reads_sample+=0.5;
        $pefq{$tmp[0]}=1;
    }
    }
    }
    %fq=(%fq,%pefq);
    %pefq =();
    close PE;
}elsif($pese =~ /\.se$/ or $pese =~ /\.se\.gz$/){
    if($pese =~ /\.se$/){
    open SE,$pese or die "can't read $pese: $!\n";
    }elsif($pese =~ /\.se\.gz$/){
    open SE,"gzip -dc $pese |" or die "can't read $pese: $!\n";
    }
    print $pese;
    my(@tmp,%past);
    while(<SE>){
    chomp;
    @tmp=split;
    if (not exists($fq{$tmp[0]})){
    if($tmp[3] == 1){
        if($tmp[6] eq '+'){
        if(($lengths{$tmp[7]}-$tmp[8])<$insert+100){
            $past{$tmp[7]}{$tmp[0]}=1;
            $sefq{$tmp[0]} = 1;
        }
        }elsif($tmp[6] eq '-'){
        if(($tmp[8]+$tmp[5])<$insert+100){
            $past{$tmp[7]}{$tmp[0]}=1;
            $sefq{$tmp[0]} = 1;
        }
        }
    }
    }
    }
    %fq=(%fq,%sefq);
    %sefq=();
    close SE;
    foreach my $k(sort keys %past){
    foreach my $read(sort keys %{$past{$k}}){
        if($past{$k}{$read}){
        $reads{$k}++;
        $reads_sample ++;
        }
    }
    }
}
}
%fq = ();
foreach my $gen(@genename){
    $reads{$gen} = 0 unless(defined $reads{$gen});
    $reads_div{$gen} = $reads{$gen}/$lengths{$gen};
    $reads_sum += $reads{$gen};
    $reads_div_sum += $reads_div{$gen};
    $rpkm{$gen} = ($reads{$gen}*(10**9))/($reads_sample*$lengths{$gen});
    $rpkm_sum += $rpkm{$gen};
}
open O2,">$out/$samp\_reads_number.xls" or die "can't write $out/$samp\_reads_number.xls: $!\n";
open O3,">$out/$samp\_reads_number_relative.xls" or die "can't write $out/$samp\_reads_number_relative.xls: $!\n";
open O4,">$out/$samp\_reads_length_ratio_relative.xls" or die "can't write $out/$samp\_reads_length_ratio_relative.xls: $!\n";
open O5,">$out/$samp\_reads_length_ratio.xls" or die "can't write $out/$samp\_reads_length_ratio.xls: $!\n";
open O6,">$out/$samp\_RPKM.xls" or die "can't write $out/$samp\_RPKM.xls: $!\n";
open O7,">$out/$samp\_TPM.xls" or die "can't write $out/$samp\_TPM.xls: $!\n";
if($ppm_arg  eq "T"){
    open O8,">$out/$samp\_PPM.xls" or die "can't write $out/$samp\_PPM.xls: $!\n";
    print O8 "$samp\n";
}

print O2 "$samp\n";
print O3 "$samp\n";
print O4 "$samp\n";
print O5 "$samp\n";
print O6 "$samp\n";
print O7 "$samp\n";
foreach my $gen(@genename){
    my $t2 = $reads{$gen}*2;
    my $t3 = $reads{$gen}/$reads_sum;
    my $t4 = $reads_div{$gen}/$reads_div_sum;
    my $t5 = $reads_div{$gen}*2;
    my $t6 = $rpkm{$gen};
    my $t7 = $t4*(10**6);
    print O2 "$t2\n";
    print O3 "$t3\n";
    print O4 "$t4\n";
    print O5 "$t5\n";
    print O6 "$t6\n";
    print O7 "$t7\n";
    if($ppm_arg eq "T"){
        my $t8 = $reads{$gen}*(10**6) /$reads_sum;
        print O8 "$t8\n";
    }

}
close O2;
close O3;
close O4;
close O5;
close O6;
close O7;
close O8;
%reads=(),%reads_div=(),%rpkm=();

sub print_usage {
        print <<EOD;
Usage: perl $script_name options
    -i1 input gene file, required
    -i2 input information of insertsize
    -i3 soap result, required ,splited by ,
    -n  sample name
     -o output directory, required
     -ppm add ppm method,T or F,default: F
     -h print this help

EOD
exit;
}

#!/usr/bin/perl -w
use List::Util qw/sum/;
use strict;

my $script_name = $0;

my $arg;
my $in1;
my $in2;
my $out;

while ($arg=shift) {
        if    ($arg eq "-h" ) { print_usage(); }
	elsif ($arg eq "-i1" ) { $in1 = shift; }
	elsif ($arg eq "-i2" ) { $in2 = shift; }
	elsif ($arg eq "-o" ) { $out = shift; }
}

($in1 and $in2 and $out) || print_usage();


my (%lengths,%reads,$reads_sum,$reads_sample,%rpkm,$rpkm_sum,%reads_div,$reads_div_sum,@genename,@samps,%fq,%pefq,%sefq);

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

open O1,">$out/head.txt" or die "can't write $out/head.txt: $!\n";
print O1 "GeneID\n";
foreach my $gen(@genename){
   print O1 "$gen\n";
}
close O1;

open I2,$in2 or die "can't read $in2: $!\n";
<I2>;
while(<I2>){
    chomp;
    my ($samp,$insert,$soap)=(split /\s+/)[0,1,2];
    push @samps,$samp;
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
    }
    close O2;
    close O3;
    close O4;
    close O5;
    close O6;
    close O7;
    %reads=(),%reads_div=(),%rpkm=();
}
@genename =();
open F3,">$out/$in1name.length.txt" or die "can't write $out/$in1name.length.txt: $!\n";
print F3 "GeneID\tgene_length\n";
my @files = ("reads_number.xls","reads_number_relative.xls","reads_length_ratio_relative.xls","reads_length_ratio.xls","RPKM.xls","TPM.xls");
my $n =1;
foreach my $file(@files){
    my $paste = "$out/head.txt";
    foreach my $samp(@samps){
        $paste = $paste." $out/$samp\_$file";
        }
    `paste $paste > $out/bak_$file`;
    open F1,"$out/bak_$file" or die "can't read $out/bak_$file: $!\n";
    open F2,">$out/$file" or die "can't write$out/$file: $!\n";
    print F2 "GeneID\t".join("\t", @samps)."\tTotal\n";
    readline F1; # skip the first line
    while (my $line=<F1>){
    chomp($line);
    my @list = split /\s+/,$line;
    my $sum = sum @list[1..$#list];
    if ($sum <= 0){
         next;
    }else{
         print F2 "$line\t$sum\n";
         if($n == 1){
         print F3 "$list[0]\t$lengths{$list[0]}\n";
         }
    }
    }
    $n += 1;
    close F1, close F2,close F3;
}

my $num = @samps + 2;
`head -n 1 $out/reads_number.xls > $out/top100_reads_number.xls`;
`head -n 1 $out/reads_number_relative.xls > $out/top100_reads_number_relative.xls`;
`awk '{if(NR!=1) print }' $out/reads_number.xls | sort -gr -k $num | head -n 100 >> $out/top100_reads_number.xls`;
`awk '{if(NR!=1) print }' $out/reads_number_relative.xls |sort -gr -k $num  | head -n 100 >> $out/top100_reads_number_relative.xls`;

sub print_usage {
        print <<EOD;
Usage: perl $script_name options
    -i1 input gene file, required
    -i2 input information of insertsize and soap result, required
     -o output directory, required
     -h print this help

EOD
exit;
}

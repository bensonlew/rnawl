#!/usr/bin perl 
use strict;
use warnings;

my ($abunTable,$corrTable,$outpre);
my $opt;
while($opt = shift){
        if($opt eq "-i"){
                $abunTable = shift;
        }elsif($opt eq "-m"){
                $corrTable = shift;
        }elsif($opt eq "-o"){
                $outpre = shift;
        }elsif($opt eq "-v"){
                print "Version: V1.201509\n";
                exit;
        }elsif($opt eq "-h"){
                &usage();
                exit;
        }else{
                print STDERR "unknown parameter $opt\n";
                &usage();exit;
        }
}

unless($abunTable and $corrTable and $outpre){
        &usage();
        exit;
}

if ($outpre =~/\//){
my $outDir= $outpre  ;
$outDir=~s/\/[^\/]*$//g;
`mkdir -p $outDir`;
}

my (%go1_gene, %go2_gene, %go3_gene, %go4_gene, %go2_go1, %go3_go1, %go4_go1);
open INA, "< $corrTable" or die "$!\n";
while (<INA>){
	chomp;
	my @line=split(/\t/);
	unless (exists $go2_go1{$line[1]}){
		$go2_go1{$line[1]}=$line[0];
	}elsif ($go2_go1{$line[1]} ne $line[0]){
		print "warning: $line[1] has more than one upper level: $go2_go1{$line[1]};$line[0]\n";
	}
	unless (exists $go3_go1{$line[3]}){
		$go3_go1{$line[3]}=$line[0];
	}elsif ($go3_go1{$line[3]} ne $line[0]){
		print "warning: $line[3] has more than one upper level: $go3_go1{$line[3]};$line[0]\n";
	}
	unless (exists $go4_go1{$line[5]}){
		$go4_go1{$line[5]}=$line[0];
	}elsif ($go4_go1{$line[5]} ne $line[0]){
		print "warning: $line[5] has more than one upper level: $go4_go1{$line[5]};$line[0]\n";
	}

	$line[-1]=~s/_1;/;/g;
	$line[-1]=~s/_1$//;
	my @t=split(/;/, $line[-1]);
	foreach (@t){
		$go1_gene{$line[0]}{$_}=1;
		$go2_gene{$line[1]}{$_}=1;
		$go3_gene{$line[3]}{$_}=1;
		$go4_gene{$line[5]}{$_}=1;
	}
}
close INA;
print "aaaa\n";

my (%go1_sam_num, %go2_sam_num, %go3_sam_num, %go4_sam_num);
open INB, "< $abunTable" or die "$!\n";
my $head = <INB>;
chomp $head;
$head =~ s/^GeneID//g;
while (<INB>){
	chomp;
#	my @temp=split(/\t/);
print "bbb\n";
	%go1_sam_num = &sam_num(\%go1_gene, $_, %go1_sam_num);
print "ccc\n";
	%go2_sam_num = &sam_num(\%go2_gene, $_, %go2_sam_num);
print "dddd\n";
	%go3_sam_num = &sam_num(\%go3_gene, $_, %go3_sam_num);
print "eeee\n";
	%go4_sam_num = &sam_num(\%go4_gene, $_, %go4_sam_num);
print "fffff\n";
}
close INB;


open O1, "> $outpre.level1.xls" or die "Fail to create file $outpre.level1.xls\n";
open O2, "> $outpre.level2.xls" or die "Fail to create file $outpre.level1.xls\n";
open O3, "> $outpre.level3.xls" or die "Fail to create file $outpre.level1.xls\n";
open O4, "> $outpre.level4.xls" or die "Fail to create file $outpre.level1.xls\n";
open O5, "> $outpre.level2.info.xls" or die "Fail to create file $outpre.level2.info.xls\n";
open O6, "> $outpre.level3.info.xls" or die "Fail to create file $outpre.level3.info.xls\n";
open O7, "> $outpre.level4.info.xls" or die "Fail to create file $outpre.level4.info.xls\n";
print "aa11\n";
my $out1 = &trav_print("",%go1_sam_num);
my $out2 = &trav_print("",%go2_sam_num);
my $out3 = &trav_print("",%go3_sam_num);
my $out4 = &trav_print("",%go4_sam_num);
my $out5 = &trav_print(\%go2_go1,%go2_sam_num);
my $out6 = &trav_print(\%go3_go1,%go3_sam_num);
my $out7 = &trav_print(\%go4_go1,%go4_sam_num);
print "aa22\n";
print O1 "Level1$head\n$out1";
print O2 "Level2$head\n$out2";
print O3 "Level3$head\n$out3";
print O4 "Level4$head\n$out4";
print O5 "Level1\tLevel2$head\n$out5";
print O6 "Level1\tLevel3$head\n$out6";
print O7 "Level1\tLevel4$head\n$out7";
print "aa33\n";
close O1;
close O2;
close O3;
close O4;
close O5;
close O6;
close O7;


sub trav_print{
	my ($corr, %hash)=@_;
	my $str="";
	foreach my $key1 (sort keys %hash){
        	if($corr ne ""){
			$str.= "$$corr{$key1}\t$key1";
		}else{
			$str.= "$key1";
		}
	
		foreach my $key2 (sort {$a<=>$b} keys %{$hash{$key1}}){
			my $num=sprintf("%.3f", $hash{$key1}{$key2});
                       	$str.= "\t".$num;
		}
		$str.= "\n";   
	}
	return $str;
}


sub sam_num{
	my ($hash, $s,%sam_num)=@_;
	my @temp=split(/\t/,$s);
         pop @temp;
	my $i=0;
	#while(my ($key,$value) = each(%hash)){print "$key => $value\n";}
	foreach my $g (keys %$hash){
		if (exists $hash->{$g}->{$temp[0]}){
                if ($hash->{$g}->{$temp[0]} == 1){
                        for ($i = 1 ; $i < @temp; $i++){
				unless(exists $sam_num{$g}{$i}){
					$sam_num{$g}{$i}=$temp[$i];
				}else{
                                	$sam_num{$g}{$i}+=$temp[$i];
				}
                        }
                }
		}
        }
	return %sam_num;
}


sub usage{
        print <<EOD
        Description: generate GO profile of 4 GO levels
        Version: V1.201509
        Contact: linmeng.liu\@majorbio.com

        usage: perl $0 -i gene.TMM.fpkm.matrix -m GO.list.level234.xls  -o output
                -i      * table of gene abundance (FPKM),required
                -m      * indo of 4 levels Go and the related genes
                -o      * output prefix,required
EOD
}

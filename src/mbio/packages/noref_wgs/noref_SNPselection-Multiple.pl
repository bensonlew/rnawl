#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($vcf,$out,$tag,$key,$config);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"vcf:s"=>\$vcf,
	"out:s"=>\$out,
	"tag:s"=>\$tag,
	"key:s"=>\$key,
	"config:s"=>\$config,
			) or &USAGE;
&USAGE unless ($vcf and $out and $key);
mkdir $out if(!-d $out);

my(@Alle_Number,$Model,@group,$gnum);
my %sample;
my %group;
my $sround;#sample's rounds
my $ground;#group's rounds
my $subgroup;
my $test;
open IN,$config;
while(<IN>){
	chomp;
	next if($_ eq ""|| /^$/|| /^#/|| /^\{/|| /^\}/);
	if($_=~/Alle_Number/){
		@Alle_Number = split(/,/,(split(/\:/,$_))[-1]);
	}
	next if($_=~/Model/);
	next if($_=~/Selecttion/);
	if($_=~/sampleinfo/){
		$sround++;
		#$test++;
		#print $sround,"\n";
		next;
	}
		if($_ =~ /sample1_name:(\S+)/){
			$test++;
			$sample{$sround}{s1name}=$1;next;
		}
		if($_ =~ /sample1_type:(\w+)/){
			$sample{$sround}{s1type}=$1;next;
		}
		if($_ =~ /sample1_depth:(.*),(.*)/){
			$sample{$sround}{s1min}=$1;
			$sample{$sround}{s1max}=$2;
			next;
		}
		if($_ =~ /sample2_name:(\S+)/){
			$test++;
			$sample{$sround}{s2name}=$1;next;
		}
		if($_ =~ /sample2_type:(\w+)/){
			$sample{$sround}{s2type}=$1;next;
		}
		if($_ =~ /sample2_depth:(.*),(.*)/){
			$sample{$sround}{s2min}=$1;
			$sample{$sround}{s2max}=$2;
			next;
		}
		if($_ =~ /sampleresult:(\w+)/){
			$sample{$sround}{result}=$1;
		}

	if($_=~/groupinfo/){
		$ground++;
		$subgroup=0;
		next;
	}
		if($_ =~ /group_name:(\S+)/){
			$subgroup++;
			$test++;
			$group{$ground}{$subgroup}{name}=$1;
			next;
		}
		if($_ =~ /group_sample:(\S+)/){
			my @groups=split(/,/,$1);
			$gnum+=scalar@groups;
			$group{$ground}{$subgroup}{samples}=$1;next;
		}
		if($_ =~ /group_depth:(.*),(.*)/){
			$group{$ground}{$subgroup}{mindp}=$1;
			$group{$ground}{$subgroup}{maxdp}=$2;
			next;
		}
		if($_ =~ /group_af:(.*),(.*)/){
			$group{$ground}{$subgroup}{minaf}=$1;
			$group{$ground}{$subgroup}{maxaf}=$2;
			next;
		}
		if($_ =~ /group_miss:(.*)/){
			$group{$ground}{$subgroup}{miss}=$1;next;
		}
}
close IN;
#print $sround,"\n",$gnum,"\n";die;
#foreach my$sround(sort{$a<=>$b}keys%sample){
#	print join("\n","$sround\nname:$sample{$sround}{s1name}","type:$sample{$sround}{s1type}","dp:$sample{$sround}{s1min}\t$sample{$sround}{s2max}"),"\n";
#}
#foreach my$ground(sort{$a<=>$b}keys%group){
#	print "$ground\n";
#	foreach my$subgroup(sort{$a<=>$b}keys%{$group{$ground}}){
#		print join("\n","name:$group{$ground}{$subgroup}{name}","samples:$group{$ground}{$subgroup}{samples}","dp:$group{$ground}{$subgroup}{mindp}\t$group{$ground}{$subgroup}{maxdp}","af:$group{$ground}{$subgroup}{minaf}\t$group{$ground}{$subgroup}{maxaf}","miss:$group{$ground}{$subgroup}{miss}"),"\n\n";
#	}
#}

my %seq;
$/=">";
open SEQ,$tag;
while(<SEQ>){
	chomp;
	next if($_ eq ""|| /^$/|| /^#/);
	my($seqid,$squence,$undi)=split(/\n/,$_);
	$seq{$seqid}=$squence;
}
close SEQ;

$/="\n";
my%stat;
my @groups;
my @indi;
my %snp;
my $fIn;
open In,$vcf;
open OUT,">$out/$key.table.xls";
print OUT "Tag ID\tPos\tSNP ID\tREF\tALT\t";
open OUT2,">$out/$key.vcf";
foreach my$sround(sort{$a<=>$b}keys%sample){
	print OUT "$sample{$sround}{s1name} Genotype\t$sample{$sround}{s1name} Allele Depth\t$sample{$sround}{s2name} Genotype\t$sample{$sround}{s2name} Allele Depth";
}
foreach my$ground(sort{$a<=>$b}keys%group){
	foreach my$subgroup(sort{$a<=>$b}keys%{$group{$ground}}){
		print OUT "\t$group{$ground}{$subgroup}{name} Genotype Frequency\t$group{$ground}{$subgroup}{name} Average Depth";
	}
}
print OUT "\tTag seq\n";
#open STAT,">$out/$key.stat";
#print STAT "Analysis ID\tSNP Number\tAverage Depth\tMiss Ratio\n";
my($outsnp,$outnum,$outad);
while(<In>){
	chomp;
	next if($_ eq ""|| /^$/|| /^##/);
	my ($chr,$pos,$id,$ref,$alt,$qual,$Filter,$indo,$format,@info)=split(/\t/,$_);#CHROM  POS     ID      REF     ALT     QUAL	FILTER  INFO    FORMAT  A01     A02     A03
	if (scalar @indi ==0) {
		push @indi,@info;
		#print join("\n",scalar@indi,@info),"\n";die;
		for(my$i=0;$i<scalar@indi;$i++){
			foreach my$sround(keys%sample){
				$sample{$sround}{s1pos}=$i if($indi[$i] eq $sample{$sround}{s1name});
				$sample{$sround}{s2pos}=$i if($indi[$i] eq $sample{$sround}{s2name});
			}
		}
		print OUT2 "$_\n";
		next;
	}else{
		my @alle=split(/\,/,join(",",$ref,$alt));
		my $allenumber=scalar@alle;
		my $d = grep {$_ eq $allenumber }@Alle_Number;#检测是否符合需要的等位基因数
		#print "alle number:$allenumber\n",join("\,",@Alle_Number),"\n",$d,"\n";die;
		next if($d eq "0");
		my $gt;
		my @out;
		my $nowdp;
		my $nownum;
		foreach my$sround(sort{$a<=>$b}keys%sample){
			my @info1=split(/\:/,$info[$sample{$sround}{s1pos}]);
			my @info2=split(/\:/,$info[$sample{$sround}{s2pos}]);
			next if($info1[0] eq "./.");#样本间计较，缺失直接过滤掉
			next if($info2[0] eq "./.");
			next if(($info1[0] eq $info2[0]) && ($sample{$sround}{result} eq "diff"));
			next if(($info1[0] ne $info2[0]) && ($sample{$sround}{result} eq "same"));
			next if($info1[1] < $sample{$sround}{s1min});
			next if($info1[1] > $sample{$sround}{s1max});
			next if($info2[1] < $sample{$sround}{s2min});
			next if($info2[1] > $sample{$sround}{s2max});
			
			my @geno1=split(/\//,$info1[0]);
			my @geno2=split(/\//,$info2[0]);
			next if(($geno1[0] eq $geno1[1]) && ($sample{$sround}{s1type} eq "hete"));
			next if(($geno1[0] ne $geno1[1]) && ($sample{$sround}{s1type} eq "homo"));
			next if(($geno2[0] eq $geno2[1]) && ($sample{$sround}{s2type} eq "hete"));
			next if(($geno2[0] ne $geno2[1]) && ($sample{$sround}{s2type} eq "homo"));
			
			my $genotype1=join("",$alle[$geno1[0]],$alle[$geno1[1]]);
			my $genotype2=join("",$alle[$geno2[0]],$alle[$geno2[1]]);
					
			$nowdp+=$info1[1];
			$nowdp+=$info2[1];
			push @out,"$genotype1\t$info1[2]";
			push @out,"$genotype2\t$info2[2]";#print $indi[$i],"\n",$gt,"\n";die;
			
			foreach my$ground(sort{$a<=>$b}keys%group){
				foreach my$subgroup(sort{$a<=>$b}keys%{$group{$ground}}){
					my $subgroupnumber=0;
					my @samples=split(/,/,$group{$ground}{$subgroup}{samples});#print $subgroup,":",join("\t",@samples),"\n";die;
					my ($tdp,$adp);#miss是组内缺失样本数，tdp是组内样本总深度，adp是组内样本的0基因型的深度
					my $miss||=0;
					for (my $n=0;$n<scalar@info;$n++) {
						my$e = grep { $_ eq $indi[$n] } @samples;
						next if($e eq "0");
						my @infos=split(/\:/,$info[$n]);
						if($infos[0] eq "./."){$miss++;next;}
						my $dp=$infos[1];#print $dp,"\n";die;
						my $alle0=(split(/\//,$infos[0]))[0];
                        my $af=0;
                        if($dp < $group{$ground}{$subgroup}{mindp}){$miss++;next;}
                        if($dp > $group{$ground}{$subgroup}{maxdp}){$miss++;next;}
                        $af=(split(/\,/,$infos[2]))[0] if($alle0 eq "0");
                        $af=$af/$dp;#print $af,"\n";die;
                        if($af < $group{$ground}{$subgroup}{minaf}){$miss++;next;}
                        if($af > $group{$ground}{$subgroup}{maxaf}){$miss++;next;}
                        $adp+=$af if($alle0 eq "0");
                        $tdp+=$dp;
                       	$subgroupnumber++;
                        $nownum++;
					}
					next if($subgroupnumber eq "0");
                    my $nowmiss=(scalar@samples)*$group{$ground}{$subgroup}{miss};#print $nowmiss,"\n";die;
                    #$miss=$miss/scalar@samples;#print $miss,"\n";die;
                    #print "$miss\t$nowmiss\n";die;
                    next if($miss > $nowmiss);
                    $adp=sprintf("%.2f",$adp/$tdp);
                    $tdp=sprintf("%.2f",$tdp/$subgroupnumber);
                    push @out,"$adp\t$tdp";
				}
			}
            next if(scalar@out ne $test);
            my ($tagid,$snppos)=split(/\_/,$id);
            #print join("\t",$tagid,$snppos,$id,@out),"\n";die;
            print OUT join("\t",$tagid,$snppos,$id,$ref,$alt,@out,$seq{$tagid}),"\n";
			print OUT2 $_,"\n";
            #$outsnp++;
            #$outad+=$nowdp;
            #$outnum+=$nownum;
		}
	}
}
close In;
close OUT;
close OUT2;
#print STAT join("\t",$key,$outsnp,sprintf("%.2f",$outad/($outsnp*$sround*2)),sprintf("%.2f",$outnum/($outsnp*$gnum))),"\n";
#close STAT;
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
#my @uniq_times = grep { ++$count{ $_ } < 2; } @array;#数组元素去重
#my %hash_a = map{$_=>1} @a;
#my %hash_b = map{$_=>1} @b;
#my %merge_all = map {$_ => 1} @a,@b;
#my @common = grep {$hash_a{$_}} @a;#交集
#my @a_only = grep {!$hash_b{$_}} @a;#a独有的
#my @b_only = grep {!$hash_a{$_}} @b;#b独有的
#my @merge = keys (%merge_all);     #并集
sub uniq {
  my %seen;
  return grep { !$seen{$_}++ } @_;
  #grep { ++$count{ $_ } < 2; };
}
sub USAGE {#
        my $usage=<<"USAGE";
Contact:        long.huang\@majorbio.com;
Script:			$Script
Description:

Usage:
  Options:
  -vcf	<file>	input vcf list file
  -out	<dir>	output dir
  -tag	<file>	tag file
  -key	<str>	analysis id
  -config	<file>	config file
  -h         Help

USAGE
        print $usage;
        exit;
}

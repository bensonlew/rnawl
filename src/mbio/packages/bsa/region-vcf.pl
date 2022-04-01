#!/usr/bin/env perl -w 
use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Data::Dumper;
use List::Util qw(sum);
my $BEGIN_TIME=time();
my $version="1.0.0";
my ($fIn,$fOut,$wp,$mp,$wb,$mb,$ftype,$select);
GetOptions(
                                "help|?" =>\&USAGE,
                                "i|vcf:s"=>\$fIn,
				"r:s"=>\$select,
                                "out:s"=>\$fOut,
#                               "pid:s"=>\$PID,
#                               "bid:s"=>\$BID,
                                "wp:s"=>\$wp,
                                "mp:s"=>\$mp,
                                "wb:s"=>\$wb,
				"mb:s"=>\$mb,
#                                "pdep:s"=>\$Pdep,
#                                "Bdep:s"=>\$Bdep, # not filter depth; vcf all infomation;
                                "type:s"=>\$ftype, # $ftype=SNP/ALL;
                                ) or &USAGE;
&USAGE unless ($fIn and $fOut and $mb);
$wp||="";$mp||="";$wb||="";
$ftype||="ALL";
$fIn=&ABSOLUTE_DIR($fIn);
$select=&ABSOLUTE_DIR($select);
#############################################################################
my %Indi;
if($mp ne ""){
	$Indi{$mp}="P1";
}
if($wp ne ""){
	$Indi{$wp}="P2";
}
	$Indi{$mb}="B1";
if($wb ne ""){
	$Indi{$wb}="B2";
}
#############################################################################
my %region;
open In,$select;
while (<In>) {
        chomp;
        next if ($_ eq ""||/^$/ ||/#/);
        s/\"//g;
        my ($chr,$pos1,$pos2,undef)=split(/\s+/,$_);
        next if ($chr eq "chr");
        my $regioned=0;
        foreach my $region (sort keys %{$region{$chr}}) {
                my ($pos3,$pos4)=split(/\s+/,$region);
                #print $pos1,"\t",$pos2,"\t",$pos3,"\t",$pos4,"\n";
                if ($pos1 >= $pos3 && $pos1 <= $pos4) {
                        my ($p1,$p2,$p3,$p4)=sort {$a<=>$b} ($pos1,$pos2,$pos3,$pos4);
       	                my $newregion=join("\t",$p1,$p4);
               	        delete $region{$chr}{$region};
                       	$region{$chr}{$newregion}++;
                        $regioned=1;
		}
        }
        if ($regioned == 0) {
                $region{$chr}{join("\t",$pos1,$pos2)}++;
        }
}
close In;
###########################################################################
my @indi;
my %Total;
my %Eff;
#open In,$fIn;
###########################################
if ($fIn =~ /gz$/) {
        open In,"gunzip -dc $fIn|";
}else{
        open In,$fIn;
}
##########################################
open Total,">$fOut.total";
open Eff,">$fOut.eff";
while (<In>) {
        chomp;
        next if ($_ eq ""|| /^$/ || /^##/);
        my ($chr,$pos,$ids,$ref,$alt,$qual,$filter,$info,$format,@geno)=split(/\s+/,$_);
        if (/^#/) {                
		push @indi,@geno;
                my @out;
                foreach my $indi (@indi) {
                        next if (!exists $Indi{$indi});
                        push @out,$indi."-GT"; # save P1-genotype; P1 != P2, save P2 in %Indi;
                        push @out,$indi."-AD"; # save P1-allel-depth;P2
		}
#		print Total "#\@chr\tpos1\tpos2\tsnp\teffsnp\tindel\teffindel\n"; 
		print Total join("\t","#chr","pos","type","reference",@out,"ANNOTATION","HIGH","MODERATE","LOW","MODIFIER"),"\n"; 
#		print Eff "#\@chr\tpos1\tpos2\n";
		print Eff join("\t","#chr","pos","type","reference",@out,"ANNOTATION","HIGH","MODERATE","LOW","MODIFIER"),"\n"; 
        }else{
		my %info;
		my @alle=split(",",join(",",$ref,$alt));
                my %len;
                for (my $i=0;$i<@alle;$i++) {
                        $len{length($alle[$i])}=1;
                }
                my $type="SNP";
                if (scalar keys %len > 1) {
                        $type="INDEL";
                }
                my @format=split(/:/,$format);
                my %ginfo;
                my @outvariant;
                for (my $i=0;$i<@indi;$i++) {
                        next if (!exists $Indi{$indi[$i]});  # values %Indv=P1 P2 B; keys=sample name;
			my $sample=$indi[$i];
                        my @gt=split(/\//,$geno[$i]); # need $i;because @geno=sample 0/0:0,3:...
                        my $id=$Indi{$indi[$i]}; # indi=geno=sample.list; Indv{}=P1 P2 B1 B2;
                        my @info=split(/:/,$geno[$i]); #0/1 col array; @info
                        for (my $j=0;$j<@info;$j++) {
                                $info{$id}{gt}=$info[$j] if ($format[$j] eq "GT");
                                $info{$id}{ad}=$info[$j] if ($format[$j] eq "AD");
                                $info{$id}{dp}=$info[$j] if ($format[$j] eq "DP");
                        }
                        if ($info{$id}{gt} eq "./.") { # ./.
                                $ginfo{$indi[$i]}{gt}="--";
                                $ginfo{$indi[$i]}{ad}="0";
#				print "$ginfo{$indi[$i]}{gt}\t$ginfo{$indi[$i]}{ad}";
                        }else{
                                my ($g1,$g2)=split(/\//,$info{$id}{gt});
                                my @ad=split(/\,/,$info{$id}{ad});  #prepare @ad equal @alle;
                                $ginfo{$sample}{gt}=join("/",$alle[$g1],$alle[$g2]);
                                if ($g1 eq $g2) {
                                        $ginfo{$sample}{ad}=$ad[$g1];
                                }else{
                                        $ginfo{$sample}{ad}=join(",",$ad[$g1],$ad[$g2]);
                                }
                        }
                        push @outvariant,$ginfo{$sample}{gt}; # print sample-gt;
                        push @outvariant,$ginfo{$sample}{ad}; # print sample-ad;
                }
                my @out;
                my %ann;
		my @eff;
                if($info=~/ANN=([^\;]*)/g){
                        my @ann=split(/\,/,$1);
			if($1=~ /\|HIGH\|/ || /\|MODERATE\|/){
				for (my $i=0;$i<@ann;$i++) {
					my @str=split(/\|/,$ann[$i]);
					$str[1]||="";$str[2]||="";$str[3]||="";$str[4]||="";
					my $ann=join("|",$str[1],$str[2],$str[3],$str[4]);
					$ann{$str[2]}++; ##
					push @eff,$ann;
				}
					 
			}
                        for (my $i=0;$i<@ann;$i++) {
                                my @str=split(/\|/,$ann[$i]);
				$str[1]||="";$str[2]||="";$str[3]||="";$str[4]||="";
                                my $ann=join("|",$str[1],$str[2],$str[3],$str[4]);
				$ann{$str[2]}++;
                                push @out,$ann; # save  anno,anno,anno 
			}
                }
		$ann{HIGH}||=0;
		$ann{MODERATE}||=0;
		$ann{LOW}||=0;
		$ann{MODIFIER}||=0;

		next if($ftype =~ /SNP/i && $type =~ /INDEL/i);
		$Total{$chr}{$pos}=join("\t",$type,"$ref/$ref",join("\t",@outvariant),join(";",@out),$ann{HIGH},$ann{MODERATE},$ann{LOW},$ann{MODIFIER});
		$Eff{$chr}{$pos}=join("\t",$type,"$ref/$ref",join("\t",@outvariant),join(";",@eff),$ann{HIGH},$ann{MODERATE},$ann{LOW},$ann{MODIFIER}) if(scalar @eff != 0);
        }

}
close In;
foreach my $chr (sort keys %region) {
	foreach my $region (sort keys %{$region{$chr}}) {
		my ($pos1,$pos2)=split /\t/,$region;
		foreach my $pos (sort {$a<=>$b} keys %{$Total{$chr}}){
			if($pos>=$pos1 && $pos <=$pos2){
				print Total "$chr\t$pos\t$Total{$chr}{$pos}\n";
			}
		}
	}

}
close Total;
foreach my $chr (sort keys %region) {
	foreach my $region (sort keys %{$region{$chr}}) { # region =po1\tpos2;
		my ($pos1,$pos2)=split /\t/,$region;
		foreach my $pos (sort {$a<=>$b} keys %{$Eff{$chr}}){
			if($pos>=$pos1 && $pos <=$pos2){
				print Eff "$chr\t$pos\t$Eff{$chr}{$pos}\n";
			}
		}
	}
}
close Eff;
#######################################################################################
print "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
sub ABSOLUTE_DIR #$pavfile=&ABSOLUTE_DIR($pavfile);
{
        my $cur_dir=`pwd`;chomp($cur_dir);
        my ($in)=@_;
        my $return="";
        if(-f $in){
                my $dir=dirname($in);
                my $file=basename($in);
                chdir $dir;$dir=`pwd`;chomp $dir;
                $return="$dir/$file";
        }elsif(-d $in){
                chdir $in;$return=`pwd`;chomp $return;
        }else{
                warn "Warning just for file and dir \n$in";
                exit;
        }
        chdir $cur_dir;
        return $return;
}

sub USAGE {#
        my $usage=<<"USAGE";
Description: 
Version:  $Script
Contact: long.huang

Usage:
  Options:
        -vcf    <file>  input file 
	-r	select file
        -out    <file>  output file
	-mb	<str>	mut bulk
        -wp	<str>   wild parent (perhaps not have)
	-mp	<str>	mut parent (perhaps not have)
        -wb	<str>   wild bulk 
USAGE
        print $usage;
        exit;
}


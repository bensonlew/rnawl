#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
my $in;
my $in2;
#my $tree="/mnt/clustre/users/danni.li/tree/programe/20170728.txt";
#my $tree="/mnt/ilustre/users/changliang.qiu/danni.li/source/20181027_sourse_tree.txt";//0310
#my $tree="/mnt/ilustre/users/changliang.qiu/danni.li/source/20200127tree.txt";//0505
# my $tree="/mnt/ilustre/users/changliang.qiu/danni.li/source/20200505_tree.txt";//0606
# my $tree="/mnt/ilustre/users/changliang.qiu/danni.li/source/20200606_tree.txt";//0606fix
# my $tree="/mnt/ilustre/users/changliang.qiu/danni.li/source/20200801_tree.txt";//0801
my $tree="/mnt/ilustre/users/changliang.qiu/danni.li/source/20200801_tree_fix.txt";
#my $list1800="/mnt/clustre/users/danni.li/tree/newtree/Y_list1826.txt";
#my $list1800="/mnt/ilustre/users/changliang.qiu/danni.li/source/Y_2205_list.txt";//0310
#my $list1800="/mnt/ilustre/users/changliang.qiu/danni.li/source/Y_20200127_newlist.txt";//0505
# my $list1800="/mnt/ilustre/users/changliang.qiu/danni.li/source/Y_20200505_list.txt";//0606
# my $list1800="/mnt/ilustre/users/changliang.qiu/danni.li/source/Y_20200606_list.txt";//0606fix
# my $list1800="/mnt/ilustre/users/changliang.qiu/danni.li/source/Y_20200801_list.txt";//0801
my $list1800="/mnt/ilustre/users/changliang.qiu/danni.li/source/Y_20200801_list_fix.txt";
GetOptions(
'in=s'=>\$in,
'm=s'=>\$in2,
'tree=s'=>\$tree,
'pos=s' =>\$list1800
);
my $usage=<<USAGE;
-in	<str> input  YFULL_SNP_reslt.xls
-m	<str>	input 1826 SNP result eg:/mnt/clustre/users/hongyu.chen/project/source/Y_full/20171027_re_run_1026/1/database_result/201600068.txt
	-tree *	<str> tree 
	-pos	<str> pos_list
USAGE
if(!$in ||!$tree){
	print "$usage\n";
	exit;
	}
my %list_pos_result;
my %list_pos_SNP;
get_1826list_result();
#merge_result();
sub merge_result{
	 my %sample_list_pos;
	open(IN,"$in");
	open(OUT,">final_all_result.txt");
	print OUT "pos\tSNP\tresult\n";
	while(<IN>){
		chomp;
		my @array=split"\t",$_;
		if(defined $list_pos_result{$array[0]}){
                        $sample_list_pos{$array[0]}++;
                        my @tt=split"\t",$list_pos_result{$array[0]};
			my %name;
			my @name_Yfull=split"/",$array[1];
			my @name_1826=split"/",$list_pos_SNP{$array[0]};
			my $match=0;
			foreach my $k(@name_Yfull){
				$name{$k}++;
				}
			foreach my $k2(@name_1826){
                                if(defined $name{$k2}){
					$match=1;
					}
                                }
			if($match==1){
				print OUT "$array[0]\t$array[1]\t$tt[2]\n";
				}
			elsif($tt[2] ne "not detected" && $array[2] eq "not detected"){
				print OUT "$array[0]\t$array[1]\t-\n";
				}
			else{
				
				print OUT "$array[0]\t$array[1]\t$array[2]\n";
				}
			#my @fname=keys(%name);
			#my $Fname=join"/",@fname;
			#if($array[1]=~/$list_pos_SNP{$array[0]}/){}
			#else{$array[1]="$array[1]/$list_pos_SNP{$array[0]}"}
			#print OUT "$array[0]\t$Fname\t$tt[2]\n";#########both known SNP and 1826 result,print 1826 result####
				
			
			}
		else{
			#######only known SNP#####
			#print OUT "$array[0]\t$array[1]\t$array[2]\t$array[3]\t\t\n";
			print OUT "$array[0]\t$array[1]\t$array[2]\n";
			}
		}
	close IN;
	########only 1826 snp ####
	foreach my $key (keys %list_pos_result){
                if(defined $sample_list_pos{$key}){}
                else{
                         my @tt=split"\t",$list_pos_result{$key};
			print OUT "$key\t$tt[0]\t$tt[2]\n";
			}		
		}
	close OUT;
	}
sub get_1826list_result{
	 my %list_SNP_pos;
        open(IN,"$list1800"||die $!);
         while(<IN>){
                chomp;
                my @array=split"\t",$_;
		$array[1]=~s/\s//g;
                $list_SNP_pos{$array[1]}=$array[0];
		$list_pos_SNP{$array[0]}=$array[1];
                }
        close IN;
	open(OUT2,">final_list_result.txt");
	 print OUT2 "pos\tSNP\tclade\tsample_result\tsample_clade_result\n";
	open(IN,$in2||die $!);
	<IN>;
        while(<IN>){
                chomp;
                my @array=split"\t",$_;
		if(defined $list_SNP_pos{$array[0]}){
			print OUT2 "$list_SNP_pos{$array[0]}\t$_\n";
			}
		else{
			print OUT2 "\t$_\n";
			}
		$list_pos_result{$list_SNP_pos{$array[0]}}=$_;
                }
        close IN;
	close OUT2;
	}

#!usr/bin/perl -w
#use Excel::Writer::XLSX;
use strict;
use warnings;
use Getopt::Long;
use Encode;
no strict "refs";
my $list;
my $dep_Y=1;
my $dep=5;
my $delete_file="yes";
# my $make_vcf="no";
my $sample_type="Y_standard";
my $tree_list;
my $posdb_list;
my $sample;
my $sample_path;
GetOptions(
	# 'list=s'=>\$list,
	'sample_name=s'=>\$sample,
	'sample_path=s'=>\$sample_path,
	'dep_Y=i' =>\$dep_Y,
	'dep=i' =>\$dep,
	'delete_file=s' =>\$delete_file,
	# 'make_vcf=s' =>\$make_vcf,
	'type=s'=>\$sample_type,
	'tree=s'=>\$tree_list,
	'pos=s' =>\$posdb_list,
);
my $usage=<<USAGE;
-sample_name <str> **must
-sample_path <str> **must
-delete_file	<str>	no/yes (yes)
-tree <str> tree_list **must
-pos <str> pos_list **must
USAGE
if(!$sample||!$sample_path||!$posdb_list||!$tree_list){
	print "$usage\n";
	exit;
}
# my @samples;
# my %sample_path;
# my $sample;
#my $p_path="/mnt/clustre/users/danni.li/Programe/Y_sta";
# my $p_path="/mnt/ilustre/users/changliang.qiu/danni.li/Programe/yfull_updata_script";
# open(IN,"$list");
# while(<IN>){
# 	chomp;
# 	my @array=split"\t",$_;
# 	push(@samples,$array[1]);
# 	$sample_path{$array[1]}=$array[0];
# }
# close IN;
my %hash;
my %tree_sourse;
my %tree_snp;
my %snp_tree;
my @major_trees;
my %tree_child;
gettree();
my %pos_snp;
my %pos_mut;
my %snp_pos;
my %pos_ref;
my @all_list_snp;
my @snp_notree;
my @snp_intree;
my %listsnp_tree;
my %list_tree_snp;
getlist();
my %dep1_pos;
my %sample_tree_snp_num;
my %sample_tree_snp;
my %sample_tree_difsnp_num;
my %sample_tree_difsnp;
my $son_nn=0;
my @printout1;
my @inputtrees1;
mkdir($dep_Y);
chdir($dep_Y);
mkdir("database_result");
open(OO,">sample_dep_statistic.txt");
print OO "Y_dep\tdep\tsample\tY_pos_num\tY>dep_pos_num\tchr1-22_pos_num\tchr1-22>dep_pos_num\tsamplematchlistsnpnum\tsnp_in_tree\tsnp_no_tree\tsample_tree_num\tsample_difmut_num\n";
# foreach $sample(@samples){
mkdir($sample);
chdir($sample);
print OO "$dep_Y\t$dep\t$sample\t";
	# if($make_vcf eq "yes"){######if only bam#
	# 	system("/mnt/clustre/users/kefei.huang/bin/samtools mpileup -d 30000 -uf /mnt/clustre/users/kefei.huang/database/HG19/ucsc.hg19.fasta $sample_path{$sample}/2_bwa/align_final_sort.bam -l  $p_path/autosnp2_2.bed -o $sample_path{$sample}/align_final_sort.bcf");
	# 	# system("/mnt/clustre/users/kefei.huang/bin/bcftools call -mA --skip-variants indels -o $sample_path{$sample}/align_final_sort.vcf $sample_path{$sample}/align_final_sort.bcf");//0526(20200505v2)
	# 	###########add at 20200526############
	# 	system("/mnt/clustre/users/kefei.huang/bin/bcftools call -mA -o $sample_path{$sample}/align_final_sort.vcf $sample_path{$sample}/align_final_sort.bcf");
	# }
%dep1_pos=();
if(-e $sample_path){
	depfilter($sample_path,$sample);
	if(-z "raw_variants_fileter_Y.xls"){
	}
	else{
		sample_ana($sample);
		list_dep();
		database_result($sample);
	}
}
else{
	print "$sample no vcf!\n";
}
chdir("../");
all_result();
# }
close OO;
# delete_file();
chdir("../");

	# system("rm list_pos_mut.txt");
	# system("rm major_tree.txt");
	# system("rm notreesnp_list.xls");
# system("perl $p_path/autosome_result_20181228.pl -list $list");//0512
# system("perl $p_path/autosome_result_20200512.pl -list $list");###########add at 20200526############
# system("perl $p_path/autosome_result_20200512v2.pl -list $list");#//20201105
#######################20201105 set for yfull workflow######################
# system("perl  $p_path/add_sourse_database.pl -list $list -result ./$dep_Y/all_sample_result.txt -dep ./$dep_Y/all_sample_ydep_statistc.txt -type $sample_type");
# system("/mnt/ilustre/users/changliang.qiu/Software/miniconda3/envs/python3/bin/python /mnt/ilustre/users/changliang.qiu/Yoogene/Programe/Yfull/add_sourse_database.py -i ./1/all_sample_result.txt");#//20201105

#######################20201105 set for yfull workflow######################
########sub#########
sub database_result{
	my $sample=shift;
	my %clade_SNP;
	my %clade_misSNP;
	my $linenum=0;
	open(OUT,">../database_result/$sample.txt");
	print OUT "SNP\tclade\tsample_all_result\tsample_clade_result\n";
	open(IN,"final_sourse_snps.txt");
	while(<IN>){
		chomp;
		my @array=split"\t",$_;
		if($#array==2 && $array[2] eq "+"){
			$clade_SNP{$array[0]}++;
		}
		elsif($#array==2 && $array[2] eq "-"){
			$clade_misSNP{$array[0]}++;
		}
	}
	close IN;
	open(IN,"final_sample_1826list.txt");
	my @res_list = ();
	while(<IN>){
		chomp;
		my @array=split"\t",$_;
		if($#array==3 && $array[2] ne "SNP"){
			my$SNP=$array[2];
			$array[2]=~s/\*/-/g;
			$SNP=~s/-/*/g;

			$linenum++;
			if(defined $clade_SNP{$array[2]}){
				###########add at 20200526############
				if (grep { $array[2] eq $_ } @res_list){}
				else{
					if(exists $snp_tree{$SNP}){
						print OUT "$array[2]\t$tree_sourse{$snp_tree{$SNP}}\t$array[3]\t+\n";
					}
					else{
						print OUT "$array[2]\t\t$array[3]\t+\n";
					}
				}
				push(@res_list, $array[2]);
			}
			elsif(defined $clade_misSNP{$SNP}){
				###########add at 20200526############s
				if (grep { $array[2] eq $_ } @res_list){}
				else{
					if(exists $snp_tree{$SNP}){
						print OUT "$array[2]\t$tree_sourse{$snp_tree{$SNP}}\t$array[3]\t-\n";
					}
					else{
						print OUT "$array[2]\t\t$array[3]\t-\n";
					}
					push(@res_list, $array[2]);
				}
			}
			else{
				###########add at 20200526############
				if(exists $snp_tree{$SNP}){
					if (grep { $array[2] eq $_ } @res_list){}
					else{
						my $SNP_tree = $snp_tree{$SNP};
						my @SNP_tree_list = split",", $SNP_tree;
						my $SNP_tree_length = scalar(@SNP_tree_list);
						my $SNP_tree_result = '';
						if($SNP_tree_length eq 1){
							$SNP_tree_result = $tree_sourse{$snp_tree{$SNP}};
						}else{
							foreach my$SN (@SNP_tree_list){
								if($SNP_tree_result eq ''){
									$SNP_tree_result = $tree_sourse{$SN};
								}else{
									$SNP_tree_result = join(',',$SNP_tree_result, $tree_sourse{$SN});
								}
							}
						}
						print OUT "$array[2]\t$SNP_tree_result\t$array[3]\t\n";
						push(@res_list, $array[2]);
					}
				}else{
					print OUT "$array[2]\t\t$array[3]\t\n";
				}
					#          print OUT "$array[2]\t$tree_sourse{$snp_tree{$SNP}}\t$array[3]\t\n";
					#         }
					# else{
					#         print OUT "$array[2]\t\t$array[3]\t\n";
					#         }
			}
		}
	}
	close IN;
	close OUT;
}
sub all_result{
	open(OUT,">all_sample_result.txt");
	open(OUT1,">all_sample_ydep_statistc.txt");
	print OUT1 "sample\ty_pos\tdep>=1\tdep>=2\tdep>=5\tdep>=7\tdep>=10\n";
	# foreach my $sample(@samples){
		my $filepath="$ENV{PWD}/$dep_Y/$sample/";
		if(-e "$filepath/final_sourse.txt"){
			open(IN,"$filepath/final_sourse.txt");
			my $line=<IN>;
			chomp($line);
			my @array=split"\t",$line;
			print OUT "$sample\t$array[2]\t$array[1]\n";
			close IN;
		}
		else{
			print OUT "no $sample final_sample_1826list.txt!\n";
		}
		if(-e "$filepath/sample_ydep_statistic.xls"){
			open(IN,"$filepath/sample_ydep_statistic.xls");
			<IN>;
			my $line=<IN>;
			print OUT1 "$line";
			close IN;
		}
		else{
			print OUT1 "no $sample final_sample_1826list.txt!\n";
		}
	# }
	close OUT;
	close OUT1;
}
sub list_dep{
	open(IN,"final_sample_1826list.txt");
	open(OUT,">final_sample_1826list_dep.txt");
	<IN>;
	<IN>;
	my $line=<IN>;
	chomp($line);
	print OUT "$line\tdep\n";
	while(<IN>){
		chomp;
		my @array=split"\t",$_;
		if(defined $dep1_pos{$array[0]}){
			print OUT "$_\t1\n"
		}
		else{
			print OUT "$_\t\n";
		}
	}
	close IN;
	close OUT;
}
# sub delete_file{
# 	my $sample=shift;
# 	if($delete_file eq "yes"){
# 		foreach my$sample(@samples){
# 			my $filepath="$ENV{PWD}/$dep_Y/$sample";
# 			if(-e "$filepath/newfound.xls"){
# 				# system("rm $filepath/newfound.xls");
# 				# system("rm $filepath/new_mut.xls");
# 				# system("rm $filepath/pos-mut_not_same_list.xls");
# 				# system("rm $filepath/pos_match_snp.xls");
# 				# system("rm $filepath/pos_nolist.xls");
# 				# system("rm $filepath/raw_variants_fileter_Y.xls");
# 				# system("rm $filepath/sample1826list.xls");
# 				# system("rm $filepath/sample_difsnp_no_tree.xls");
# 				# system("rm $filepath/sample_difsnp_tree.xls");
# 				# system("rm $filepath/final_sourse.txt");##
# 				# system("rm $filepath/sample_majortree.txt");##
# 				# system("rm $filepath/sample_ydep_statistic.xls");##
# 				# system("rm $filepath/sample_pos_test_result.txt");
# 				# system("rm $filepath/sample_snp_no_tree.xls");
# 				# system("rm $filepath/sample_snp_tree.xls");
# 				# system("rm $filepath/sample_tree_statistic.xls");
# 				# system("rm $filepath/tree_difsnpnum.xls");
# 				# system("rm $filepath/final_tree_result.xls");
# 			}
# 		}
# 	}
# 	else{}
# }

sub sample_ana{
	my $sample_N=shift;
	my @snps=();
	my @allpos=();
	my @sample_alllist_difmut=();
	my %sample_pos_mut=();
	my %sample_inlist_difmut_snp=();
	my $file="raw_variants_fileter_Y.xls";
	open(FILE,$file);
	open(OUT,">pos-mut_not_same_list.xls");
	open(OUT3,">pos_match_snp.xls");
	open(OUT2,">pos_nolist.xls");
	open(OUT1,">new_mut.xls");
	my $pos_without_list=0;
	while(<FILE>){
		chomp;
		my @array=split"\t";
		##############get dep#################
		my %vcf_name_value;
		my @form_name=split";",$array[7];
		foreach my $fn(@form_name){
			my @form_value=split"=",$fn;
			$vcf_name_value{$form_value[0]}=$form_value[1];
		}
		my @form_name2=split":",$array[8];
		my @form_value2=split":",$array[9];
		for(my $num=0;$num<=$#form_name2;$num++){
			$vcf_name_value{$form_name2[$num]}=$form_value2[$num];
		}
		my $nn="DP";
		my $gt="GT";
		##############��������ͻ������
		if($vcf_name_value{$nn}==1){
			if($array[4] eq "<NON_REF>" || $array[4] eq "."){
				push(@{$sample_pos_mut{$array[1]}},$array[3]);
			}
			else{
				my @alt=split",",$array[4];
				foreach my $al(@alt){
					if($al eq "<NON_REF>" || $al eq "."){
						push(@{$sample_pos_mut{$array[1]}},$al);
					}
					else{
						push(@{$sample_pos_mut{$array[1]}},$al);
					}
				}
			}
		}
		##########������>5 ���������Ľ�����ͻ��######
		else{

			if($vcf_name_value{$gt} eq "0/0"){
				push(@{$sample_pos_mut{$array[1]}},$array[3]);
			}
			elsif($vcf_name_value{$gt} eq "0/1"){
				push(@{$sample_pos_mut{$array[1]}},$array[3]);
				my @alt=split",",$array[4];
				foreach my $al(@alt){
					if($al eq "<NON_REF>" || $al eq "."){
					}
					else{
						push(@{$sample_pos_mut{$array[1]}},$al);
					}
				}
			}
			else{
				my @alt=split",",$array[4];
				foreach my $al(@alt){
					if($al eq "<NON_REF>" || $al eq "."){
							}
					else{
						push(@{$sample_pos_mut{$array[1]}},$al);
					}
				}
			}
		}
		my $match=0;
		foreach my $sample_alt(@{$sample_pos_mut{$array[1]}}){
			if(defined($pos_mut{$array[1]})){
				my $pos_length=scalar(@{$pos_mut{$array[1]}});
				for(my $i=0;$i<$pos_length;$i++){
					if ($sample_alt eq @{$pos_mut{$array[1]}}[$i]){
						$match++;
					}
				}
			}
			else{
				print OUT2 "pos:$array[1] no list \n";
			}
		}

		if($match>0){
			 push(@snps,$pos_snp{$array[1]});
			 print OUT3 "$array[1]\t$pos_snp{$array[1]}\t$array[3]\t$array[4]\t@{$sample_pos_mut{$array[1]}}\t@{$pos_mut{$array[1]}}\n";
		}
		else{
			if(defined($pos_snp{$array[1]})){
				@{$sample_inlist_difmut_snp{$array[1]}}=$pos_snp{$array[1]};
				push(@sample_alllist_difmut,$pos_snp{$array[1]});
				print OUT "pos:\t$array[1]\t@{$sample_inlist_difmut_snp{$array[1]}}\t$array[3]\t@{$sample_pos_mut{$array[1]}}\t@{$pos_mut{$array[1]}}\n";
			}
			else{
				$pos_without_list++;
			}
			foreach my $sample_alt(@{$sample_pos_mut{$array[1]}}){
				if($sample_alt eq $array[3]){
				}
				else{
				}
			}
		}

	}
	close FILE;
	close OUT;
	close OUT1;
	close OUT2;
	close OUT3;
	my $snpnum=scalar(@snps);
	my $difmutsnp=scalar(@sample_alllist_difmut);
	print"pos_without_list_num:$pos_without_list\n";
	print"sample snp in list match snp number is $snpnum\n sample snp in list dif nut number is $difmutsnp \n";
	print OO "$snpnum\t";
	##########print 1826 result######
	open(O1,">sample1826list.xls");
	my %detect;
	my %detect_hash;
	@detect{@snps,@sample_alllist_difmut}++;
	my @snpdetect=keys(%detect);
	my $snpdetectnum=scalar(@snpdetect);
	print "detectnum:$snpdetectnum\n";
	foreach my $kk(@snpdetect){
		$detect_hash{$kk}++;
	}

	foreach my$kk( @all_list_snp){
		if(defined($detect_hash{$kk})){
		}
		else {
			print  O1 "$snp_pos{$kk}\t\t$kk\tnot detected\n";
		}
	}

	foreach my $kk(@snps){
		# if(defined ($snp_tree{$kk})){}
		###########add at 20200526############
		my @snp_tree_key = keys(%snp_tree);
		my $snp_n = 0;
		foreach my $ss(@snp_tree_key){
			my @ll_list = split",", $ss;
			if (grep { $kk eq $_ } @ll_list){
				$snp_n = 1;
			}
		}
		if($snp_n eq 0){
			print O1 "$snp_pos{$kk}";
			$kk=~s/IMS\*JST/IMS-JST/g;
			print O1 "\t\t$kk\t+\n";
		}
	}
	foreach my $kk(@sample_alllist_difmut){
		# if(defined($snp_tree{$kk})){}
		###########add at 20200526############
		my @snp_tree_key1 = keys(%snp_tree);
		my $snp_n1 = 0;
		foreach my $ss1(@snp_tree_key1){
			my @ll_list1 = split",", $ss1;
			if (grep { $kk eq $_ } @ll_list1){
				$snp_n1 = 1;
			}
		}
		if($snp_n1 eq 0){
			print O1 "$snp_pos{$kk}";
			$kk=~s/IMS\*JST/IMS-JST/g;
			print O1 "\t\t$kk\t-\n";
		}
	}
	close O1;
	########dif_snp
	my %sample_tree_difsnp=();
	my %sample_tree_difscore=();
	my %sample_tree_difsnp_num=();
	open(OU,">sample_difsnp_tree.xls");
	open(OU1,">sample_difsnp_no_tree.xls");
	my $difsnp_in_tree_num=0;
	my $difsnp_no_tree_num=0;
	foreach my $sample_difsnp (@sample_alllist_difmut){
		if(defined ($snp_tree{$sample_difsnp})){
			$difsnp_in_tree_num++;
			print OU "$sample_difsnp\t$snp_tree{$sample_difsnp}\n";
			my $tree1=$snp_tree{$sample_difsnp};
			###########add at 20200526############
			my @tree1_list = split",",$tree1;
			foreach my $tre(@tree1_list){
				# push(@{$sample_tree_difsnp{$tree1}},$sample_difsnp);
				push(@{$sample_tree_difsnp{$tre}},$sample_difsnp);
			}
		}
		else {
			$difsnp_no_tree_num++;
			print OU1 "$sample_difsnp\tno tree\n";
		}
	}
	close OU;
	close OU1;
	print"sample difsnp $difsnp_in_tree_num in tree , $difsnp_no_tree_num no tree\n";
	###
	open(OUU,">tree_difsnpnum.xls");
	my $sample_dif_tree_num=0;
	foreach my $key(keys(%sample_tree_difsnp)){
		$sample_dif_tree_num++;
		$sample_tree_difsnp_num{$key}=scalar(@{$sample_tree_difsnp{$key}});
		my $dif_score=-($sample_tree_difsnp_num{$key});
		$sample_tree_difscore{$key}=$dif_score;
		print OUU "$key\t$sample_tree_difsnp_num{$key}\t$sample_tree_difscore{$key}\n";
	}
	close OUU;
	###############match_snp########
	my %sample_tree_snp=();
	my %sample_tree_score=();
	my %sample_tree_snp_num=();
	my %tree_snp_num=();
	my %list_tree_snp_num=();
	open(OUT,">sample_snp_tree.xls");
	open(OUT1,">sample_snp_no_tree.xls");
	my $snp_in_tree_num=0;
	my $snp_no_tree_num=0;
	foreach my $sample_snp (@snps){
		if(defined ($snp_tree{$sample_snp})){
			$snp_in_tree_num++;
			print OUT "$sample_snp\t$snp_tree{$sample_snp}\n";
			my $tree1=$snp_tree{$sample_snp};
			push(@{$sample_tree_snp{$tree1}},$sample_snp);
		}
		else {
			$snp_no_tree_num++;
			print OUT1 "$sample_snp\tno tree\n";
		}
	}
	close OUT;
	close OUT1;
	print OO "$snp_in_tree_num\t$snp_no_tree_num\t";
	print"sample snp $snp_in_tree_num in tree, $snp_no_tree_num no tree\n";
	##########################
	open(OUT,">sample_tree_statistic.xls");
	print OUT "tree\tsourse\tmatch_snp\tmatch_snpnum\tmismath_snp_num\n";
	my $sample_tree_num=0;
	foreach my $key(keys(%sample_tree_snp)){
		$sample_tree_num++;
		$sample_tree_snp_num{$key}=scalar(@{$sample_tree_snp{$key}});
		$list_tree_snp_num{$key}=scalar(@{$list_tree_snp{$key}});
		my $score=$sample_tree_snp_num{$key};
		$sample_tree_score{$key}=$score;
		if(defined ($sample_tree_difscore{$key})){}
		else{$sample_tree_difscore{$key}=0;}
		print OUT "$key\t$tree_sourse{$key}\t@{$sample_tree_snp{$key}}\t$sample_tree_snp_num{$key}\t$sample_tree_difscore{$key}\n";
	}
	close OUT;
	print OO "$sample_tree_num\t";
	my $sample_difmut_num=scalar(@sample_alllist_difmut);
	print OO "$sample_difmut_num\n";
	####################################
	open(OUT,">final_tree_result.xls");
	my %sourse_score;
	my %sourse_array;
	my %sourse_level;
	my %sourse_array11=();

	my $maxscore=0;
	my $final_sourse="";
	my $final_tree="";
	print OUT "tree\tmajor_tree_snps\tY_list_snps\tmajor_tree_snps_num\tY_list_snps_num\tsample_match_snps\tsample_match_snps_num\tscore\tsample_mismatch_snp\tdif_snp_num\tdif_snp_score\n";
	foreach my $key(keys(%sample_tree_snp)){
		print OUT "$tree_sourse{$key}\n";
		my $level=0;
		my @sampletree=split"_",$key;
		my $treelength=scalar(@sampletree);
		my $sourse=$tree_sourse{$key};
		my $score=0;
		for(my $i=0;$i<$treelength;$i++){
			my $smalltree=join"_",@sampletree;
			push(@{$sourse_array{$sourse}},$smalltree);

			@{$sourse_array11{$sourse}}=@{$sourse_array{$sourse}};
			if(defined($sample_tree_snp_num{$smalltree})){

			}
			else {
				$sample_tree_snp_num{$smalltree}=0;
			}
			if($sample_tree_snp_num{$smalltree}>0){
				$level++;
			}
			if(defined($sample_tree_difsnp_num{$smalltree})){
			}
			else{
				$sample_tree_difsnp_num{$smalltree}=0;
			}
			$score=$score+$sample_tree_snp_num{$smalltree}-$sample_tree_difsnp_num{$smalltree};

			splice (@sampletree, -1);
		}
		$sourse_level{$sourse}=$level;
		for(my $i=0;$i<$treelength;$i++){
			my $soure_smalltree=pop(@{$sourse_array11{$sourse}});
			my @major_tree_snps=split",",$tree_snp{$soure_smalltree};
			$tree_snp_num{$soure_smalltree}=scalar(@major_tree_snps);
			my @treelistsnps=();
			my $listsnpnum=0;
			my $nolistsnpnum=0;
			foreach my $major_snp(@major_tree_snps){
				if(defined($listsnp_tree{$major_snp})){
					$listsnpnum++;
					push(@treelistsnps,$major_snp);
				}
				else{
					$nolistsnpnum++;
				}
			}
			print OUT ("$soure_smalltree\t$tree_snp{$soure_smalltree}\t@treelistsnps\t$tree_snp_num{$soure_smalltree}\t$listsnpnum\t");
			if (defined($sample_tree_snp{$soure_smalltree})){
				print OUT "@{$sample_tree_snp{$soure_smalltree}}\t$sample_tree_snp_num{$soure_smalltree}\tscore:$sample_tree_score{$soure_smalltree}\t"
			}
			else{
				print OUT "\t0\t0\t";
			}
			if (defined($sample_tree_difsnp{$soure_smalltree})){
				print OUT "@{$sample_tree_difsnp{$soure_smalltree}}\t$sample_tree_difsnp_num{$soure_smalltree}\t $sample_tree_difscore{$soure_smalltree}\n"
			}
			else{
				print OUT "\t0\t0\n";
			}
		}

		$sourse_score{$sourse}=$score;
		if($sourse_score{$sourse}>$maxscore){
			$maxscore=$sourse_score{$sourse};
			$final_sourse=$sourse;
			$final_tree=$key;
		}
		else{
			$maxscore=$maxscore;
			$final_sourse=$final_sourse;
			$final_tree=$final_tree;
		}

	}
	close OUT;
	open(OUT1,">final_sourse.txt");
	open(OUT2,">final_sourse_snps.txt");
	print OUT1 "machine out:$final_tree\t@{$sample_tree_snp{$final_tree}}\t$final_sourse\t$maxscore\tmatchlevel:$sourse_level{$final_sourse}\n\n";
	print"$final_sourse\n";
	print OUT2 "$final_sourse\t\n\n";
	$final_sourse=~/^([A-Z]).*/;
	my$type=$1;
	my @final_sourse_snps;
	print OUT1 "tree\tsourse\tsample_match_snp\tmatch_snp_num(+)\tmismatch_snp_num(-)\n";
	print OUT2 ("SNP\t����Ⱥ\t����\n");
	my $arrleng=scalar(@{$sourse_array{$final_sourse}});
	for(my $i=0;$i<$arrleng;$i++){
		print OUT1 "@{$sourse_array{$final_sourse}}[$i]\t$tree_sourse{@{$sourse_array{$final_sourse}}[$i]}\t";
		if(defined($sample_tree_snp{@{$sourse_array{$final_sourse}}[$i]})){
			print OUT1 "@{$sample_tree_snp{@{$sourse_array{$final_sourse}}[$i]}}\t";
			foreach my $snp1(@{$sample_tree_snp{@{$sourse_array{$final_sourse}}[$i]}}){
				my $matchsnp="$snp1\t$tree_sourse{$snp_tree{$snp1}}\t+\n";
				push(@final_sourse_snps,$matchsnp);
			}
			foreach my $snp1(@{$sample_tree_difsnp{@{$sourse_array{$final_sourse}}[$i]}}){
				my $matchsnp="$snp1\t$tree_sourse{$snp_tree{$snp1}}\t-\n";
				push(@final_sourse_snps,$matchsnp);
			}
		}
		else{
			print OUT1 "\t";
		}
		if(defined($sample_tree_snp_num{@{$sourse_array{$final_sourse}}[$i]})){
			print OUT1 "$sample_tree_snp_num{@{$sourse_array{$final_sourse}}[$i]}\t";
		}
		else{
			print OUT1 "0\t";
		}
		if(defined($sample_tree_difsnp_num{@{$sourse_array{$final_sourse}}[$i]})){
			print OUT1 "$sample_tree_difsnp_num{@{$sourse_array{$final_sourse}}[$i]}\t";
		}
		else{
			print OUT1 "0\t";
		}
		print OUT1 "\n";
	}
	close OUT1 ;
	my $final_sourse_snps_mun=scalar(@final_sourse_snps);
	for(my $i=0;$i<$final_sourse_snps_mun;$i++){
		my $sn=pop(@final_sourse_snps);
		$sn=~s/IMS\*JST/IMS-JST/g;;
		print OUT2 "$sn";
		print "$sn";
	}
	###out put son#

	no strict "refs";
	*get_son_out2 = sub{

		my $final_tree1=shift;
		my @printout_small;
		my $nn1=0;
		@inputtrees1="";
		@printout1="";
		if(@{$tree_child{$final_tree1}}){
			foreach my $son_tree(@{$tree_child{$final_tree1}}){
				if(defined($sample_tree_snp_num{$son_tree}) && $sample_tree_snp_num{$son_tree}>0){
					foreach my $snp1(@{$sample_tree_snp{$son_tree}}){
						push @printout_small,"$snp1\t$tree_sourse{$son_tree}\t+\n";
						$nn1++;
					}
				}
				else{
				}

				if(defined($sample_tree_difsnp_num{$son_tree})){
					foreach my $snp1(@{$sample_tree_difsnp{$son_tree}}){
						push @printout_small,"$snp1\t$tree_sourse{$son_tree}\t-\n";
						$nn1++;
					}
				}
				else{
				}
				push(@inputtrees1,$son_tree);
			}
			push(@printout1,@printout_small);
		}
		else{
			$nn1++;
		}
		$son_nn=$nn1;
	};


	my @inputtrees=$final_tree;
	my @printout="";
	my  $nn=0;
	@printout1="";
	until($nn>0){
		foreach my $tr(@inputtrees){

			get_son_out2($tr);
			$nn=$son_nn+$nn;
			push(@printout,@printout1);
			if($nn==0){
				@inputtrees=@inputtrees1;
			}
			else{}
		}
	}
	if($nn>0){
		print OUT2 join "",@printout;
		print  join "",@printout;
	}
	close OUT2;

	open(OUT,">sample_majortree.txt");
	open(OUT3,">sample_pos_test_result.txt");
	print OUT3 "$sample_N\t$final_sourse\t\n\nposition(hg19)\tmajortree\tSNP\tresult\n";
	open(OUT2,">newfound.xls");
	print OUT2 "tree\ttype\tmj_tree_snp_num\tmatch_snp\tmismatch_snp\tmatch_snp_num\tmismatch_snp_num\t\n";
	print OUT "tree\ttype\tmj_tree_snp_num\tmatch_snp\ttmismatch_snp\tmatch_snp_num(+)\tmismatch_snp_num(-)\tscore\n";
	my $majortree_length=scalar(@major_trees);
	my @major_tree=@major_trees;
	for(my $i=0;$i<$majortree_length;$i++){
		my $key=shift(@major_tree);
		if($key ne ""){

			print OUT "$key\t$tree_sourse{$key}\t";
			my @major_tree_snps=split",",$tree_snp{$key};
			$tree_snp_num{$key}=scalar(@major_tree_snps);
			print OUT "$tree_snp_num{$key}\t";

			if(defined($sample_tree_snp{$key})&& defined($sample_tree_difsnp{$key})){
				print OUT "@{$sample_tree_snp{$key}}\t@{$sample_tree_difsnp{$key}}\t$sample_tree_snp_num{$key}\t-$sample_tree_difsnp_num{$key}\t";
				if($sample_tree_snp_num{$key}>0 && $sample_tree_difsnp_num{$key} > 0){
					print OUT2 "$key\t$tree_sourse{$key}\t$tree_snp_num{$key}\t@{$sample_tree_snp{$key}}\t@{$sample_tree_difsnp{$key}}\t$sample_tree_snp_num{$key}\t$sample_tree_difsnp_num{$key}\t\n";
				}
			}
			elsif(defined($sample_tree_snp{$key})){
				print OUT "@{$sample_tree_snp{$key}}\t-\t$sample_tree_snp_num{$key}\t0\t";
			}
			elsif(defined($sample_tree_difsnp{$key})){
				print OUT "-\t@{$sample_tree_difsnp{$key}}\t0\t-$sample_tree_difsnp_num{$key}\t";
			}
			else{
				print OUT "-\t-\t0\t0\t";
			}
			my $tree_finalscore=0;
			my @mj_trees=split"_",$key;
			my $tree_length=scalar(@mj_trees);
			for(my $j=0;$j<$tree_length;$j++){
				my $small_mjtree=join"_",@mj_trees;
				if(defined($sample_tree_snp{$small_mjtree} && defined($sample_tree_difsnp{$small_mjtree}))){
					$tree_finalscore=$tree_finalscore+$sample_tree_snp_num{$small_mjtree}-$sample_tree_difsnp_num{$small_mjtree};
					}
				elsif(defined($sample_tree_snp{$small_mjtree})){
					$tree_finalscore=$tree_finalscore+$sample_tree_snp_num{$small_mjtree};
					}
				elsif(defined($sample_tree_difsnp{$small_mjtree})){
					$tree_finalscore=$tree_finalscore-$sample_tree_difsnp_num{$small_mjtree};
					}
				else{
					$tree_finalscore=$tree_finalscore;
					}
				splice (@mj_trees, -1);
			}
			print OUT "$tree_finalscore\t\n";
			#######out3#######
			if(defined($sample_tree_snp{$key})){
				foreach my $snp1(@{$sample_tree_snp{$key}}){
					print OUT3 "$snp_pos{$snp1}\t$tree_sourse{$key}\t$snp1\t+\n";
				}
			}
			###########add at 20200526############
			my @key_lists = keys(%sample_tree_difsnp);
			my $fin_res = '';
			my $match_num = 0;
			foreach my $ll(@key_lists){
				my @ll_list = split",", $ll;
				if (grep { $key eq $_ } @ll_list){
					$match_num = 1;
					foreach my$tl(@ll_list){
						if($fin_res eq ''){
							$fin_res = $tree_sourse{$tl};
						}else{
							$fin_res = join(',',$fin_res, $tree_sourse{$tl});
						}
					}
				}
			}
			if($match_num eq 1){
				foreach my $snp1(@{$sample_tree_difsnp{$key}}){
					print OUT3 "$snp_pos{$snp1}\t$fin_res\t$snp1\t-\n";
				}
			}
		}
	}
	close OUT3;
	close OUT;
	close OUT2;
	###########add file 1826list####
	open(IN,"sample_pos_test_result.txt");
	open(OUT,">final_sample_1826list.txt");
	local $/=undef;
	my $scalar = <IN>;
	$scalar=~s/IMS\*JST/IMS-JST/g;
	close IN;
	print OUT "$scalar";
	open(IN,"sample1826list.xls");
	local $/=undef;
	$scalar = <IN>;
	close IN;
	print OUT "$scalar";
	close OUT;
}


sub depfilter{
	my $path=shift;
	my $sample_name=shift;
	open(OUT,">raw_variants_fileter_Y.xls");
	open(OUT1,">sample_ydep_statistic.xls");
	print OUT1 "sample_name\ty_pos\ty_dep>=1\ty_dep>=2\ty_dep>=5\ty_dep>=7\ty_dep>=10\n";
	my $input;
	$input = $sample_path;

	if(-e $input){
		print "input:$input\n";
		open(FI,$input);
		print"$input\n";
		my $depget=0;
		my $depget_Y=0;
		my $chrY=0;
		my $noXY=0;
		my $CHR=0;
		my $chrX=0;
		my $y_1=0;
		my $y_2=0;
		my $y_5=0;
		my $y_7=0;
		my $y_10=0;
		my $chrydep1num=0;
		while(<FI>){
        	chomp;
        	if (/^#/){}
			else{
				$CHR++;
				my @array=split"\t",$_;
				my %sample_vcf;
				my @txt=split";",$array[7];
				foreach my $k(@txt){
					my @tt=split"=",$k;
					$sample_vcf{$tt[0]}=$tt[1];
				}
				my $nn="DP";

				if($array[0] eq "chrY" && $sample_vcf{$nn}>=$dep_Y){
					$depget_Y++;
					print OUT "$_\n";
				}
				elsif($array[0] ne "chrY" && $array[0] ne "chrX" && $sample_vcf{$nn}>=$dep){
					$depget++;
				}
				if($array[0] eq "chrY" && $sample_vcf{$nn}==1){
					$chrydep1num++;
					$dep1_pos{$array[1]}++;
					}
				if($array[0] eq "chrY" && $sample_vcf{$nn}>=1){
					$y_1++;
				}
				if($array[0] eq "chrY" && $sample_vcf{$nn}>=2){
					$y_2++;
					}
				if($array[0] eq "chrY" && $sample_vcf{$nn}>=5){
					$y_5++;
					}
				if($array[0] eq "chrY" && $sample_vcf{$nn}>=7){
					$y_7++;
					}
				if($array[0] eq "chrY" && $sample_vcf{$nn}>=10){
					$y_10++;
					}
				if($array[0] eq "chrY"){
					$chrY++;
					}
				if($array[0] ne "chrY" && $array[0] ne "chrX"){
					$noXY++;
					}
				if($array[0] eq "chrX"){
					$chrX++;
					}
			}
		}
		print OUT1 "$sample_name\t$chrY++\t$y_1\t$y_2\t$y_5\t$y_7\t$y_10\n";
		print"chrY:$chrY\t$dep_Y\t$depget_Y++\nchr1-22:$noXY\t$dep\t$depget\nchrX:$chrX\nchrall:$CHR\n";
		print "sample chrY 1X pos num:$chrydep1num\n";
		print OO "$chrY\t$depget_Y\t$noXY\t$depget\t";
		close FI;
	}
	else{
		print "$input vcf not found!\n";
		}

	close OUT;
	close OUT1;

}



sub gettree{
	my $count=0;
	open(OUT,">major_tree.txt");
	# open(FILE,"/mnt/clustre/users/danni.li/tree/programe/20170510.txt");
	#open(FILE,"$p_path/20170728.txt");
	# open(FILE,"/mnt/ilustre/users/changliang.qiu/danni.li/source/20181027_sourse_tree.txt");
	# open(FILE,"/mnt/ilustre/users/changliang.qiu/danni.li/source/20200127tree.txt");//0505
    # open(FILE,"/mnt/ilustre/users/changliang.qiu/danni.li/source/20200505_tree.txt");//0606
    # open(FILE,"/mnt/ilustre/users/changliang.qiu/danni.li/source/20200606_tree.txt");//0606fix
    # open(FILE,"/mnt/ilustre/users/changliang.qiu/danni.li/source/20200801_tree.txt");//0801
    # open(FILE,"/mnt/ilustre/users/changliang.qiu/danni.li/source/20200801_tree_fix.txt");//20201105
	open(FILE,$tree_list); 
	#############20201105 set for yfull workflow######################
	<FILE>;
	while(<FILE>){
		my $snp_num=0;
		chomp;
		$_ =~ s/\.//g;
		$_ =~ s/\r//g;
		$_=~s/"//g;
		$_=~/(\s+)(.*)/;
		my $space=$1;
		my $str=$2;
		if($str ne ""){
			my @txt=split"\t",$str;
			$str=$txt[0];
			$count = ($space =~ s/\t//g);
			$str=~s/\s//g;
			$str=~s/IMS-JST/IMS*JST/g;
			my @array=split"-",$str;
			my $sourse=$array[0];
			if(defined($array[1])){
				my $snp=$array[1];
				###add at 20200602###
				$snp =~ s/&//g;
				my @snps=split",",$snp;
				$snp_num=scalar(@snps);
				if(exists($hash{$count})){
					my @bear=split"_",$hash{$count};
					my $bear_length=scalar(@bear);
					if($bear_length>=2){
						my @for_tree=split"_",$hash{$count-1};
						my $end=$bear[$count-1]+1;
						splice(@bear,-1);
						my $head=join"_",@bear;
						if($head eq  $hash{$count-1}){
							my $id=join"_",$head,$end;
							$hash{$count}=$id;
							print OUT join"#",$hash{$count},$str;
							}
						else{
							my $end=1;
							my $head=$hash{$count-1};
							my $id=join"_",$head,$end;
							$hash{$count}=$id;
							print OUT join"#",$hash{$count},$str;
							}
					}
					else{
						my $id=$hash{$count}+1;
						$hash{$count}=$id;
						print OUT join"#",$hash{$count},$str;
						}
				}
				else{
					if(exists($hash{$count-1})){
						my $id=join"_",$hash{$count-1},1;
						$hash{$count}=$id;
						print OUT join "#",$hash{$count},$str;
						}
					else{
						my $id=$count;
						$hash{$count}=$id;
						print OUT join"#",$hash{$count},$str;
					}
				}
                my $tree=$hash{$count};
                $tree_sourse{$tree}=$sourse;
				push(@major_trees,$tree);
                $tree_snp{$tree}=$snp;
                foreach my $i(@snps){
					if(exists $snp_tree{$i}){
						###########add at 20200526############
						my $last_tree = $snp_tree{$i};
						my @last_tree_list = split",",$last_tree;
						# if(grep /^$tree$/, @last_tree_list){
						if(grep { $tree eq $_ } @last_tree_list){
							my $new_tree = $last_tree;
							$snp_tree{$i} = $new_tree;
						}
						else{
							my $new_tree = join(",",$last_tree,$tree);
							$snp_tree{$i} = $new_tree;
						}
					}else{ 
						$snp_tree{$i}=$tree;
					}
				}
                print OUT "\n";
			}
	 	}
	}
	close FILE;
	close OUT;
	my @majortreeallsnp=keys(%snp_tree);
	my $majortree_allsnp_num=scalar(@majortreeallsnp);
	print"majortree_allsnp_num:$majortree_allsnp_num\n";
	##############get son nodes######
	foreach my$k(keys(%tree_sourse)){
		my @fa=split"_",$k;
		splice(@fa,-1);
		my $father=join"_",@fa;
		push(@{$tree_child{$father}},$k);
		}
}

sub getlist{
	my %snp_list;
	open(O,">list_pos_mut.txt");
	print O "pos\tSNP\tref\tmut\n";
	#open(IN,"$p_path/Y_list1826.txt");
	# open(IN,"/mnt/ilustre/users/changliang.qiu/danni.li/source/Y_2205_list.txt");//0309
	# open(IN,"/mnt/ilustre/users/changliang.qiu/danni.li/source/Y_20200127_newlist.txt");//0505
    # open(IN,"/mnt/ilustre/users/changliang.qiu/danni.li/source/Y_20200505_list.txt");//0606
    # open(IN,"/mnt/ilustre/users/changliang.qiu/danni.li/source/Y_20200606_list.txt");//0606fix
    # open(IN,"/mnt/ilustre/users/changliang.qiu/danni.li/source/Y_20200801_list.txt");//0801
    # open(IN,"/mnt/ilustre/users/changliang.qiu/danni.li/source/Y_20200801_list_fix.txt");//20201105
	open(IN, $posdb_list) or die "Can' t open file ' $posdb_list'! $!";
	############20201105 set for yfull workflow######################
	<IN>;
	while(<IN>){
		chomp;
		my @array=split"\t",$_;
		$array[1]=~s/\s//g;
		$array[2]=~tr/[a-z]/[A-Z]/;
		$array[1]=~s/IMS-JST/IMS*JST/g;
		$snp_list{$array[1]}=$_;
		if($array[4] eq "*"){
			$array[2]=~tr/[a-z]/[A-Z]/;
			push @all_list_snp,$array[1];
			$pos_snp{$array[0]}=$array[1];
			$snp_pos{$array[1]}=$array[0];
			# push @{$pos_mut{$array[0]}},$array[2];  ##del at 20200602
			push(@{$pos_ref{$array[0]}},$array[2]);   ##add at 20200602
			my@TT=("A","T","C","G");
			foreach my $alt(@TT){
				if($alt eq $array[2]){}
				else{
					# push(@{$pos_ref{$array[0]}},$alt);   ##del at 20200602
					push(@{$pos_mut{$array[0]}},$alt);   ##add at 20200602
					}
			}
		}

		else{
			$array[2]=~tr/[a-z]/[A-Z]/;
			push @all_list_snp,$array[1];
			$pos_snp{$array[0]}=$array[1];
			$snp_pos{$array[1]}=$array[0];
			push(@{$pos_ref{$array[0]}},$array[2]);
			if($array[3] eq "-"){
				my @li=("A","T","C","G");
				foreach my $alt(@li){
					if($alt eq $array[2]){}
					else{
						push@{$pos_mut{$array[0]}},$alt;
						}
					}
				}
			else{
				push(@{$pos_mut{$array[0]}},$array[3]);
				}

			}
			print O "$array[0]\t$pos_snp{$array[0]}\t@{$pos_ref{$array[0]}}\t@{$pos_mut{$array[0]}}\n";
		}
	close IN;
	close O;
	my $all_list_snp_num=scalar(@all_list_snp);
    print"all_list_snp_num:$all_list_snp_num\n";
    my $listsnpnotree=0;
	my $listsnointree=0;
	foreach my $snp(@all_list_snp){
		if(defined($snp_tree{$snp})){
			$listsnp_tree{$snp}=$snp_tree{$snp};
			###########add at 20200526############
			my $str2 = ',';
			if($snp_tree{$snp} =~ /$str2/){
				my @snp_array = split ",",$snp_tree{$snp};
				foreach my$aa(@snp_array){
					push(@{$list_tree_snp{$aa}},$snp);
				}
			}else{
				push(@{$list_tree_snp{$snp_tree{$snp}}},$snp);
			}
			push(@snp_intree,$snp);
			$listsnointree++;
		}
		else{
			push(@snp_notree,$snp);
			$listsnpnotree++;
			}

		}
	print"listsnpnotree:$listsnpnotree\nlistsnpintree:$listsnointree\n";
	#########print notree snp#
	open(OUTT,">notreesnp_list.xls");
	foreach my $ss(@snp_notree){
		if(defined($snp_list{$ss})){
			print OUTT"$snp_list{$ss}\n";
			}
		}
	close OUTT;
}

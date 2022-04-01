#!usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;

# my $list;
my $vcf_file;
my $sample_name;
my $snp_file; #/mnt/ilustre/users/changliang.qiu/danni.li/Programe/Y_sta/chrom_pos_20200512v2.txt
GetOptions(
'vcf=s'=>\$vcf_file,
'sample_name=s'=>\$sample_name,
'snp_file=s'=>\$snp_file,
);

my %pos_infor;
# my $snp_file="/mnt/ilustre/users/changliang.qiu/danni.li/Programe/Y_sta/chrom_pos_20200512v2.txt";
# if(!$list){
# 	print "no in put\n";
# 	exit;
# 	}
# mkdir ("autosome_result");
refget();
# open(IN,"$list");
# while(<IN>){
# 	chomp;
# 	my @array=split"\t",$_;
mkdir ("$sample_name");
chdir("$sample_name");
if(-e "$vcf_file"){
	autosome_snp($vcf_file,$sample_name);
		# system("/mnt/ilustre/users/changliang.qiu/Software/miniconda3/envs/python3/bin/python /mnt/ilustre/users/changliang.qiu/danni.li/Programe/yfull_updata_script/combine_result.py -i *antosome_snp_result.txt -n $array[1]");
		# system("rm *all_antosome_snp_result.txt");
		# system("mv *antosome_snp_result.txt ../autosome_result/");
	}
else{
	print "no $vcf_file";
	}
	# chdir("../");
	# }
# close IN;


sub autosome_snp{
	my $input=shift;
	my $name=shift;
	my $outname=$name."all_antosome_snp_result.txt";
	my %sample_pos_alt;
	open(OUT,">$outname");
	print OUT "chr\tpos\trs\tsample_type\tDP\n";
	open(FILE,"$input");
	while(<FILE>){
		if($_=~/^#/){}
                else{
			chomp;
			my @array=split"\t",$_;
			$array[0]=~/chr(.*)/;
                        my $chrnum=$1;
                        $array[7]=~/(DP=\d+);.*(DP4=\d+,\d+,\d+,\d+)/;
                        my $dp="$1;$2";
                        my @name=split":",$array[8];
                        my @value=split":",$array[9];
                        my %format;
			my @alts=split",",$array[4];
			$array[4]=$alts[0];
                        for(my $i=0;$i<=$#name;$i++){
                                $format{$name[$i]}=$value[$i];
                                }
			# if $array[7].split(";")[0] == 'INDEL':##判断INDEL类型
			if(defined $pos_infor{$chrnum}{$array[1]}){
				#print "$_\n";
				my @inf=split"\t",$pos_infor{$chrnum}{$array[1]};
               		         if($format{"GT"} eq "0/0"){
                        	        $sample_pos_alt{$chrnum}{$array[1]}="$array[3]$array[3]";
                        	        }
                        	elsif($format{"GT"} eq "1/1"){
                        	        $sample_pos_alt{$chrnum}{$array[1]}="$array[4]$array[4]";
                        	        }
                        	elsif($format{"GT"} eq "0/1" ){
                        	        $sample_pos_alt{$chrnum}{$array[1]}="$array[3]$array[4]";
                        	        }
				            elsif($format{"GT"} eq "./." ){
                                        $sample_pos_alt{$chrnum}{$array[1]}="not detected";
                                        }
				            else{
                                        $sample_pos_alt{$chrnum}{$array[1]}="$alts[0]$alts[1]";
                                        } 
				print OUT "$inf[0]\t";
				print OUT "$inf[1]\t$inf[2]\t";
				print OUT "$sample_pos_alt{$chrnum}{$array[1]}\t";
				 print OUT "$dp\n";
			}
			}
		}
	close FILE;
	close OUT;
	}

sub refget{
        open(FILE,$snp_file);
        while(<FILE>){
                chomp;
                my @array=split"\t",$_;
                if($array[2] ne ""){
                        $pos_infor{$array[0]}{$array[1]}="$array[0]\t$array[1]\t$array[2]";
                        }
                }
        close FILE;
	}

	

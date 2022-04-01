#!usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
my %opts;
GetOptions(\%opts,"f=s","o=s","label=s","target=s","n=s","type=s","a=s","b=s","r=s","g=s","method=s");

my $usage = << "USAGE";
	Program : $0
	Discription : Picking Interesting Taxonomic Abundance
	Contact : sen.yang\@majorbio.com
	Usage : perl $0 [options]
		-f	* abundance file
		-type	* tap format  or biom format ,[X/B]default :X	
		-label  ** first column name  of abundance file,like "Taxon" 
			must be selected with type "X"
		-a	alpha diverdity metrics ,[observed_species/margalef/menhinick/dominance/reciprocal_simpson/shannon/equitability/berger_parker_d/mcintosh_d/brillouin_d/strong/fisher_alpha/simpson/mcintosh_e/heip_e/simpson_e/robbins/michaelis_menten_fit/chao1/ACE]default : simpson
			set with type "X" 
		-b	beta diverdity metrics,[braycurtis/canberra/chebyshev/cityblock/correlation/cosine/euclidean/hamming/sqeuclidean/unifrac_unweighted/unifrac_weighted]default : braycurtis
		-target	  target feature list , e.g[TaxonA
						    TaxonB
						    TaxonC
						       ...]
		-r	target's metrics,[rank/abundance] default : rank
		-g	group file , e.g [SampleID Sample1 Sample2 Sample3 ...
				          Group	   GroupA  GroupB  GroupC  ...]
		-method	supervised methods,[distinct/discriminant] default : distinct
        -method2 discriminant
		-o 	* outdir
		-n	number of sample selectd,default: 10
		paramenter with * must be selected
	example:perl $0  -f phylum.xls  -label Taxon -o ./
USAGE

die $usage if (!($opts{f}&&$opts{o}));
$opts{target}=defined $opts{target}?$opts{target}:"F";
$opts{method}=defined $opts{method}?$opts{method}:"distinct";
$opts{method2}=defined $opts{method2}?$opts{method2}:"discriminant";
$opts{type}=defined $opts{type}?$opts{type}:"X";
$opts{n}=defined $opts{n}?$opts{n}:10;
$opts{r}=defined $opts{r}?$opts{r}:"rank";
$opts{a}=defined $opts{a}?$opts{a}:"simpson";
$opts{b}=defined $opts{b}?$opts{b}:"braycurtis";
$opts{g}=defined $opts{g}?$opts{g}:"none";

my $file = $opts{f};
my $label2;
my %hash;
my(@c,@d,@key,@value);

if (! -e "$opts{o}"){
        system("mkdir $opts{o}");
}

### PCL format With Group file###	
if("$opts{g}" ne  "none" and  "$opts{type}" eq "X" ){
	open(A,"$opts{f}" ) || die $!;
	open(B,">tmp.xls" ) || die $!;
	open(C,"$opts{g}" ) || die $!;
	while(<C>){
        	chomp;
        	push @c,$_;
	}
	@key = split/\s+/,$c[0];
        @value = split/\s+/,$c[1];
	for(my $i=1;$i<=$#key;$i++){
	$hash{$key[$i]} = $value[$i];
	}	
	my @key = split/\s+/,$c[0];
	my @value = split/\s+/,$c[1];
	for(my $i=1;$i<=$#key;$i++){
        	$hash{$key[$i]} = $value[$i] ;
	}
	my @a = <A>;	
	shift @a;
	print B @a;
	my @tmp = split/\s+/,$c[1];
	$label2 = $tmp[0];
	system("cat $opts{g} tmp.xls > supervised.shuju.pcl ; rm tmp.xls");
}

### Biom format With Group file### 
if("$opts{o}" ne  "none" and  "$opts{type}" eq "B" ){
	system(" biom convert -i $opts{f} -o b2p.txt  --table-type \"OTU table\" --to-tsv ");
	system(" less b2p.txt |grep -v \"Constructed from biom file\" |sed 's/#//g' > biom2pcl.xls ");
	open(A,"biom2pcl.xls") || die $!;
        open(B,">tmp.xls" ) || die $!;
        open(C,"$opts{g}" ) || die $!;
        while(<C>){
                chomp;
                push @c,$_;
        }
        @key = split/\s+/,$c[0];
        @value = split/\s+/,$c[1];
        for(my $i=1;$i<=$#key;$i++){
	        $hash{$key[$i]} = $value[$i];
        }
        my @key = split/\s+/,$c[0];
        my @value = split/\s+/,$c[1];
        for(my $i=1;$i<=$#key;$i++){
                $hash{$key[$i]} = $value[$i] ;
        }
        my @a = <A>;
        shift @a;
        print B @a;
        my @tmp = split/\s+/,$c[1];
        $label2 = $tmp[0];
        system(" cat $opts{g} tmp.xls > b2p_group.pcl ; rm tmp.xls;rm b2p.txt ");
	
	
}

if  (-e "$opts{g}"){
	open(OUT,">$opts{o}/hash.xls") || die $!;
	print OUT "Sample\tGroup\n";
	foreach (sort keys %hash){
		print OUT  $_."\t".$hash{$_}."\n";
	}
}
open SH,">$opts{o}/PITA.sh";
print SH "

#################### Picking Interesting Taxonomic Abundance ###############################

#### Format PCL ####
if [ \"$opts{type}\" == \"X\" ];then cp $opts{f} shuju.pcl
MicroPITA.py --lastmeta $opts{label} -b $opts{b} -m representative -n $opts{n} shuju.pcl  $opts{o}/representative.xls >>$opts{o}/PITA.log 2>&1
MicroPITA.py --lastmeta $opts{label} -a $opts{a} -m diverse -n $opts{n} shuju.pcl  $opts{o}/diverse.xls >>$opts{o}/PITA.log 2>&1
MicroPITA.py --lastmeta $opts{label} -b $opts{b} -m extreme -n $opts{n}  shuju.pcl  $opts{o}/extreme.xls >>$opts{o}/PITA.log 2>&1
if [ \"$opts{target}\" != \"F\" ];then MicroPITA.py --lastmeta $opts{label} -r $opts{r} -m features -n $opts{n} --targets $opts{target}  shuju.pcl  $opts{o}/feature.xls ; fi >>$opts{o}/PITA.log 2>&1
if [ \"$opts{target}\" != \"F\" ];then cat $opts{o}/diverse.xls $opts{o}/extreme.xls $opts{o}/representative.xls  $opts{o}/feature.xls >$opts{o}/select.xls ; fi >>$opts{o}/PITA.log 2>&1
if [ \"$opts{target}\" == \"F\" ];then cat $opts{o}/diverse.xls $opts{o}/extreme.xls $opts{o}/representative.xls  >$opts{o}/select.xls ; fi   >>$opts{o}/PITA.log 2>&1 
if [ \"$opts{g}\" != \"none\" ];then
        if [ \"$opts{method}\" == \"distinct\" ]; then 
	MicroPITA.py --lastmeta $label2  -n $opts{n} --label $label2  -m distinct supervised.shuju.pcl  $opts{o}/distinct.xls ; fi >>$opts{o}/PITA.log 2>&1
	    if [ \"$opts{method2}\" == \"discriminant\" ]; then 
    MicroPITA.py --lastmeta $label2  -n $opts{n} --label $label2  -m discriminant  supervised.shuju.pcl  $opts{o}/discriminant.xls ;fi ; fi; fi >>$opts{o}/PITA.log 2>&1
#### Format Biom ####
if [ \"$opts{type}\" == \"B\" ];then biom convert -i $opts{f} -o tmp.txt  --table-type \"OTU table\" --to-tsv   >>$opts{o}/PITA.log 2>&1
source /mnt/ilustre/users/linmeng.liu/anaconda.sh ; source activate root; unset PYTHONPATH  >>$opts{o}/PITA.log 2>&1
biom convert -i tmp.txt -o tmp.biom  --table-type \"OTU table\";source deactivate;source /mnt/ilustre/users/meta/profile_meta.sh >>$opts{o}/PITA.log 2>&1
MicroPITA.py -b $opts{b} -n $opts{n} -m representative tmp.biom  $opts{o}/representative.xls >>$opts{o}/PITA.log 2>&1 
MicroPITA.py -a $opts{a} -n $opts{n} -m diverse tmp.biom  $opts{o}/diverse.xls >>$opts{o}/PITA.log 2>&1
MicroPITA.py -b $opts{b} -n $opts{n} -m extreme tmp.biom  $opts{o}/extreme.xls >>$opts{o}/PITA.log 2>&1
if [ \"$opts{target}\" != \"F\" ];then MicroPITA.py -r $opts{r} -m features -n $opts{n} --targets $opts{target}  tmp.biom  $opts{o}/feature.xls ; fi   >>$opts{o}/PITA.log 2>&1
if [ \"$opts{target}\" != \"F\" ];then cat $opts{o}/diverse.xls $opts{o}/extreme.xls $opts{o}/representative.xls  $opts{o}/feature.xls >$opts{o}/select.xls ; fi >>$opts{o}/PITA.log 2>&1
if [ \"$opts{target}\" == \"F\" ];then cat $opts{o}/diverse.xls $opts{o}/extreme.xls $opts{o}/representative.xls  >$opts{o}/select.xls ; fi   >>$opts{o}/PITA.log 2>&1
if [ \"$opts{g}\" != \"none\" ];then
        if [ \"$opts{method}\" == \"distinct\" ]; then 
	MicroPITA.py   --lastmeta   $label2 -n $opts{n}  --label $label2  -m distinct b2p_group.pcl  $opts{o}/distinct.xls >>$opts{o}/PITA.log 2>&1
        else MicroPITA.py  --lastmeta  $label2 -n $opts{n}  --label $label2  -m b2p_group.pcl  $opts{o}/discriminant.xls ;fi ; fi; rm  tmp* ; rm  biom2pcl.xls b2p_group.pcl  ; fi  >>$opts{o}/PITA.log 2>&1
#### Beta Diversity Matrix ####
if [ \"$opts{type}\" == \"X\" ];then 
biom convert -i $opts{f} -o shuju.biom --table-type \"OTU table\"  --to-hdf5 >>$opts{o}/PITA.log 2>&1 
beta_diversity.py -i shuju.biom   -m abund_jaccard,bray_curtis,euclidean,hellinger -o shuju_beta_diversity >>$opts{o}/PITA.log 2>&1
else beta_diversity.py -i $opts{f}   -m abund_jaccard,bray_curtis,euclidean,hellinger -o shuju_beta_diversity ; fi  >>$opts{o}/PITA.log 2>&1

";

if (-e "shuju.pcl"){
	system("mv shuju.pcl $opts{o}/$opts{f}.pcl");
}
system("sh $opts{o}/PITA.sh");
if (-e "shuju.biom"){
	system("rm shuju.biom");
}
chdir "$opts{o}";

###File Processing###
open(IN,"select.xls") ||die $!;
open(OUT,">Picked.xls") ||die $!;
if (-e "distinct.xls" ){
	my @d;
	open(S,"distinct.xls") || die $!;
	open(Y,">G.xls") || die $!;
	print Y "Sample\tGroup\n";
	while(<S>){
        	chomp;
        	@d=split/\s+/,$_;
	}
        shift @d;
      #  pop @d;
        for (sort @d){
        print Y  "$_\t$hash{$_}\n";

	}
}	
if (-e "discriminant.xls" ){
        my @d;
        open(S,"discriminant.xls") || die $!;
        open(Y,">G.xls") || die $!;
        print Y "Sample\tGroup\n";
        while(<S>){
                chomp;
                @d=split/\s+/,$_;
        }
        shift @d;
      #  pop @d;
        for (sort @d){
        print Y  "$_\t$hash{$_}\n";

        }
}
print OUT "Method\t";
my $n=$opts{n};
for (my$i=1;$i<$n;$i++){
	print OUT "Sample".$i."\t";
}
print OUT "Sample".$n."\n";
while(<IN>){
	chomp;
	print OUT $_."\n";
}

system ("dos2unix Picked.xls");
system ("rm select.xls");
chdir "../";
open RCMD,">$opts{o}/pita.cmd.r";
print RCMD "

	#library(dplyr)
	#library(ggplot2)
	list.files(\"./shuju_beta_diversity/\")->files
	files[grep(\"^bray.*.\",files)]->file
	paste0(\"shuju_beta_diversity/\",file)->file
	da <-read.table(file,sep=\"\t\",head=T,check.names=F)
	rownames(da) <-da[,1]
	da <-da[,-1]

#################### PCOA Analysis ##########################################################

	pca <- prcomp(da)
	pc_num =as.numeric(unlist(strsplit(\"1-2\",\"-\")))
	pc_x =pc_num[1]
	pc_y =pc_num[2]
	pc12 <- pca\$x[,c(1,2)]
	pc <-summary(pca)\$importance[2,]*100
	as.data.frame(pc12)->pc12
    pc12\$name<-row.names(pc12)
    write.table(pc12,\"$opts{o}/pcoa.xls\",sep=\"\t\",col.names=T,row.names=F,quote=F)
    ";
#    `~/app/program/R-3.3.3/bin/R --restore --no-save < $opts{o}/pita.cmd.r`;


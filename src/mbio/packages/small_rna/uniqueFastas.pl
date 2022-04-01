#! /usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use List::Util qw ( sum);
use Config::IniFiles;
use Algorithm::Combinatorics qw(combinations);
use Array::Compare;
use Math::Round qw(:all);
# use Venn::Chart;
my %opts;
my $VERSION="2.0";

GetOptions( \%opts, "i=s","u=s","a!","table=s","uniq=s","uniq_all!","list=s", "h!");

my $usage = <<"USAGE";
       Program : $0
       Version : $VERSION
       Contact : quan.guo\@majorbio.com  yuntao.guo\@majorbio.com
       Discription:uniqe given fastq fomart seqs 
       Last updata:2016-9-16
       Usage :perl $0 [options]
                -i*		fasta.ini  config file
[FASTA]
gf1=/data1/wangyan/projets/gaofeng_SR/trim/3960E-SAM-SR/3960E-SAM-SR_trim_adpter.fq.trimmed.single.fa
gf2=/data1/wangyan/projets/gaofeng_SR/trim/YB-SAM-SR/YB-SAM-SR_trim_adpter.fq.trimmed.single.fa


[NAME]
gf1=3960E-SAM-SR
gf2=YB-SAM-SR

[Venn all]
path=/data1/wangyan/projets/gaofeng_SR/all

                
                The first field 'hb1' shows which sample the sequence is from. This prefix is given on
				the command line, and must consist of exactly three alphabetic letters,cannot contain 'x'.
				
                -u		unique.fasta		unique output file name,Defualt:unique.fasta
                -uniq           uniq.fa                 uniq fasta foreach sample
                -uniq_all					output unique file for each sample
                -table		table.xls		unique info table file name ,Defualt:table.xls
                -list		listfile
                -a		outputs progress        
                -h				Display this usage information
                * 				must be given Argument  
                
       example:uniqueFastas.pl 
USAGE

die $usage if ( !(  $opts{i}) || $opts{h} );

$opts{u}= $opts{u}? $opts{u}:"unique.fasta";
$opts{uniq}= $opts{uniq}? $opts{uniq}:"uniq.fasta";
$opts{table}=$opts{table}?$opts{table}:"table.xls";


my $cfg = Config::IniFiles->new(-file => $opts{i});

die "Config file not exits 'FASTA' Section!\n" unless $cfg->SectionExists("FASTA");
my @names=$cfg->Parameters("FASTA");

foreach (@names){
	test_prefix($_);	
}

my %seqs;
my %fa_names;
foreach my $f (@names){
	my $filename=$cfg->val('FASTA',$f);
	warn "Reading Fasta file $filename ...\n" if($opts{a});
	open FASTA, "< $filename" or die "Cannot open file $filename : $! \n";
	my ($fastaname,$fastaseq);
	while (<FASTA>)
	{
	chomp;
	
		if (/^>\s*(\S+)/){
			if($fastaseq && $fastaseq ne ""){
				$fastaseq=~tr/[acgtun\.]/[ACGTTNN]/;
				if(exists($seqs{$fastaseq}{$f})){
					$seqs{$fastaseq}{$f}++;
					push(@{$fa_names{$fastaseq}{$f}},$fastaname);
				}else{
					$seqs{$fastaseq}{$f}=1;
					push(@{$fa_names{$fastaseq}{$f}},$fastaname);
				}
			}
	 		$fastaname=$1;
	 		$fastaseq="";
	    }else{
	   		if ($fastaname ne ""){
				$fastaseq.=$_;
	  		}
		}
	}
	$fastaseq=~tr/[acgtun\.]/[ACGTTNN]/;
	if(exists($seqs{$fastaseq}{$f})){
		$seqs{$fastaseq}{$f}++;
		push(@{$fa_names{$fastaseq}{$f}},$fastaname);
	}else{
		$seqs{$fastaseq}{$f}=1;
		push(@{$fa_names{$fastaseq}{$f}},$fastaname);
	}
	close FASTA;
}

warn "Creating uniq fasta file and unique info table ...\n" if($opts{a});
open TABLE, "> $opts{table}";
	print TABLE "uniq_name";

my %titles;
foreach my $f (@names){
	my $name;
	if($cfg->SectionExists("NAME")&& $cfg->val('NAME',$f)){
		$name=$cfg->val('NAME',$f);
		$titles{$f}=$name;
	}else{
		$name=$f;
		$titles{$f}=$f;
	}
	print TABLE "\t$name";
}

my $filenum =scalar(@names);
	if($filenum>1){
		print TABLE "\ttotal\n";
	}else{
		print TABLE "\n";
	}
	
my $i=0;
my $count=0;
my $commit_interval=1000;
open UNIQ, "> $opts{u}";
open UNIQ2, "> $opts{uniq}";
open SEQ, "> seq.list";
	if($opts{list}){
		open LIST, "> $opts{list}";
	}
	foreach my $s (sort keys(%seqs)){
            $i++;
            print SEQ $s."\n";
		warn("$count entries\n") if (++$count % $commit_interval == 0 && $opts{a});
		my $total=sum(values(%{$seqs{$s}}));
		print TABLE "unq_".$i."_x".$total;
		if($opts{list}){
			print LIST "unq_".$i."_x".$total;
		}
		my $info="";
		foreach my $n (@names){
			if(exists($seqs{$s}{$n})){
				print TABLE "\t".$seqs{$s}{$n};
				$info.="\t".$titles{$n}."=".$seqs{$s}{$n};
                                print UNIQ2 ">".$n."_".$i."_x".$seqs{$s}{$n}."\n".$s."\n";
				if($opts{list}){
					print LIST "\t".$titles{$n}.":".join(",",@{$fa_names{$s}{$n}});
				}
			}else{
				print TABLE "\t0";
			}
		}
		if($opts{list}){
			print LIST "\n";
		}
		if($filenum>1){
			print TABLE "\t$total\n";
		}else{
			print TABLE "\n";
		}
		
		print UNIQ ">unq_".$i."_x".$total.$info."\n";
		print UNIQ $s."\n";
	}
close UNIQ;
close UNIQ2;
	close TABLE;	
	if($opts{list}){
		close LIST;
	}
	# system("cat $opts{table} |autogb -i utf8 -o gb >$opts{table}.1");
	# system("rm $opts{table}");
	# system("mv $opts{table}.1 $opts{table}");
	
if($opts{uniq_all}){
	foreach my $f (@names){
		my $filename=$cfg->val('FASTA',$f);
		if($filename=~/\.(\w+)$/){
			$filename=~s/\.(\w+)$/\.uniq\.$1/;
		}else{
			$filename.= ".uniq";
		}
		warn "Creating uniq fasta file for file ".$cfg->val('FASTA',$f)." ...\n" if($opts{a});
		open SUNIQ, "> $filename";
		$i=0;
		foreach my $s (keys(%seqs)){
			$i++;
			if(exists($seqs{$s}{$f}) && $seqs{$s}{$f} >0){
				print SUNIQ ">".$f."_".$i."_x".$seqs{$s}{$f}."\n";
				print SUNIQ $s."\n";
			}
		}
		close SUNIQ;
	}	
}


my @venn_group=$cfg->GroupMembers('Venn');
if(scalar(@venn_group)>0){



foreach my $g (@venn_group){
	#print "$g\n";
	my $members;
	if($cfg->val($g,'members')){	#jump to else as it's false	
		my $me=$cfg->val($g,'members');
		$me=~s/\s//g;
		my @a=split(",",$me);
		
		my @b;
		foreach my $m (@a){
			if(scalar(grep($m,@names))>0){
				push(@b,$m);
			}else{
				warn "member $m of $g in config file is not exit ,skiping ...\n";
			}
		}
		@b=&uniq(\@b);
		if(scalar(@b)<2){
			warn "members of $g in config file is can't less then 2,skiping ...\n";
			next;
		}else{
			$members=\@b;
		}
	}else{
		$members=\@names;
	}
	my $x=$g;
	$x=~s/^Venn //g;
	my $path=$cfg->val($g,'path')?$cfg->val($g,'path')."/":"./$x";
	mkdir($path,"493") or die "Can't create dir at $path\n" unless -e $path;
        
	# my $vennnum=scalar(@$members);
	# if($vennnum>3||$vennnum==2){
	# 	my $iter= combinations($members,2);
	# 	while(my $c = $iter->next){
	# 		&vennchartR($c,$path);
	# 	}
	# }elsif($vennnum==3){
	# 	my $iter= combinations($members,3);
	# 	while(my $c = $iter->next){
	# 		&vennchartR($c,$path);
	# 	}
	# }
        # =cut
	# &make_sumarytable($members,$path,$g);
}
}

sub vennchartR(){
	my $arry=shift;
	my $path=shift;
		my (%uniq_num,%total_num);
	foreach my $s (keys(%seqs)){
				my @q;
				foreach my $k (sort(keys(%{$seqs{$s}}))){
					foreach my $z (@$arry){
						if($k eq $z){
							push(@q,$z);
						}						
					}					
				}
			if(scalar(@q)>0){
			    my 	$qq=join('&',sort(@q));
			    if(exists($uniq_num{$qq})){
			    	$uniq_num{$qq}++;			    	
			    }else{
			    	$uniq_num{$qq}=1;			    	
			    }
			    foreach my $x (@q){
			    	$total_num{$qq}+=$seqs{$s}{$x};
			    }
			}
			undef @q;
	}
	my (@weight,@weightuniq);
	if(scalar(@$arry)==3){
		
		$weight[0]=0;
		$weight[1]=$total_num{$arry->[0]};
		$weightuniq[0]=0;
		$weightuniq[1]=$uniq_num{$arry->[0]};
		$weight[2]=$total_num{$arry->[1]};
		$weightuniq[2]=$uniq_num{$arry->[1]};
		my @z=($arry->[0],$arry->[1]);
		my $y=join("&",sort @z);
		$weight[3]=$total_num{$y};
		$weightuniq[3]=$uniq_num{$y};
		$weight[4]=$total_num{$arry->[2]};
		$weightuniq[4]=$uniq_num{$arry->[2]};
		@z=($arry->[0],$arry->[2]);
		 $y=join("&",sort @z);
		$weight[5]=$total_num{$y};
		$weightuniq[5]=$uniq_num{$y};
		@z=($arry->[1],$arry->[2]);
		$y=join("&",sort @z);
		$weight[6]=$total_num{$y};
		$weightuniq[6]=$uniq_num{$y};
		$y=join("&",sort @$arry);
		$weight[7]=$total_num{$y};
		$weightuniq[7]=$uniq_num{$y};

	}elsif(scalar(@$arry)==2){
		$weight[0]=0;
		$weight[1]=$total_num{$arry->[0]};
		$weightuniq[0]=0;
		$weightuniq[1]=$uniq_num{$arry->[0]};
		$weight[2]=$total_num{$arry->[1]};
		$weightuniq[2]=$uniq_num{$arry->[1]};
		my @z=($arry->[0],$arry->[1]);
		my $y=join("&",sort @z);
		$weight[3]=$total_num{$y};
		$weightuniq[3]=$uniq_num{$y};
	}
		my (@t,@ss);
		foreach my $x (@$arry){
			#my $xx=$titles{$x};
			#$xx=~s/[^\w]//g;
			push(@ss,$titles{$x});
			push(@t,'"'.$titles{$x}.'"');
		}
		open(R,"> $path/".join("_vs_",@ss).".R");
		print R "library(Vennerable)\n";
		print R "pdf(\"$path/".join("_vs_",@ss).".uniq.pdf\")\n";
		print R "myVenn<-Venn(SetNames=c(".join(",",@t)."),Weight=c(".join(",",@weightuniq)."))\n";
		print R "plot(myVenn,doWeights=FALSE,show=list(Faces=TRUE))\n";
		print R "dev.off()\n";
		print R "pdf(\"$path/".join("_vs_",@ss).".total.pdf\")\n";
		print R "myVenn<-Venn(SetNames=c(".join(",",@t)."),Weight=c(".join(",",@weight)."))\n";
		print R "plot(myVenn,doWeights=FALSE,show=list(Faces=TRUE))\n";
		print R "dev.off()\n";
		close R;
		system("R --no-save < $path/".join("_vs_",@ss).".R");
}

sub make_sumarytable(){
	
	my $mem=shift;
	my $path=shift;
	my $grouname=shift;
	print 
	warn "Creating sumary table for $grouname ...\n" if($opts{a});
	my $table_file=$cfg->val($grouname,'sumary')?$path.$cfg->val($grouname,'sumary'):$path."sumary.xls";
	open SUMARY, "> $table_file";
	print SUMARY "type\tUnique_sRNA\tpercent_of_Unique\ttotal_sRNA\tpercent_of_Total\n";
	my $num=scalar(@$mem);
	my $uniq=0;
	my $total=0;
	my %uniqseqs;
	foreach my $s (keys(%seqs)){
		foreach my $m (@$mem){
			if(exists($seqs{$s}{$m})){
				$uniqseqs{$s}=1;
				$total+=$seqs{$s}{$m};
			}
		}
	}
	$uniq=scalar(keys(%uniqseqs));
	my @ss;
	foreach my $x(@$mem){
		push (@ss,$titles{$x});
	}
	print SUMARY "Total_sRNAs(".join(",",@ss).")\t$uniq\t100%\t$total\t100%\n";
	my (%uniq_num,%total_num);
	foreach my $s (keys(%seqs)){
				my @q;
				foreach my $k (sort(keys(%{$seqs{$s}}))){
					foreach my $z (@$mem){
						if($k eq $z){
							push(@q,$z);
						}						
					}					
				}
			if(scalar(@q)>0){
			    my 	$qq=join('&',sort(@q));
			    #print $qq."\n";
			    if(exists($uniq_num{$qq})){
			    	$uniq_num{$qq}++;			    	
			    }else{
			    	$uniq_num{$qq}=1;			    	
			    }
			    foreach my $x (@q){
			    	$total_num{$qq}+=$seqs{$s}{$x};
			    }
				#if($comp->perm(\@a,\@b)){
#				if($aa eq $bb){
#					$uniq_num++;
#					$total_num+=sum(values(%{$seqs{$s}}));
#				}
			}
			undef @q;
	}
	#my $comp=Array::Compare->new();
	foreach(my $i=1;$i<=$num;$i++){		
		my $iter= combinations($mem,$i);
		while(my $c = $iter->next){
			#print join(",",@$c)."\n";
			#print SUMARY "Total_sRNAs\t$uniq_num\t100%\t$total_num\t100%\n";
			#my $diff=Array::Diff->diff(\@old,\@new);	
			my @a=sort(@$c);	
			my $aa=join('&',@a);
			warn "count for ".join("&",@a)." ...\n" if($opts{a});
			$uniq_num{$aa} =0 unless exists($uniq_num{$aa});
			$total_num{$aa} =0 unless exists($total_num{$aa});
			my $percent1=nearest(.01, $uniq_num{$aa}/$uniq * 100);
			my $percent2=nearest(.01, $total_num{$aa}/$total * 100);
			my @aaa;
			foreach my $x (@a){
				push(@aaa,$titles{$x});
			}
			print SUMARY join("&",@aaa)." specific\t$uniq_num{$aa}\t".$percent1."%\t$total_num{$aa}\t".$percent2."%\n";
		}
		
	}
	close SUMARY;
	system("cat $table_file |autogb -i utf8 -o gb >$table_file.1");
	system("rm $table_file");
	system("mv $table_file.1 $table_file");
}

sub uniq {
	my $array      = shift;
	my %hash       = map { $_ => 1 } @$array;
	my @uniq_array = sort( keys %hash );
	return @uniq_array;
}

sub test_prefix{

    my $prefix=shift;

    unless($prefix=~/^\w\w\w$/ and !($prefix=~/_/)){

        die "prefix $prefix does not contain exactly three alphabet letters\n";
    }
    if($prefix=~/x/){
    	die "prefix $prefix contain letter 'x'\n";
    }
    return;
}



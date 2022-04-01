#!/usr/bin/perl -w
use warnings;
use strict;

die "usage: perl $0 more1000.list less1000.list min_length newbler_stat out_dir\n" unless(@ARGV == 5);

my ($list1,$list2,$min_len,$newbler_stat,$out)=@ARGV;
my ($combine_file,%more_list,%less_list,%stat_file);
#my $ref_count;
open LIST1,$list1 or die "can't open file $list1 $!\n";
while(<LIST1>){
	chomp;
	my @tmp = split /\t/;
	$more_list{$tmp[1]} = $tmp[0];
}
close LIST1;
open LIST2,$list2 or die "can't open file $list2 $!\n";
while(<LIST2>){
	chomp;
	#my $file = (split /\//)[-1];
	#my $sam = (split /\./,$file)[0];
	#$sam =~ s/newbler_input_//;
	#$less_list{$sam} = $_;
	my @tmp = split /\t/;
	$less_list{$tmp[1]} = $tmp[0];
}
close LIST2;

open STAT,"$newbler_stat" or die "can't open file $newbler_stat $!\n";
while(<STAT>){
	if(/Outlier/ or /Singleton/){
		chomp;
		my $contig = (split(/\t/,$_))[0];
		my $sam = (split(/_contig/,$contig))[0];
        if(!/_contig/){
            my @tmp  = (split(/_/,$contig));
            pop(@tmp);
            pop(@tmp);
            $sam = join('_', @tmp);
        }
        if(/Megahit_Mix/){
            $sam = 'Megahit_Mix';
        }
		#$stat_file{$sam} .= $contig.",";
		$stat_file{$sam}{$contig} = 1;
	}
}
close STAT;

open FA,">$out/result.tmp.fa" or die $!;
foreach my $sample (sort keys %stat_file){
	#`mkdir -p $out/$sample`;
	open IN,">$out/$sample.contig.fa" or die "Can't write file $out/$sample.contig.fa $!\n";
	open LFA,$more_list{$sample} or die "Can't find file $more_list{$sample} $!\n";
	while(<LFA>){
		print IN "$_";
		if(/^>/){
			my $line = $_;
			$line =~ s/>//;
			print FA ">$sample"."_"."$line";
		}else{
			print FA "$_";
		}
	}
	close LFA;
	open SFA,$less_list{$sample} or die "Can't find file $less_list{$sample} $!\n";
	#print "read $sample less1000 fasta\n";
	my $boolean = 0;
	while(<SFA>){
		if(/^>/){
			chomp;
			my $line = $_;
			$line =~ s/>//;
			my @tmp = split / /,$line;
			#print "line is $line\n";
			#print "tmp1 is $tmp[0],tmp2 is $tmp[1],tmp3 is $tmp[2]\n";
			$tmp[1] =~ s/length_//;
            if (defined($tmp[3]) and $tmp[3] =~ /len=/){
                $tmp[1] = $tmp[3];  # megahit's length
                $tmp[1] =~ s/len=//;
            }
			my $id = $tmp[0];
			#my $string = $id.",";
			#if($stat_file{$sample} =~ /$string/){
			my $contig = "$sample"."_"."$id";
			if($stat_file{$sample}{$contig}){
				if($tmp[1] >= $min_len){
					$boolean = 1;
					print IN "$_\n";
					print FA ">$sample"."_"."$line\n";
				}else{
					$boolean = 0;
					#print "boolean is 0,tmp1 is $tmp[1],minlen is $min_len\n";
				}
			}
		}else{
			if($boolean == 1){
				print IN "$_";
				print FA "$_";
			}
		}
	}
    close SFA;
    close IN;
}
my $newbler_fasta = $newbler_stat;
$newbler_fasta =~ s/454ReadStatus.txt/newbler.contig.fa/;
open NEWB, $newbler_fasta or die "Can't find file $newbler_fasta\n";
my $seq = "";
while(<NEWB>){
	chomp;
	if(/>/){
		if ($seq ne ""){
			print FA "$seq\n";
			$seq = "";
		}
		my $line = $_;
		$line =~ s/>//;
		print FA ">newbler_$line\n";
	}else{
		$seq .= $_;
	}
}
print FA "$seq\n";
close NEWB;
close FA;
`mv $out/result.tmp.fa $out/Newbler_Mix.contig.fa`;
#`cat $out/result.tmp.fa $newbler_fasta > $out/Newbler_Mix.contig.fa`;
#`assembly.statistics.pl $rout/RESULT.fa $rout/RESULT.fa.stat`;
#`rm -rf $out/result.tmp.fa`;

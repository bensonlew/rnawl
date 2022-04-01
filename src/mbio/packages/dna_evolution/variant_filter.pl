#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($vcffile,$outdir,$type,$config,$gfile);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"vcffile:s"=>\$vcffile,
	"outdir:s"=>\$outdir,
	"config:s"=>\$config,
	"group:s"=>\$gfile,
			) or &USAGE;
&USAGE unless ($vcffile and $outdir and $config and $gfile);
mkdir $outdir if (!-d $outdir);
open In,$gfile;
my %groupid;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/);
	my ($id,$groupid)=split(/\s+/,$_);
	$groupid{$id}=$groupid;
}
close In;
open In,$config;
my %varianttype;
my %variantann;
my %varianteff;
my %sampledep;
my %samplehomo;
my %samplediff;
my %groupad;
my %groupmiss;
my %groupmaf;
my %position;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/||/^#/);
	if (/Variant Type/) {
		my @info=split(",",(split(/\=/,$_))[-1]);
		foreach my $info (@info) {
			$varianttype{$info}=1;
		}
	}elsif (/Variant Eff/) {
		my @info=split(",",(split(/\=/,$_))[-1]);
		foreach my $info (@info) {
			$varianteff{$info}=1;
		}
	}elsif (/Variant Ann/) {
		my @info=split(",",(split(/\=/,$_))[-1]);
		foreach my $info (@info) {
			$variantann{$info}=1;
		}
	}elsif (/Sample Diff/) {
		#sample1 dep1 dep2 hete sample2 dep3 dep4 hete diff
		my ($sample1,$dep1,$dep2,$hete1,$sample2,$dep3,$dep4,$hete2,$diff)=split(",",(split(/\=/,$_))[-1]);
		$sampledep{$sample1}{min}=$dep1;
		$sampledep{$sample1}{max}=$dep2;
		$samplehomo{$sample1}=$hete1;
		$sampledep{$sample2}{min}=$dep1;
		$sampledep{$sample2}{max}=$dep2;
		$samplehomo{$sample2}=$hete2;
		$samplediff{$sample1}{$sample2}=$diff;
		$samplediff{$sample2}{$sample1}=$diff;
	}elsif (/Group Info/) {
		#group1 ad1 ad2 miss1 miss2 maf1 maf2
		my ($group,$ad1,$ad2,$miss1,$miss2,$maf1,$maf2)=split(",",(split(/\=/,$_))[-1]);
		$groupad{$group}{min}=$ad1;
		$groupad{$group}{max}=$ad2;
		$groupmiss{$group}{min}=$miss1;
		$groupmiss{$group}{max}=$miss2;
		$groupmaf{$group}{min}=$maf1;
		$groupmaf{$group}{max}=$maf2;
	}elsif (/Region/) {
		my ($chr,$start,$end)=split(",",(split(/\=/,$_))[-1]);
		$position{region}=join(",",$chr,$start,$end);
	}
}
close In;

open In,$vcffile;
open VCF,">$outdir/pop.filtered.vcf";
open Out,">$outdir/pop.table";
my @sample;
my %stat;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/ );
	if (/^##/) {
		print VCF "$_\n";
	}elsif (/^#/) {
		(undef,undef,undef,undef,undef,undef,undef,undef,undef,@sample)=split(/\t/,$_);
		print VCF "$_\n";
		my @out;
		foreach my $sample (sort keys %sampledep) {
			push @out,$sample."_Genotype";
			push @out,$sample."_Allele_Depth";
		}
		foreach my $group (sort keys %groupad) {
			push @out,$group."_Averge_Depth";
			push @out,$group."_miss_Ratio";
			push @out,$group."_Genotype_Frequency";
		}
		print Out "#Chr\tPos\tType\tRef\tAlt\tAnnotion\t",join("\t",@out),"\n";
	}else{
		my ($chr,$pos,$id,$ref,$alt,$qual,$filter,$info,$format,@info)=split(/\t/,$_);
		my $next=0;
		my @alles=split(",",join(",",$ref,$alt));
		my %vlen;
		foreach my $alle (@alles) {
			$vlen{length($alle)}=1;
		}
		my $vtype;
		if (scalar keys %vlen >1) {
			$vtype="INDEL";
		}else{
			$vtype="SNP";
		}
		$next=1 if (!exists $varianttype{$vtype});
		if ($next ==1) {
			print $vtype;die;
		}
		next if ($next == 1);
		my %ann;
		my %eff;
		my $anntest=0;
		my $efftest=0;
		my @anns;
		if($info=~/ANN=([^\;]*)/g){
			my @ann=split(/\,/,$1);
			for (my $i=0;$i<@ann;$i++) {
				my @str=split(/\|/,$ann[$i]);
				$str[0]||="--";
				$str[1]||="--";
				$str[2]||="--";
				$str[3]||="--";
				$str[4]||="--";
				my $ann=join("|",$str[0],$str[1],$str[2],$str[3],$str[4]);
				$ann{$str[2]}++;
				if (exists $variantann{$str[1]}) {
					$anntest=1;
				}
				$eff{$str[1]}++;
				if (exists $varianteff{$str[2]}) {
					$efftest=1;
				}
				push @anns,$ann;
			}
		}
		if ($anntest == 0 || $efftest ==0) {
			$next=1;
		}
		next if ($next == 1);
		if ($next ==1) {
			print Dumper %variantann;
			print Dumper %varianteff;
			print Dumper %ann;
			print Dumper %eff;
			die;
		}
		my @format=split(/\:/,$format);
		my %statsample;
		my %statgroup;
		my %allegroup;
		my %outgroup;
		$statsample{ref}{gt}="$ref/$ref";
		$statsample{ref}{dep}=10;
		for (my $i=0;$i<@info;$i++) {
			next if (!exists $groupid{$sample[$i]});
			my $gid=$groupid{$sample[$i]};
			#	print $gid;die;
			my @ginfo=split(/\:/,$info[$i]);
			for (my $j=0;$j<@ginfo;$j++) {
				if ($format[$j] eq "GT") {
					if ($ginfo[$j] eq "./.") {
						$statgroup{$gid}{miss}++;
						$statsample{$sample[$i]}{gt}="N/N";
					}else{
						my @alle=split(/\//,$ginfo[$j]);
						my $gt=join("/",sort($alles[$alle[0]],$alles[$alle[1]]));
						$statsample{$sample[$i]}{gt}=$gt;
						$outgroup{$gid}{$gt}++;
						$allegroup{$gid}{$alles[$alle[0]]}++;
						$allegroup{$gid}{$alles[$alle[1]]}++;
					}
					$statgroup{$gid}{num}++;
				}
				if ($format[$j] eq "DP") {
					$ginfo[$j]="0" if ($ginfo[$j] eq ".");
					$statsample{$sample[$i]}{dp}=$ginfo[$j];
					$statgroup{$gid}{dp}+=$ginfo[$j];
				}
				if ($format[$j] eq "AD") {
					$ginfo[$j]="0" if ($ginfo[$j] eq ".");
					$statsample{$sample[$i]}{ad}=$ginfo[$j];
				}
			}
		}
		next if ($next == 1);
		foreach my $sample (sort keys %sampledep) {
			$statsample{$sample}{dp}||=0;
			if ($statsample{$sample}{gt} eq "N/N") {
				$next=1;
			}
			if ($statsample{$sample}{dp} < $sampledep{$sample}{min} || $statsample{$sample}{dp} > $sampledep{$sample}{max}) {
				$next=1;
			}
			my @gt=split(/\//,$statsample{$sample}{gt});
			if ($gt[0] eq $gt[1] && $samplehomo{$sample} == 1) {
				$next=1;
			}elsif ($gt[0] ne $gt[1] && $samplehomo{$sample} != 1) {
				$next=1;
			}
		}
		next if ($next == 1);
		foreach my $s1 (sort keys %samplediff) {
			foreach my $s2 (sort keys %{$samplediff{$s1}}) {
				next if (!exists $statsample{$s1}{gt} || !exists $statsample{$s2}{gt});
				if ($statsample{$s1}{gt} eq $statsample{$s2}{gt} && $samplediff{$s1}{$s2}==1) {
					$next=1;
				}elsif($statsample{$s1}{gt} ne $statsample{$s2}{gt} && $samplediff{$s1}{$s2}!=1){
					$next=1;
				}
			}
		}
		next if ($next == 1);
		foreach my $gid (sort keys %groupad) {
			$statgroup{$gid}{dp}||=0;
			$statgroup{$gid}{miss}||=0;
			$allegroup{$gid}{$ref}||=0;
			if ($statgroup{$gid}{dp}/$statgroup{$gid}{num} < $groupad{$gid}{min} || $statgroup{$gid}{dp}/$statgroup{$gid}{num} > $groupad{$gid}{max}) {
				$next=1;
			}
			if ($statgroup{$gid}{miss}/$statgroup{$gid}{num} < $groupmiss{$gid}{min} || $statgroup{$gid}{miss}/$statgroup{$gid}{num} > $groupmiss{$gid}{max}) {
				$next=1;
			}
			if ($allegroup{$gid}{$ref}/($statgroup{$gid}{num}*2) < $groupmaf{$gid}{min} || $allegroup{$gid}{$ref}/($statgroup{$gid}{num}*2) > $groupmaf{$gid}{max}) {
				$next=1;
			}
		}
		next if ($next == 1);
		foreach my $region (values %position) {
			my ($chrs,$start,$end)=split(",",$position{region});
			if ($chrs ne $chr){
				$next=1;
			}
			if ($pos <$start || $pos > $end){
				$next=1;
			};
		}
		next if ($next == 1);
		foreach my $ann (sort keys %eff) {
			$stat{eff}{$ann}{$vtype}++;
		}
		foreach my $ann (sort keys %ann) {
			$stat{ann}{$ann}{$vtype}++;
		}
		
		print VCF "$_\n";
		my @out;
		foreach my $sample (sort keys %sampledep) {
			push @out,$statsample{$sample}{gt};
			push @out,$statsample{$sample}{ad};
		}
		foreach my $gid (sort keys %groupad) {
			push @out,$statgroup{$gid}{dp}/$statgroup{$gid}{num};
			push @out,$statgroup{$gid}{miss}/$statgroup{$gid}{num};
			my @gts=sort {$outgroup{$gid}{$a}<=>$outgroup{$gid}{$b}} keys %{$outgroup{$gid}};
			my @values=sort {$a<=>$b} values %{$outgroup{$gid}};
			push @out,join(":",join(",",@gts),join(",",@values));
		}

		print Out "$chr\t$pos\t$vtype\t$ref\t$alt\t",join(";",@anns),"\t",join("\t",@out),"\n";
	}
}
close In;
open Out,">$outdir/eff.type";
foreach my $ann (sort keys %{$stat{eff}}) {
	$stat{eff}{$ann}{SNP}||=0;
	$stat{eff}{$ann}{INDEL}||=0;
	print Out join("\t",$ann,$stat{eff}{$ann}{SNP},$stat{eff}{$ann}{INDEL},$stat{eff}{$ann}{SNP}+$stat{eff}{$ann}{INDEL}),"\n";
}
close Out;
open Out,">$outdir/function.type";
foreach my $ann (sort keys %{$stat{ann}}) {
	$stat{ann}{$ann}{SNP}||=0;
	$stat{ann}{$ann}{INDEL}||=0;
	print Out join("\t",$ann,$stat{ann}{$ann}{SNP},$stat{ann}{$ann}{INDEL},$stat{ann}{$ann}{SNP}+$stat{ann}{$ann}{INDEL}),"\n";
}
close Out;

#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
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
Contact:        minghao.zhang\@majorbio.com;
Script:			$Script
Description:
	fq thanslate to fa format
	eg:
	perl $Script -i -o -k -c

Usage:
  Options:
	--vcf	<file> vcf file to filer
	--config	<file>	filter config
	--out	<dir>	output dir
	--group	<file>	group list file

USAGE
        print $usage;
        exit;
}

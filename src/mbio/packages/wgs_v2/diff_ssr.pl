#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($input,$output,$bam);
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"input:s"=>\$input,
	"bam:s"=>\$bam,
	"output:s"=>\$output,
			) or &USAGE;
&USAGE unless ($input and $output);
open In,$input;
my %ssr;
open Result,">$output.ssr.result";
print Result "#chr\tpos\tSSRbase\tUnit\tDep\tTotaldep\n";
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/ || /^#/ ||/ID/);
	my ($chr,$ssrn,$ssrtype,$ssr,$size,$start1,$end1)=split(/\s+/,$_);
	next if ($ssrtype =~ "c");
	my $pos=join("\t",$chr,$start1,$end1);
	my $ssrseq="";
	my $oldpos;
	my $ssrbase;
	while ($ssr =~ m/\((\D+)\)(\d+)/g) {
		$ssrbase=$1;
		if (length($ssrseq)==0) {
			$ssrseq.=$1 x $2;
			$oldpos=pos($ssr);
		}else{
			my $pos=pos($ssr)-$oldpos-length($1.$2)-2;
			$ssrseq.=substr($ssr,$oldpos,$pos).$1x$2;
			$oldpos=pos($ssr);
		}
	}
	$ssrseq=uc($ssrseq);
	open Out,">$output.tmp.bed";
	print Out $pos;
	close Out;
	#open Bamread,"samtools view -M $bam -L $output.tmp.bed|";
	open Bamread,"samtools view $bam -L $output.tmp.bed|";
	my %stat;
	my $read;
	my @alignInfo;
	my %reads;
	while ($read=<Bamread>) {
		chomp;
		next if ($read eq ""||$read =~ /^$/);
		my ($seqID,$len,$chr,$start,undef,undef,undef,undef,undef,$seq,undef)=split(/\s+/,$read);
		next if (150-$start1+$start < length($ssrseq));
		my $i=0;
		my $minUnitRep=2;
		my $maxUnitRep=int(length($seq)/length($ssrbase));
		next if ($minUnitRep > $maxUnitRep);
		while ($seq =~ m/(($ssrbase){$minUnitRep,$maxUnitRep})/g){
			my $ssrSeqInRead = $1; 
			my $ssrSeqStart=index($seq,$ssrSeqInRead,$i);
			my $ssrSeqStop=$ssrSeqStart+length($ssrSeqInRead);
			if (($ssrSeqStart > 0) and ($ssrSeqStop < length($seq))) {
				my $length5p=$ssrSeqStart;
				my $overhang5p=substr($seq,0,$length5p);
				my $length3p=length($seq)-$ssrSeqStop;
				my $overhang3p=substr($seq,-$length3p);
				my $ssrUnitNum;
				if (length($overhang5p) < length ($ssr)) {
					if ($overhang5p eq substr($ssr,-length($overhang5p))){
						#print $seqID." was discarded due to undefined 5' limit"."\n";
						$i=$ssrSeqStop;
						next;
					}	
					elsif (length($overhang3p) < length ($ssr)){
						if ($overhang3p eq substr($ssr,0,length($overhang3p))) {
						#	print $seqID." was discarded due to undefined 3' limit"."\n";
							$i=$ssrSeqStop;
							next;
						}
						else {
							$ssrUnitNum = int(length($ssrSeqInRead)/length($ssrbase));
							@alignInfo=($seq,$ssrUnitNum);
							$reads{$ssrSeqStart+$start}{$ssrUnitNum}++;
						}
					}
				}
				else {
						$ssrUnitNum = int(length($ssrSeqInRead)/length($ssrbase));
						@alignInfo=($seq,$ssrUnitNum);
						$reads{$ssrSeqStart+$start}{$ssrUnitNum}++;
				}
			}
			else { 
			#	print $seqID.' was discarded because the SSR is to one of two extremes '.$seq."\n";
			}
		}
	}
	close Bamread;
	foreach my $pos (sort keys %reads) {
		my @outd;
		my @outu;
		my $dep;
		foreach my $unit (sort keys %{$reads{$pos}}) {
			push @outd,$unit;
			push @outu,$reads{$pos}{$unit};
			$dep+=$reads{$pos}{$unit};
		}
		print Result join("\t",$chr,$pos,$ssrbase,join(":",@outd),join(":",@outu),$dep),"\n";
	}
}
close In;
close Result;
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
Contact:        long.huang\@majorbio.com;
Script:			$Script
Description:
Usage:
  Options:
	"input:s"=>\$input,
	"bam:s"=>\$bam,
	"output:s"=>\$output,
  -h         Help

USAGE
        print $usage;
        exit;
}

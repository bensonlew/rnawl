#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($in,$out);        #1
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"in:s"=>\$in,     #2
	"out:s"=>\$out,   #3
			) or &USAGE;
&USAGE unless ($in and $out); #4
#######################################################################################
$in = ABSOLUTE_DIR($in);
my %times;
my @files = glob "$in/*.psmc.result";
foreach my $file(@files){
	open IN,$file;
	while (<IN>){
		chomp;
		next if ($_ =~/^\"/||/^$/);
		my ($generation,undef)= split /\s+/;
		$times{$generation}++;
	}
	close IN;
}
my @outputnames; 
my %output;
my %count; #save the first non-null Ne
my $number=-1; #record pop order
foreach my $file(@files){
	open IN,$file;
	my %temp;
	push @outputnames,basename($file,".psmc.result");
	$number ++;
	while (<IN>){
		chomp;
		next if ($_ =~/^\"/||/^$/);
		my ($generation,$Ne)= split /\s+/;
		$temp{$generation} = $Ne;
	}
	close IN;
	my $lastNe='--';
	my $flag=0;
	foreach (sort {$a <=> $b} keys %times){
		if (exists $temp{$_}){
			$lastNe = $temp{$_};
			$output{$_} .= "$temp{$_}\t";
			if ($flag==0){
				$count{"$number"}= $temp{$_};
			}
			$flag++;
		}else{
			$output{$_} .= "$lastNe\t";
		}
	}
}
open OUT,">$out";
my $temp = join "\t",@outputnames;
print OUT "Generation\t$temp\n";
foreach (sort {$a <=> $b} keys %output){
	my @full = split /\t/,$output{$_};
	foreach my $num(sort {$a <=> $b} keys %count){
		$full[$num]=~s/\-\-/$count{$num}/;
	}
	my $line = join "\t",@full;
	print OUT "$_\t$line\n";
}
close OUT;
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

sub USAGE {           #5
        my $usage=<<"USAGE";
Contact:	tong.wang\@majorbio.com
Version:	$version
Script:		$Script
Description:	Merge x.psmc.result into a table
Usage:
  Options:
  -in	<dir>	input dirname of x.psmc.result 
  -out	<file>	output file; [absolute path]
  -h		Help

USAGE
        print $usage;
        exit;
}

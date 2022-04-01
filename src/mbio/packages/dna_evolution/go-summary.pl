#!/usr/bin/perl -w
use strict;
use warnings;
my $BEGIN_TIME=time();
use Getopt::Long;
my ($input,$output,$annotation,$level);
use Data::Dumper;
use GO::OntologyProvider::OboParser;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $version="1.0.0";
GetOptions(
	"help|?" =>\&USAGE,
	"obofile:s"=>\$input,
	"annotation:s"=>\$annotation,
	"output:s"=>\$output,
	"level:s"=>\$level,
			) or &USAGE;
&USAGE unless ($input and $output);
$level||=2;
my $process   = GO::OntologyProvider::OboParser->new(ontologyFile => $input,
                                                     aspect       => 'P');
my $component = GO::OntologyProvider::OboParser->new(ontologyFile => $input,
                                                     aspect       => 'C');
my $function  = GO::OntologyProvider::OboParser->new(ontologyFile => $input,
                                                     aspect       => 'F');
my @pnodes = $process->allNodes;
my @cnodes = $component->allNodes;
my @fnodes = $function->allNodes;
my %terms;
foreach my $nodes (@pnodes) {
	$terms{$nodes->goid}="biological_process";
}
foreach my $nodes (@cnodes) {
	$terms{$nodes->goid}="cellular_component";
}
foreach my $nodes (@fnodes) {
	$terms{$nodes->goid}="molecular_function";
}
$/="\n";
#print Dumper  \%terms;die;
open In,$annotation;
my %stat;
my %enrich;
my %idterm;
my $totalgene;
my %lenrich;
my %goid;
while (<In>) {
	chomp;
	next if ($_ eq ""||/^$/||/^#/);
	my ($geneinfo,$geneid,$tran_id,$biotype,$chr,$pos1,$pos2,$high,$moderate,$low,$modifier,$NR,$NRD,$Uni,$UniD,$KEGG,$KEGGD,$goID,$term,undef)=split(/\t/,$_);
	my $gene;
	if ($geneinfo ne "--") {
		$gene = $geneinfo;
	}else{
		$gene = $geneid;
	}
	$totalgene++;
	next if ($goID eq "--");
	my @goid=split(/\,/,$goID);
	foreach my $goid (@goid) {
		my $node;
		$enrich{$goid}{total}++;
		next if (!exists $terms{$goid} || $goid eq "GO:0003674" || $goid eq "GO:0005575" ||$goid eq "GO:0008150");
		if ($terms{$goid} eq "biological_process") {
			$node= $process->nodeFromId($goid);
		}elsif ($terms{$goid} eq "molecular_function") {
			$node= $function->nodeFromId($goid);
		}elsif ($terms{$goid} eq "cellular_component") {
			$node= $component->nodeFromId($goid);
		}
		foreach my $ancestor ($node->ancestors){
			next if ($ancestor->meanLengthOfPathsToRoot != $level+1);
			my $goid1=$ancestor->goid;
			$stat{$goid1}{$gene}=1;
			$goid{$goid1}{$goid}=1;
			$idterm{$goid1}=$ancestor->term;
			$lenrich{$goid1}{total}++;
		}	
	}
}
close In;
#print Dumper \%goid;die;
open Out,">$output.detail";
print Out "goid\tgoterm\tis_a\tgenes\n";
foreach my $goid (sort {$terms{$a} cmp $terms{$b}} keys %stat) {
	print Out $goid,"\t",$idterm{$goid},"\t",$terms{$goid},"\t",join(",",keys %{$stat{$goid}}),"\n";
}
close Out;
open Out,">$output.$level.enrich";
foreach my $goid(sort {$terms{$a} cmp $terms{$b}} keys %stat) {
	$lenrich{$goid}{total}||=0;
	my $key = scalar keys %{$goid{$goid}}; 
	print Out $goid,"\t",$idterm{$goid},"\t",$terms{$goid},,"\t",$lenrich{$goid}{total},"\t",$totalgene,"\t",$key,"\t",join(",",keys %{$goid{$goid}}),"\n";
}
close Out;
open Out,">$output.entrich";
foreach my $goid(sort {$terms{$a} cmp $terms{$b}} keys %stat) {
	$enrich{$goid}{total}||=0;
	print Out $goid,"\t",$idterm{$goid},"\t",$terms{$goid},"\t",$enrich{$goid}{total},"\t",$totalgene,"\n";
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
Contact:        long.huang\@majorbio.com;
Script:			$Script
Description:
Usage:
  Options:
	"obofile:s"=>\$input,  #/mnt/ilustre/users/long.huang/DataBase/2018-8-27/GO/go-basic.obo
	"annotation:s"=>\$annotation,   pop.summary
	"output:s"=>\$output,
	"level:s"=>\$level,   default 2
  -h         Help

USAGE
        print $usage;
        exit;
}

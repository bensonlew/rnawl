#! /usr/bin/perl
#guanqing.zou
use strict;
#use warnings;
use Getopt::Long;
use Bio::SeqIO;
my ($help, $gff, $transcripts, $proteins, $genename);
GetOptions(
    "help!"      	 => \$help,
    "gff=s"      	 => \$gff,
    "transcripts=s"      => \$transcripts,
    "proteins=s"      	 => \$proteins,
    "genename=s"	 => \$genename,
);
my $INFO = <<LINES;
Usage:
	 perl $0 -gff maker1.all.gff -transcripts maker1.all.maker.transcripts.fasta -proteins maker1.all.maker.proteins.fasta -genename OH
LINES
die $INFO if ($help);
if (   !defined $gff
    || !defined $transcripts
    || !defined $proteins
    || !defined $genename )
{
    print "ERROR : Please set all parameters!\n";
    exit;
}




open (PRO2, ">proteins.format") or die;
open (PLEN,">protein.info");
my %geneid;
my $j;
my $geneidkey;
#my $num=4;
my $num=5;
my $inj  = Bio::SeqIO->new(-file => "$proteins" ,-format => 'fasta');
while(my $obj=$inj->next_seq()){
    my $id=$obj->id;
    my $seq=$obj->seq;
    my $lenp=length($seq);
    if ($lenp < 30 ){
    	next;
    }
    $j++;
    my $zero=$num-length($j);
    my $zstr='0' x $zero;
    $geneidkey="$genename".$zstr."$j";
	$geneid{$id}=$geneidkey;
	
    print PRO2 ">$geneidkey $id\n$seq\n" ;
	#my $lenp=length($seq);
	print PLEN "$geneidkey\t$lenp\n";
}

close (PRO2);
close (PLEN);


my $ina  = Bio::SeqIO->new(-file => "$transcripts" ,-format => 'fasta');
open (TRAN2,">transcripts.format") or die;
open (GENEINFO,">gene.info");



while(my $obj=$ina->next_seq()){
	my $id=$obj->id;
	my $seq=$obj->seq; 

	#$geneidkey="$genename"."_"."$i";
	#$geneid{$id}=$geneidkey;
	if(not exists $geneid{$id}){
		next;
	}
	$geneidkey=$geneid{$id};
	print TRAN2 ">$geneidkey $id\n$seq\n" ;
	my $lenseq=length($seq);
	my $sstart=substr($seq,0,3);
	my $send=substr($seq,$lenseq-3,3);
	print GENEINFO "$geneidkey\t$lenseq\t$sstart\t$send\n" ;
	
}


close (TRAN2);
close (GENEINFO);




`grep maker $gff > gff.maker`;
open (GFF, "<gff.maker") or die;
open (GFF2, ">gff.format") or die;
while (<GFF>){
	chomp;
	my $gff_line = $_;
	
	my @sid=(split /\t/,$gff_line);
	my $sid0=$sid[8];
	my $sid1=(split /;/, $sid0)[0];
	my $sid2=(split /-mRNA/,$sid1)[0];
	my $sid3=(split /=/,$sid2)[1];
	if(not exists $geneid{"$sid3-mRNA-1"}){
		next;
	}
	$sid[8]=$geneid{"$sid3-mRNA-1"};
	my $newline=join("\t",@sid);
	print GFF2 "$newline\n";
	
}
close (GFF);

close (GFF2);

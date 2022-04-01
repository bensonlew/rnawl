#! /usr/bin/perl

use strict;
use warnings;

use Getopt::Long;
my %opts;
GetOptions (\%opts,"m=s","f=s","o=s");

my $usage = <<"USAGE";
        Program : $0
        Discription:extrac gene nuc seqs by metagene result.
        Usage:perl $0 [options]         
         -m     metagen result format:
                # [sequence name]
                # gc = [gc%], rbs = [rbs%]
                # self: [(b)acteria/(a)rchaea/(p)hage/unused(-)]
                [gene ID] [start pos.] [end pos.] [strand] [frame] [complete/partial]
                [gene score] [used model] [rbs start] [rbs end] [rbs score]                
        -f      fasta 
        -o      output


USAGE
die $usage if (!$opts{m}|| !$opts{f} ||!$opts{o});

open FA,"<$opts{f}";

my $seq;
my %fa;
while(<FA>){chomp;
        if(/>(\S+)/){
            $seq=$1;        
        }else{
            $fa{$seq}.=$_;        
        }
}
close FA;

open MG,"<$opts{m}";
open OUT,">$opts{o}";
my $sn;  #seq name
my $gn;  #gene name

while(<MG>){chomp;
           if($_=~/#\s+(\S+)/&&$_!~/#\s+gc|#\s+self/){
                  $sn =$1;#print $sn;
           }
           if($_=~/^(gene_\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\d+)/){
                  $gn=$1;my $s1=$2;my $s2=$3;my $sd=$4;my $f=$5;
                  $s1=$s1+$f if($sd=~/\+/);
                  $s2=$s2-$f if($sd=~/\-/);  
                  my $eseq= substr ($fa{$sn},$s1-1,$s2-$s1+1);
                  if($sd=~/-/){
         		 $eseq =~ tr/atcgATCG/tagcTAGC/;
         		 $eseq =reverse $eseq;	    
     		 }
     		 print OUT ">$sn\_$gn\n$eseq\n";                    
           }
}
close MG;
close OUT;








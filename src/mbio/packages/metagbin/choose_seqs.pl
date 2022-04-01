#! /usr/bin/perl

use strict;
use warnings;
my $VERSION = "2012-11-5";
use Getopt::Long;
my %opts;
GetOptions (\%opts,"f=s","fq=s","l=s","c=s","co=s","o=s","q=s","qa=i");

my $usage = <<"USAGE";
        Program : $0
        Version : $VERSION
        Contact : guo.yu\@majorbio.com        
        Discription:choose fasta or qual from a list reverse with "R|r|*|-" ,cut with "start end"
        Usage:perl $0 [options] 
                -f     file.fasta
		        -fq    file.fastq
                -q     qual file (choose qual file like fasta choosed)
                -l *   list of names
                -c     cat seqs to one [i/n for interContig]
                -co    one contig name
                -o *   output file name               
                -qa    min average qual needed (output a seqname list satisfied)
USAGE
die $usage if ( !($opts{f}||$opts{q}||$opts{fq})&&!($opts{l}&&$opts{o}));

open LIS,"<$opts{l}" or die "could not open file $opts{l}!";

my $seqn;
my %fas;
if($opts{f}){
   open FAS,"<$opts{f}" or die "could not open file $opts{f}!";
   while(<FAS>){ chomp;
         if(/>(\S+)/){ 
             $seqn = $1;   
         }else{
             $fas{$seqn}.=$_;
         }
  }
  close FAS;
  open OUT,">$opts{o}";
}

my (%fqhead,%fqseq,%fqinfo,%fqqual);
if($opts{fq}){
   open FQS,"<$opts{fq}" or die "could not open file $opts{fq}!";
   while(<FQS>){
   chomp;
   if(/^@(\S+)/){
   	$seqn = $1;
	$fqhead{$seqn}=$_;
	chomp($fqseq{$seqn}=<FQS>);
	chomp($fqinfo{$seqn}=<FQS>);
	chomp($fqqual{$seqn}=<FQS>);
	}
   }
   close FQS;
   open OUTFQ,">$opts{o}"
}

my %quals;
if($opts{q}){
         open QUAL,"<$opts{q}";
         while(<QUAL>){chomp;
                  if(/^>(\S+)/){
                        $seqn = $1;                                        
                  }else{
                        my @line =split(/\s+/,$_);
                        if($quals{$seqn}){
                             @{$quals{$seqn}} =(@{$quals{$seqn}},@line);
                        }else{
                             @{$quals{$seqn}} =@line;
                        }
                  }                          
         }
         close QUAL;
         open OUTQ ,">$opts{o}.qual";
}



my $cat;
my $qa ;
if($opts{qa}){
     $qa =$opts{qa};
     open OUTQL ,">$opts{o}.qual.qa.list" ;
}


while(<LIS>){
         
         /^(\S+)\s(.*)$/;
         my $li = $1;
         my $na = $1; my $info=$2;
         my $sl; #sum length
         my $seq;
         if($opts{f}){
         if(!exists $fas{$li}){print "$li is not in your fasta file,disregard it.\n";next; }
              $seq = $fas{$li};  $sl=length($seq);
         }
         my @qual; my $qprint=1;
	 if($opts{fq}){
	 if(!exists $fqhead{$li}){print "$li is not in your fastq file,disregard it.\n";next; }
		$seq = $fqseq{$li};  $sl=length($seq);
	 }
         if($opts{q}){
              if(!exists $quals{$li}){print "$li is not in your qual file,disregard it.\n";next; }
              @qual =@{$quals{$li}};$sl=scalar(@qual);
         }
         my $rev =0;
	 my @infos=split(/\s+/,$info);	
         my @pos ;
         foreach my $if (@infos){ 
               push @pos,$if if($if=~/^\d+$/);
         }
         if(@pos>0){
               #$li .="\t@pos";
               my $s;my $e;
               if(@pos==2 && $pos[0] >$pos[1]){
                    $s =$pos[1]-1;$e=$pos[0]-1;
                    $rev=1;
               }elsif(@pos==2 && $pos[0] <=$pos[1]){
                    $s =$pos[0]-1;$e=$pos[1]-1                                        
               }else{
                    $s =$pos[0]-1;
                    $e =$sl-1; 
               }               
               my $len =$e-$s+1;
               $seq= substr ($seq,$s,$len)if($opts{f});
	       if($opts{fq}){
		$fqseq{$li}=substr ($fqseq{$li},$s,$len);
		$fqqual{$li}=substr ($fqqual{$li},$s,$len);
		}
               if($opts{q}){
                       @qual =@qual[$s..$e];
                       #print OUTQ ">$li\n";
                       #print OUTQ "@qual\n";  
                       if($opts{qa}){
                               my $av=0;
                               foreach my $q(@qual){
                                    $av +=$q;
                               }$av =$av/@qual;
                               if($av>=$opts{qa}){
                                     print OUTQL "$na\n";
                               }else{	$qprint=0;}                       
                       }                                   
               }                
         }

         if(/\s+([Rr\*\-])\s/i) { $rev=1;#$li .="\t$1";
	 } ##qual file had not been reversed;
         
         if($rev==1){         
               $seq =~ tr/atcgATCG/tagcTAGC/;
               $seq =reverse $seq;	    
         }
         
         if($opts{c}){
             $cat .=$seq; 
             $cat .="NNNNNCACACACTTAATTAATTAAGTGTGTGNNNNN" if($opts{c}=~/i/);
         }         
         if($opts{f}) {
             print OUT ">$li\n$seq\n" if(!$opts{c});  
         }
	 if($opts{fq}) {
	     print OUTFQ "$fqhead{$li}\n$fqseq{$li}\n$fqinfo{$li}\n$fqqual{$li}\n";
	 }
	 if($opts{q}) {
             print OUTQ ">$li\n@qual\n" if($qprint==1);  
         }
}

if($opts{c} && !$opts{co}){ 
    print "Warn:no onecontig name,name will be set to \"allcontig\"\n";
    $opts{co}="allcontig";
}

if($opts{f}) {
   print OUT ">$opts{co}\n$cat\n" if($opts{c});
}

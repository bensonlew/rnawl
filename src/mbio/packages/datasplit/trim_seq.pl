#!/usr/bin/perl -w

use strict;
use warnings;

use Getopt::Long;
my %opts;
GetOptions (\%opts,"f=s","min=i","max=i","s=i","l=i","o=s","q=s","qa=s","lg=s");

my $usage = <<"USAGE";
        Program : $0
        Discription:trim fasta with length.
        Usage:perl $0 [options] 
                -f    file.fasta
                -q    file.qual
                -qa   min average qual needed after trimed,used with -q
                -min  trim seqs shorter than minlength default:0
                -max  trim seqs longer than maxlegth   default:length(seq)
                -s    start position substr seqs       default:1
                -l    length to substr seqs            default:length(seq)
                -o    outputfile
		-lg   logfile			       default: trim.log
USAGE
die $usage if ( !($opts{f}||$opts{q}&&$opts{o}));
die $usage if($opts{qa}&&!$opts{q});

$opts{s}=$opts{s}?$opts{s}:1;


my $f=$opts{f};

open (FASTA, "<$opts{f}") || die ("Could not open $opts{f}!\n"); 
#open (GOOD,">$f.good.fasta") or die ("Could not open $f.good.fasta!\n"); 
#open (BAD,">$f.bad.fasta") or die ("Could not open $f.bad.fasta!\n"); 



my $sn;
my $seq="";
my $seqlen=-1;
my $numgood=0;
my $numbad=0;
my $numshort=0;
my $numlong=0;

$opts{s}=$opts{s}?$opts{s}:1;
$opts{l}=$opts{l}?$opts{l}:0;
$opts{min}=$opts{min}?$opts{min}:0;
$opts{max}=$opts{max}?$opts{max}:0; 
$opts{lg}=$opts{lg}?$opts{lg}:"trim.log";

my $max=$opts{max};
my $len=$opts{l};

my %qf;
if($opts{q}){
         my $seqn;my @seq=();
         open QUAL,"<$opts{q}";
         open OUTQ ,">$opts{o}.qual";
         while(<QUAL>){chomp;
		if(/>/||eof) {	if(/>(\S+)/){$sn=$1;};    
				if(!$seqn) {$seqn=$sn;$qf{$seqn}=0;next;}
		                if(/>/) {$seqlen=scalar(@seq);}else{my @line =split(/\s+/,$_);@seq=(@seq,@line);$seqlen=scalar(@seq)};
		                $max=$seqlen if($opts{max}==0);
		                $len=$seqlen-$opts{s}+1 if($opts{l}==0);
		                                     
                            	my $sublen=($len+$opts{s}-1)<=$seqlen?$len:($seqlen-$opts{s}+1);
                                                                                            
                    		if($seqlen<$opts{min}||$seqlen>$max) {
                    			$numbad++;
                    		}else{ 	my $s1=$opts{s}-1;my $s2=$sublen+$opts{s}-2;
					my @qual=@seq[$s1..$s2];
					my $ql=join(" ",@qual);
                                	if($opts{qa}){
						my $av=0;
                               			foreach my $q(@qual){
                                    	        $av +=$q;
                                        	}$av =$av/@qual;
				        	if($av>=$opts{qa}){
						print OUTQ ">$seqn\n$ql\n";  $qf{$seqn}=1;
						}					
					}else{
						print OUTQ ">$seqn\n$ql\n";  $qf{$seqn}=1;
					}
		   		}
		    @seq=();$seqn=$sn;$qf{$seqn}=0;
	     }else {
	        	my @line =split(/\s+/,$_);@seq=(@seq,@line);	
	     }
	}
        close QUAL;close OUTQ;
}



if($opts{f}){
my $seqname;
open (GOOD,">$opts{o}") or die ("Could not open $opts{o}!\n");
while(<FASTA>) {chomp;
	if(/>/||eof) {		if(/>(\S+)/){$sn=$1}; if(!$seqname) {$seqname=$sn;next} ;				    
		                if(/>/) {$seqlen=length($seq);}else{$seq.=$_;$seqlen=length($seq)};
		                $max=$seqlen if($opts{max}==0);
		                $len=$seqlen-$opts{s}+1 if($opts{l}==0);
		                                     
                          	my $sublen=($len+$opts{s}-1)<=$seqlen?$len:($seqlen-$opts{s}+1);
                            
                                                                
                    		if($seqlen<$opts{min}) {  
                 		   	$numshort++;
					$numbad++;
                   		}elsif($seqlen>$max){
					$numlong++;
					$numbad++;
				}elsif($opts{q}&&$qf{$seqname}==0){
					$numbad++;
				}else{
                    	  		 print GOOD ">$seqname\n";
                    	   		 print GOOD substr($seq,$opts{s}-1,$sublen)."\n"; $numgood++;
                    		}                   	        				
		   		$seq="";$seqname=$sn;
	}else {
		$seq.=$_;		
	}
}
close GOOD;
}
#add 20160901
print "discarded $numshort less than $opts{min} reads.\ndiscarded $numlong bigger than $opts{max} reads.\nremain $numgood clean reads.\n";
#close BAD;


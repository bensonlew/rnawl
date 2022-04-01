#!/usr/bin/perl -w
# Author: Thomas Thiel
# Program name: prim_output.pl
# Description: converts the Primer3 output into an table

open (SRC,"<$ARGV[0]") || die ("\nError: Couldn't open Primer3 results file (*.p3out) !\n\n");
my $filename = $ARGV[0];
$filename =~ s/\.p3out//;
# open (IN,"<$ARGV[1]") || die ("\nError: Couldn't open source file containing MISA (*.misa) results ! \n\n");
open (OUT,">$filename.result") || die ("\nError: Couldn't create file !\n\n");

my ($seq_names_failed,$count,$countfailed);

print OUT "#Chr\tSSR.nr\tSSR type\tSSR\tSize\tStart\tEnd\t";
print OUT "FORWARD PRIMER1 (5'-3')\tTm\tGC(%)\tLength(bp)\tREVERSE PRIMER1 (5'-3')\tTm\tGC(%)\tLength(bp)\tPRODUCT size(bp)\t";
print OUT "FORWARD PRIMER1 (5'-3')\tTm\tGC(%)\tLength(bp)\tREVERSE PRIMER1 (5'-3')\tTm\tGC(%)\tLength(bp)\tPRODUCT size(bp)\t";
print OUT "FORWARD PRIMER1 (5'-3')\tTm\tGC(%)\tLength(bp)\tREVERSE PRIMER1 (5'-3')\tTm\tGC(%)\tLength(bp)\tPRODUCT size(bp)\n";

undef $/;
$/ = "=\n";
while (<SRC>)
  {
  my @misa= (/PRIMER_SEQUENCE_ID=(\S+)_(\S+)_(\S+)_(\S+)_(\S+)_(\S+)_(\S+)/);
  my $misa=join("\t",$misa[0],$misa[5],$misa[1],$misa[6],$misa[2],$misa[3],$misa[4]);
  /PRIMER_LEFT_0_SEQUENCE=(.*)/ || do {print OUT $misa,"\t",join("\t",split//,"-"x27),"\n";next;};
  my $info = "$1\t";                                      # sequence
  /PRIMER_LEFT_0_TM=(.*)/; $info .="$1\t";                # TM
  /PRIMER_LEFT_0_GC_PERCENT=(.*)/; $info .= "$1\t";      # GC
  /PRIMER_LEFT_0=\d+,(\d+)/; $info .= "$1\t";             # length
  /PRIMER_RIGHT_0_SEQUENCE=(.*)/;  $info .= "$1\t";       # sequence
  /PRIMER_RIGHT_0_TM=(.*)/; $info .= "$1\t";              # TM
  /PRIMER_RIGHT_0_GC_PERCENT=(.*)/; $info .= "$1\t";     # GC
  /PRIMER_RIGHT_0=\d+,(\d+)/; $info .= "$1\t";            # length
  /PRIMER_PAIR_0_PRODUCT_SIZE=(.*)/; $info .= "$1\t";     # PRODUCT_SIZE
  
  
  /PRIMER_LEFT_1_SEQUENCE=(.*)/; $info .= "$1\t";         # sequence
  /PRIMER_LEFT_1_TM=(.*)/; $info .="$1\t";                # TM
  /PRIMER_LEFT_1_GC_PERCENT=(.*)/; $info .= "$1\t";      # GC
  /PRIMER_LEFT_1=\d+,(\d+)/; $info .= "$1\t";             # length
  /PRIMER_RIGHT_1_SEQUENCE=(.*)/;  $info .= "$1\t";       # sequence
  /PRIMER_RIGHT_1_TM=(.*)/; $info .= "$1\t";              # TM
  /PRIMER_RIGHT_1_GC_PERCENT=(.*)/; $info .= "$1\t";     # GC
  /PRIMER_RIGHT_1=\d+,(\d+)/; $info .= "$1\t";            # length
  /PRIMER_PAIR_1_PRODUCT_SIZE=(.*)/; $info .= "$1\t";     # PRODUCT_SIZE

  /PRIMER_LEFT_2_SEQUENCE=(.*)/; $info .= "$1\t";         # sequence
  /PRIMER_LEFT_2_TM=(.*)/; $info .="$1\t";                # TM
  /PRIMER_LEFT_2_GC_PERCENT=(.*)/; $info .= "$1\t";      # GC
  /PRIMER_LEFT_2=\d+,(\d+)/; $info .= "$1\t";             # length
  /PRIMER_RIGHT_2_SEQUENCE=(.*)/;  $info .= "$1\t";       # sequence
  /PRIMER_RIGHT_2_TM=(.*)/; $info .= "$1\t";              # TM
  /PRIMER_RIGHT_2_GC_PERCENT=(.*)/; $info .= "$1\t";     # GC
  /PRIMER_RIGHT_2=\d+,(\d+)/; $info .= "$1\t";            # length
  /PRIMER_PAIR_2_PRODUCT_SIZE=(.*)/; $info .= "$1";     # PRODUCT_SIZE

  print OUT join("\t",$misa, $info),"\n";
  };

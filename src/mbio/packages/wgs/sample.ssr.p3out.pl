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
print OUT "FORWARD PRIMER1 (5'-3')\tTm\tGC(%)\tLength(bp)\tREVERSE PRIMER1 (5'-3')\tTm\tGC(%)\tLength(bp)\tPRODUCT size(bp)\n";

undef $/;
$/ = "=\n";
while (<SRC>)
  {
  my @misa= (/PRIMER_SEQUENCE_ID=(\S+)_(\S+)_(\S+)_(\S+)_(\S+)_(\S+)_(\S+)/);
  my $misa=join("\t",$misa[0],$misa[5],$misa[1],$misa[6],$misa[2],$misa[3],$misa[4]);
  my $yyy=$misa[6];
  my (@Forward,@Reverse,@Tm1,@GC1,@length1,@Tm2,@GC2,@length2,@product);
  for($i=0;$i<5;$i++){
  /PRIMER_LEFT_$i\_SEQUENCE=(.*)/ || do {next};push @Forward,$1;
  /PRIMER_LEFT_$i\_TM=(.*)/; push @Tm1,$1;
  /PRIMER_LEFT_$i\_GC_PERCENT=(.*)/;push @GC1,$1;
  /PRIMER_LEFT_$i=\d+,(\d+)/; push @length1,$1;

  /PRIMER_RIGHT_$i\_SEQUENCE=(.*)/;push @Reverse,$1;
  /PRIMER_RIGHT_$i\_TM=(.*)/; push @Tm2,$1;
  /PRIMER_LEFT_$i\_GC_PERCENT=(.*)/;push @GC2,$1;
  /PRIMER_RIGHT_$i=\d+,(\d+)/; push @length2,$1;

  /PRIMER_PAIR_$i\_PRODUCT_SIZE=(.*)/; push @product,$1;
#  /PRIMER_LEFT_$i=(\d+),\d+/; push @Start,$npos-$1;my $s=$npos-$1;
#  /PRIMER_RIGHT_$i=(\d+),\d+/; push @End,$s;
#  push @End,$s+length($ref)-1;
  }
  print OUT join("\t",$misa,join(";",@Forward),join(";",@Tm1),join(";",@GC1),join(";",@length1),join(";",@Reverse),join(";",@Tm2),join(";",@GC2),join(";",@length2),join(";",@product)),"\n";
  };

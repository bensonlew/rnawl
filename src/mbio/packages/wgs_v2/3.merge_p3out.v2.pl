#!/usr/bin/perl -w
# Author: Thomas Thiel
# Program name: prim_output.pl
# Description: converts the Primer3 output into an table

open (SRC,"<$ARGV[0]") || die ("\nError: Couldn't open Primer3 results file (*.p3out) !\n\n");
my $filename = $ARGV[0];
$filename =~ s/\.p3out//;
open (OUT,">$filename.result") || die ("nError: Couldn't create file !\n\n");
my $number = $ARGV[1];
my ($count,$count_failed,$num,$info);

print OUT "#CHROM\tPOS\tTotal number\tType\tRef\tAlt\tMarker size(bp)\tMarker start(bp)\tMarker end(bp)\t";
for (my$i=0;$i<$number;$i++){
	my $return = $i + 1;
	print OUT "FORWARD PRIMER$return (5'-3')\tTm(.C)\tGC(%)\tsize\tREVERSE PRIMER$return (5'-3')\tTm(.C)\tGC(%)\tsize\tPRODUCT$return size (bp)\tstart (bp)\tend (bp)\t";
}
print OUT "\n";


$/ = "=\n";
while (<SRC>)
  {
  my ($id,$pos,$type,$ref,$alt,$start) = (/PRIMER_SEQUENCE_ID=(\S+)_(\d+)_(\S+)_(\S+)_(\S+)_(\d+)/);
  $num++;
  my %stat;
  $end=(length$ref) + $pos - 1;
  for (my$n=0;$n<$number;$n++){
	  my $return = $n + 1;
	  /PRIMER_LEFT_$n\_SEQUENCE=(.*)/ || do{$stat{$return}{failed}++;$count_failed++ ; next;};
	  #print $1,"\n";die;
	  $stat{$return}{left_seq} = $1;
	  /PRIMER_LEFT_$n\_TM=(.*)/; $stat{$return}{left_tm} = $1;
	  /PRIMER_LEFT_$n\_GC_PERCENT=(.*)/; $stat{$return}{left_gc} = $1;
	  /PRIMER_LEFT_$n=(\d+),(\d+)/;
	  $stat{$return}{left_start} = $1 + $start;
	  $stat{$return}{left_size} = $2;
	  /PRIMER_RIGHT_$n\_SEQUENCE=(.*)/; $stat{$return}{right_seq} = $1; 
	  /PRIMER_RIGHT_$n\_TM=(.*)/; $stat{$return}{right_tm} = $1;
	  /PRIMER_RIGHT_$n\_GC_PERCENT=(.*)/;  $stat{$return}{right_gc} = $1; 
	  /PRIMER_RIGHT_$n=(\d+),(\d+)/;
	  $stat{$return}{right_end} = $1 + $start;
	  $stat{$return}{right_size} = $2;
	  /PRIMER_PAIR_$n\_PRODUCT_SIZE=(.*)/; $stat{$return}{product} = $1; 
  }

  $count++;
  print OUT join("\t",$id,$pos,$num,$type,$ref,$alt,length$ref,$pos,$end),"\t";
  foreach my$return(sort {$a<=>$b} keys%stat){
	  $stat{$return}{failed}||=0;
	  next if($stat{$return}{failed} ne "0");
	  print OUT join("\t",$stat{$return}{left_seq},$stat{$return}{left_tm},$stat{$return}{left_gc},$stat{$return}{left_size},$stat{$return}{right_seq},$stat{$return}{right_tm},$stat{$return}{right_gc},$stat{$return}{right_size},$stat{$return}{product},$stat{$return}{left_start},$stat{$return}{right_end}),"\t";
  }
  print OUT "\n";
}

print "\nPrimer modelling was successful for $count sequences.\n";
print "Primer modelling failed for $count_failed sequences.\n";

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

#print OUT "#CHROM\tPOS\tTotal number\tType\tRef\tAlt\tMarker size(bp)\tMarker start(bp)\tMarker end(bp)\t";
print $number;
print "*****************************";
#for (my$i=0;$i<$number;$i++){
#    print "----------------";
#	my $return = $i + 1;
#	print OUT "FORWARD PRIMER$return (5'-3')\tTm(.C)\tGC(%)\tsize\tself_end_th\t self_any_th\thairpin_th\tend_stability\tpenalty\tREVERSE PRIMER$return (5'-3')\tTm(.C)\tGC(%)\tsize\tself_end_th\t self_any_th\thairpin_th\tend_stability\tpenalty\tPRODUCT$return size (bp)\tstart (bp)\tend (bp)\tpenalty pair\tcompl_any_th\tcompl_end_th\tSequence\t";
#}
print OUT "FORWARD PRIMER$return (5'-3')\tTm(.C)\tGC(%)\tsize\tself_end_th\t self_any_th\thairpin_th\tend_stability\tpenalty\tREVERSE PRIMER$return (5'-3')\tTm(.C)\tGC(%)\tsize\tself_end_th\t self_any_th\thairpin_th\tend_stability\tpenalty\tPRODUCT$return size (bp)\tstart (bp)\tend (bp)\tpenalty pair\tcompl_any_th\tcompl_end_th\tSequence\t";
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

	  # added by binbinzhao@20200506
	  /SEQUENCE_TEMPLATE=(.*)/;$stat{$return}{sequence} = $1;
	  /PRIMER_LEFT_$n\_PENALTY=(.*)/;$stat{$return}{left_penalty} = $1;
	  /PRIMER_LEFT_$n\_SELF_END_TH=(.*)/;$stat{$return}{left_self_end} = $1;
	  /PRIMER_LEFT_$n\_SELF_ANY_TH=(.*)/;$stat{$return}{left_self_any} = $1;
	  /PRIMER_LEFT_$n\_HAIRPIN_TH=(.*)/;$stat{$return}{left_hairpin} = $1;
	  /PRIMER_LEFT_$n\_END_STABILITY=(.*)/;$stat{$return}{left_end_stability} = $1;

      # added by binbinzhao@20200506
	  /PRIMER_PAIR_$n\_PENALTY=(.*)/;$stat{$return}{pair_penalty} = $1;
	  /PRIMER_PAIR_$n\_COMPL_ANY_TH=(.*)/;$stat{$return}{pair_compl_any} = $1;
	  /PRIMER_PAIR_$n\_COMPL_END_TH=(.*)/;$stat{$return}{pair_end_any} = $1;

      # added by binbinzhao@20200506
	  /PRIMER_RIGHT_$n\_SELF_END_TH=(.*)/;$stat{$return}{right_self_end} = $1;
	  /PRIMER_RIGHT_$n\_SELF_ANY_TH=(.*)/;$stat{$return}{right_self_any} = $1;
	  /PRIMER_RIGHT_$n\_HAIRPIN_TH=(.*)/;$stat{$return}{right_hairpin} = $1;
	  /PRIMER_RIGHT_$n\_END_STABILITY=(.*)/;$stat{$return}{right_end_stability} = $1;
	  /PRIMER_RIGHT_$n\_PENALTY=(.*)/;$stat{$return}{right_penalty} = $1;


	  /PRIMER_RIGHT_$n\_SEQUENCE=(.*)/; $stat{$return}{right_seq} = $1;
	  /PRIMER_RIGHT_$n\_TM=(.*)/; $stat{$return}{right_tm} = $1;
	  /PRIMER_RIGHT_$n\_GC_PERCENT=(.*)/;  $stat{$return}{right_gc} = $1; 
	  /PRIMER_RIGHT_$n=(\d+),(\d+)/;
	  $stat{$return}{right_end} = $1 + $start;
	  $stat{$return}{right_size} = $2;
	  /PRIMER_PAIR_$n\_PRODUCT_SIZE=(.*)/; $stat{$return}{product} = $1; 
  }

  $count++;
#  print OUT join("\t",$id,$pos,$num,$type,$ref,$alt,length$ref,$pos,$end),"\t";
  foreach my$return(sort {$a<=>$b} keys%stat){
	  $stat{$return}{failed}||=0;
	  next if($stat{$return}{failed} ne "0");
	  print OUT join("\t",$stat{$return}{left_seq},$stat{$return}{left_tm},$stat{$return}{left_gc},$stat{$return}{left_size},
	  $stat{$return}{left_self_end}, $stat{$return}{left_self_any}, $stat{$return}{left_hairpin},$stat{$return}{left_end_stability},
	  $stat{$return}{left_penalty},$stat{$return}{right_seq},$stat{$return}{right_tm},$stat{$return}{right_gc},$stat{$return}{right_size},
	  $stat{$return}{right_self_end}, $stat{$return}{right_self_any}, $stat{$return}{right_hairpin},$stat{$return}{right_end_stability},
	  $stat{$return}{right_penalty}, $stat{$return}{product},$stat{$return}{left_start},$stat{$return}{right_end},
	  $stat{$return}{pair_penalty}, $stat{$return}{pair_compl_any}, $stat{$return}{pair_end_any}, uc($stat{$return}{sequence}),),"\n";
  }
  print OUT "\n";
}

print "\nPrimer modelling was successful for $count sequences.\n";
print "Primer modelling failed for $count_failed sequences.\n";

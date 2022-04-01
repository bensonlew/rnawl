#!/usr/bin/perl -w
use strict;
use Bio::Seq;
use Bio::SeqIO;

my $input = $ARGV[0];

my $seqio_obj = Bio::SeqIO->new(-file => $input, -format => "fasta" );

my $line_count;

# process multi-fasta sequences
while(my $seq_obj = $seqio_obj->next_seq){
  $line_count++;

#  print "===========================================================\n";
  # obtain id of the nt sequence
  my $id = $seq_obj->display_id;
  # print id of nt sequence
#  print "SEQ ID\t>\t", $id, "\n";
#  print "===========================================================\n";

  # use translate to convert amino acid squence to protein sequence from frame 0 (default)
  # 'complete' - do some checks for ORF
#  print ">>>>>>>> Frame 0 <<<<<<<\n";
  print ">", $id, "-0", "\n";
  my $prot_obj = $seq_obj->translate();
  # print the protein sequence
  print $prot_obj->seq,"\n";
#  print "\n";
  
#  print ">>>>>>>> Frame 1 <<<<<<<\n";
  # translation starting from the second nucleotide (frame 1)
  print ">", $id, "-1", "\n";
  $prot_obj = $seq_obj->translate(-frame => 1);
  # print the protein sequence
  print $prot_obj->seq,"\n";
#  print "\n";

#  print ">>>>>>>> Frame 2 <<<<<<<\n";
  print ">", $id, "-2", "\n";
  # translation starting from the second nucleotide (frame 2)
  $prot_obj = $seq_obj->translate(-frame => 2);
  # print the protein sequence
  print $prot_obj->seq,"\n";
  
#  print "\n";
#  if($line_count == 1){last;}
}


#!perl -w
use strict;
use Bio::SeqIO;
use Bio::SeqFeature::Generic;

if(@ARGV != 5)  {
    print "perl GBK_generation.pl gene.gff3 tRNA.gff rRNA.gff Protein.faa Genome.fna \n";
    exit;
}
my($gene_gff3, $trna_gff3,$rrna_gff3,$protein_faa, $genome_faa,$output) = @ARGV;
my %protein_seqs;
my $all = "all.gff3";
`cat $trna_gff3 $rrna_gff3 $gene_gff3 > $all`;
open (IN2,$all) || die $!;
my %hash;
while(<IN2>){
   chomp;
   next if (/^Gene ID/);
   my @temp=split /\t/;
   my @seq_id = split /_ORF|_tRNA|_rRNA/,$temp[1];
   $hash{$seq_id[0]}{$temp[2]}=$_;
}
close IN2;
my $in = Bio::SeqIO->new(-file => " $protein_faa", -format =>"fasta");
while(my $seq = $in->next_seq())  {
    my $id = $seq->id;
    $protein_seqs{$id} = $seq->seq;
}
my $in1 = Bio::SeqIO->new(-file => "$genome_faa", -format => "fasta");
`mkdir -p gbk`;
while(my $seq = $in1->next_seq()) {
   my $id=$seq->id;
   my $out = Bio::SeqIO->new(-file => ">gbk/$id.gbk", -format => "genbank") ;
   my $feat3 = Bio::SeqFeature::Generic->new(
         -start     => 1,
         -end      => $seq->length,
         -primary => "source",
         -tag => {locus_tag =>$seq->id}
       );
   $seq->add_SeqFeature($feat3);
   foreach my $key (sort {$a cmp $b}keys %hash){
       if ($key eq $id){
           foreach my $key2 (sort {$a <=> $b} keys %{$hash{$key}}){
               my @temp = split /\t/,$hash{$key}{$key2};
               if ($temp[1] =~/_ORF/){
                   my $start = $temp[2];
                   my $end = $temp[3];
                   my $strand = $temp[4];
                   my $primary = "CDS";
                   my $locus_tag = $temp[0];
                   my $translation = $protein_seqs{$temp[0]};
                   my $codon_start = 1;
                   my $trans_table = 11;
                   my $feat1 = Bio::SeqFeature::Generic->new(
                      -start     => $start,
                      -end      => $end,
                      -primary => "gene",
                      -tag => {
                           locus_tag =>$locus_tag,
                              }
                      );
                   $seq->add_SeqFeature($feat1);
                   my $feat = Bio::SeqFeature::Generic->new(
                      -start     => $start,
                      -end      => $end,
                      -strand  => $strand,
                      -primary => $primary,
                      -tag => {
                      locus_tag =>$locus_tag,
                      conon_start => $codon_start,
                      transl_table => $trans_table,
                      translation => $translation,
                              }
                      );
                   $seq->add_SeqFeature($feat);
               }elsif($temp[1] =~/_rRNA/)  {
                   my $start = $temp[2];
	               my $end = $temp[3];
                   my $primary = 'rRNA';
                   my $locus_tag = $temp[0];
                   my $product;
                   my ($rRNA, $rRNA_type) = split /=/, $temp[7];
                   $product =$rRNA_type;
       	           my $feat = Bio::SeqFeature::Generic->new(
                       -start     => $start,
                       -end      => $end,
                       -primary => $primary,
                       -tag => {
                       locus_tag =>$locus_tag,
                       product => $product,
                    }
                   );
         	       $seq->add_SeqFeature($feat);
               }elsif($temp[1] =~/_tRNA/) {
    		       my $start = $temp[2];
    		       my $end = $temp[3];
    		       my $primary = 'tRNA';
    		       my $locus_tag = $temp[0];
    		       my $product;
    		       $product = "tRNA" . "-" . $temp[4];
   		           my $feat = Bio::SeqFeature::Generic->new(
                      -start     => $start,
                      -end      => $end,
                      -primary => $primary,
                      -tag => {
                         locus_tag =>$locus_tag,
                         product => $product,
                       }
                   );
        	       $seq->add_SeqFeature($feat);
               }
            }
         }
       }
      $out->write_seq($seq);
}


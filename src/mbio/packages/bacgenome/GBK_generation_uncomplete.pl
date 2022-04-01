#!perl -w
use strict;
use Bio::SeqIO;
use Bio::SeqFeature::Generic;

if(@ARGV != 7)  {
    print "perl GBK_generation.pl gene.gff3 tRNA.gff rRNA.gff Protein.faa Genome.fna anno OUTPUT \n";
    exit;
}
my($gene_gff3, $trna_gff3,$rrna_gff3,$protein_faa, $genome_faa, $anno, $output) = @ARGV;
my %protein_seqs;


open(OUT, ">$output.cgview.cog") || die "could not open cgvew.cog \n";
print OUT "seqname\tsource\tfeature\tstart\tend\tscore\tstrand\tframe\n";
my (%kegg,%besthit);
open(IN, $anno) || die "could not open $anno \n";
<IN>;
while(<IN>)  {
    chomp;
    my @temp = split /\t+/;
    if($temp[11] ne '-'){
       $kegg{$temp[0]}=$temp[11];
     }
       $besthit{$temp[0]}=$temp[6];
    if($temp[10] ne '-'){
        print OUT "$temp[0]\t\t$temp[10]\n";
     }
}
close IN;
`cat $trna_gff3 $rrna_gff3 $gene_gff3 > all.gff3`;
my $in = Bio::SeqIO->new(-file => " $protein_faa", -format =>"fasta");
while(my $seq = $in->next_seq())  {
    my $id = $seq->id;
    $protein_seqs{$id} = $seq;
}


my $in2 = Bio::SeqIO->new(-file => "$genome_faa", -format => "fasta");
my $ss ="all.fna";
open (DD,">$ss") || die $!;
my $seqs;
while(my $seq = $in2->next_seq()) {
   $seqs .=$seq->seq;
}
print DD ">Scaffold\n$seqs\n";
close DD;
my $in1 = Bio::SeqIO->new(-file => "$ss", -format => "fasta");
my $out = Bio::SeqIO->new(-file => ">$output.gbk", -format => "genbank") ;


while(my $seq = $in1->next_seq()) {
    my $id=$seq->id;
    open(IN, "all.gff3")  || die $!;
            my $feat3 = Bio::SeqFeature::Generic->new(
                                                                            -start     => 1,
                                                                            -end      => $seq->length,
                                                                            -primary => "source",
                                                                            -tag => {
                                                                                            locus_tag =>$seq->id,
          
                                                                                         }
                                                                            );
     $seq->add_SeqFeature($feat3);
    while(<IN>)  {
        chomp;
        next if (/^Gene ID/);
        my @temp = split /\t+/;
        my ($seq_id, $gene) = split /_/, $temp[1];
            if($gene =~/ORF/) {
                if(exists $besthit{$temp[0]})  {
                    my $start = $temp[7];
                    my $end = $temp[8];
                    my $strand = $temp[4];
                    my $primary = "CDS";
                    my $locus_tag = $temp[0];
                    my $translation = $protein_seqs{$temp[0]}->seq;
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
                       my $product = $besthit{$temp[0]};
                    my $feat = Bio::SeqFeature::Generic->new(
                                                                            -start     => $start,
                                                                            -end      => $end,
                                                                            -strand  => $strand,
                                                                            -primary => $primary,
                                                                            -tag => {
                                                                                            locus_tag =>$locus_tag,
                                                                                            conon_start => $codon_start,
                                                                                            transl_table => $trans_table,
                                                                                            product => $product,
                                                                                            translation => $translation,           
                                                                                         }
                                                                            );
                $seq->add_SeqFeature($feat);
               }
               else  {
                    my $start = $temp[7];
                    my $end = $temp[8];
                    my $strand = $temp[4];
                    my $primary = "CDS";
                    my $locus_tag = $temp[0];
                    my $translation = $protein_seqs{$temp[0]}->seq;
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
               }
            }
            elsif($gene=~/rRNA/)  {
                 my $start = $temp[8];
	          my $end = $temp[9];
        	my $primary = 'rRNA';
        	my $locus_tag = $temp[0];
        	my $product;
       	         my @temp1 = split /;/, $temp[7];
        	my ($rRNA, $rRNA_type) = split /=/, $temp1[1];
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
             }
             elsif($gene =~/tRNA/) {
    		my $start = $temp[9];
    		my $end = $temp[10];
    		my $primary = 'tRNA';
    		my $locus_tag = $temp[0];
    		my $product;
    		$product = "$gene" . "-" . $temp[4];
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
            $out->write_seq($seq);
}


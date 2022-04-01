#!perl -w
use strict;
use Bio::SeqIO;
use Bio::SeqFeature::Generic;

if(@ARGV != 6)  {
    print "perl GBK_generation.pl gene.gff3 tRNA.gff rRNA.gff Protein.faa Genome.fna anno \n";
    exit;
}
my($gene_gff3, $trna_gff3,$rrna_gff3,$protein_faa, $genome_faa, $anno, $output) = @ARGV;
my %protein_seqs;
my %hash_chr=('Chromosome'=>'Chr','Chromosome1'=>'Chr1',,'Chromosome2'=>'Chr2','Chromosome3'=>'Chr3','Chromosome4'=>'Chr4');
my %hash_pla=('Plasmid'=>'p','PlasmidA'=>'pA',,'PlasmidB'=>'pB','PlasmidC'=>'pC','PlasmidD'=>'pD','PlasmidE'=>'pE',,'PlasmidF'=>'pF','PlasmidG'=>'pG','PlasmidH'=>'pH');
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
   }
close IN;
my $all = "all.gff3";
my $all2 = "all2.gff3";
`cat $trna_gff3 $rrna_gff3 $gene_gff3 > $all`;
`cat  $all|sort -k 1 > $all2`;
my $in = Bio::SeqIO->new(-file => " $protein_faa", -format =>"fasta");
while(my $seq = $in->next_seq())  {
    my $id = $seq->id;
    $protein_seqs{$id} = $seq;
}

my $in1 = Bio::SeqIO->new(-file => "$genome_faa", -format => "fasta");

`mkdir -p gbk`;
my @num_id;
while(my $seq = $in1->next_seq()) {
   my $id=$seq->id;
    my ($des)=$id=~/[a-zA-Z]*([0-9]*)/;
    push @num_id,$des;
my $out = Bio::SeqIO->new(-file => ">gbk/$id.gbk", -format => "genbank") ;
    open(IN, "all2.gff3")  || die $!;
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
              if($seq_id eq $hash_chr{$id} || $seq_id eq $hash_pla{$id} || $seq_id eq $id){
                if($gene =~/ORF/) {
                if(exists $besthit{$temp[0]})  {
                    my $start = $temp[2];
                    my $end = $temp[3];
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
                    my $start = $temp[2];
                    my $end = $temp[3];
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
                 my $start = $temp[2];
	          my $end = $temp[3];
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
    		my $start = $temp[2];
    		my $end = $temp[3];
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
   }

              $out->write_seq($seq);
}
close IN;

my @files =glob "gbk/*.gbk";
 for my $dee (sort {$a <=> $b} @num_id){
  for my $file (sort @files){
    my @temp=split '\/',$file;
    my ($name2,$de) =split '\.',$temp[1];
    my ($des)=$name2=~/[a-zA-Z]*([0-9]*)/;
     if($dee == $des){
 `sed 's/ACCESSION   unknown/ACCESSION   $name2/g' $file -i`;
 `cat $file >>all.antismash.gbk`;
}      
}
}

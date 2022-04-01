#! perl -w
use strict;
use Bio::SeqIO;

my ($scaffolds,$contigs)=@ARGV;

my $output = "scaffolds.stat.xls";
&stat('scaffold','GH',$scaffolds,$output);
my $output2 = "contigs.stat.xls";
&stat('contig','GH',$contigs,$output2);


sub stat{
 my ($type,$sample,$input,$output)=@_;
 open (OUT,">$output") || die $!;
 my $in = Bio::SeqIO->new(-file => "$input", -format => "fasta");
 my $total;
 if($type eq  "scaffold"){
           `mkdir -p scaffold` ;
      }elsif($type eq "contig"){
           `mkdir -p contig` ;
      }elsif($type eq "gene"){
            `mkdir -p gene` ;
      }
 while(my $seq = $in->next_seq())  {
    my $id=$seq->id;
    my $seqs=$seq->seq;
    my $len=$seq->length;
    my $GC_number += $seq->seq =~ tr/cgCG/cgCG/;
    my $g_num += $seq->seq =~ tr/gG/gG/;
    my $c_num += $seq->seq =~ tr/cC/cC/;
    my $t_num += $seq->seq =~ tr/tT/tT/;
    my $a_num += $seq->seq =~ tr/aA/aA/;
    my $n_num += $seq->seq =~ tr/nN/nN/;
    my $total_sequence_length += $seq->seq =~ tr/atcgATCG/atcgATCG/;
    my $gc_rate = sprintf("%.2f",($GC_number/$total_sequence_length)*100);
    my $des;
    if($n_num != 0){
        $des="A:" . $a_num . "; T:" . $t_num . "; G:" . $g_num . "; C:" . $c_num . "; N:" . $n_num;
    }else{
        $des="A:" . $a_num . "; T:" . $t_num . "; G:" . $g_num . "; C:" . $c_num;
    }
    print OUT "$id\t$len\t$gc_rate\t$des\n";
    my $seq1 =Bio::Seq->new(-id =>$id,-seq =>$seqs);
    if($type eq  "scaffold"){
            my $out=Bio::SeqIO->new(-file =>">scaffold/$id.fasta",-format =>"fasta");
		    $out->write_seq($seq1);
    }elsif($type eq "contig"){
            my $out=Bio::SeqIO->new(-file =>">contig/$id.fasta",-format =>"fasta");
		    $out->write_seq($seq1);
    }elsif($type eq "gene"){
            my $out=Bio::SeqIO->new(-file =>">gene/$id.fasta",-format =>"fasta");
		    $out->write_seq($seq1);
      }
   }
 }

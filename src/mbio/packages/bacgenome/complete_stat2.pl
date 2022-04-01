#! perl -w
use Bio::SeqIO;

my ($scaf, $sample_name)=@ARGV;

my $output = "assemble";
&stat($scaf,$sample_name,$output);

sub stat{
 my ($input,$sample_name,$output)=@_;
 open (OUT,">$output.summary.xls") || die $!;
 open (OUT2,">$output.stat.xls") || die $!;
 print OUT2 "Chromosome No.\tPlasmid No.\tGenome Size (bp)\tG+C (%)\n";
 my $in = Bio::SeqIO->new(-file => "$input", -format => "fasta");
 my $total;
 `mkdir -p seq_dir` ;
 open (OUTL, ">seq_dir/list.txt") || die $!;
 print OUTL "sample\tfna\ttype\n";
 my $chr_num =0;
 my $pla_num =0;
 my $GC;
 while(my $seq = $in->next_seq())  {
    my $id=$seq->id;
    $id =~ s/(_pilon)+$//;
    my $desc=$seq->desc;
    if($id =~/chromosome/ || $id =~/Chromosome/){
      $chr_num++;
    }
     if($id =~/plasmid/ || $id =~/Plasmid/){
      $pla_num++;
    }
    my $seqs=$seq->seq;
    my $len=$seq->length;
    my $GC_number += $seq->seq =~ tr/cgCG/cgCG/;
     $GC +=$GC_number;
    my $g_num += $seq->seq =~ tr/gG/gG/;
    my $c_num += $seq->seq =~ tr/cC/cC/;
    my $t_num += $seq->seq =~ tr/tT/tT/;
    my $a_num += $seq->seq =~ tr/aA/aA/;
    my $n_num += $seq->seq =~ tr/nN/nN/;
    my $total_sequence_length += $seq->seq =~ tr/atcgATCG/atcgATCG/;
     $total +=$total_sequence_length;
    my $gc_rate = sprintf("%.2f",($GC_number/$total_sequence_length)*100);
    #my $des;
    #if($n_num != 0){
    #    $des="A:" . $a_num . "; T:" . $t_num . "; G:" . $g_num . "; C:" . $c_num . "; N:" . $n_num;
    #}else{
    #    $des="A:" . $a_num . "; T:" . $t_num . "; G:" . $g_num . "; C:" . $c_num;
    #}
    my $ss=ucfirst($id);
    print OUT "$id\t$len\t$gc_rate\t$desc\n";
    my $de = $sample_name . "." . $ss . ".fasta";
    print OUTL "$sample_name\t$de\t";
    if($id =~/chromosome/ || $id =~/Chromosome/){
      print OUTL "chromosome\n";
    }elsif($id =~/plasmid/ || $id =~/Plasmid/){
      print OUTL "plasmid\n";
    }else{
      print OUTL "\t\n";
    }
    my $seq1 =Bio::Seq->new(-id =>$ss,-seq =>$seqs);
    my $out=Bio::SeqIO->new(-file =>">seq_dir/$de",-format =>"fasta");
	$out->write_seq($seq1);
   }
   my $gc_rate2 = sprintf("%.2f",($GC/$total)*100);
   print OUT2 "$chr_num\t$pla_num\t$total\t$gc_rate2\n";
 }

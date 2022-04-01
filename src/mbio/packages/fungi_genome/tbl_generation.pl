#!perl -w
use Bio::SeqIO;

if(@ARGV != 7)  {
    print "perl $0 type gene_gff3 tRNA_gff3 rRNA_gff3 Protein.faa Genome.fna nr_anno \n";
    exit;
}
my($type_l,$gene_gff3,$tRNA_gff3, $rRNA_gff3, $protein_faa,$genome_faa,$nr_anno) = @ARGV;
my (%kegg,%nr_des);
open(IN, $nr_anno) || die "could not open $nr_anno \n";
<IN>;
while(<IN>)  {
    chomp;
    my @temp = split /\t+/;
    if($temp[14] ne '-'){
       $kegg{$temp[0]}=$temp[14];  ##[11] -> [14]
     }
       $nr_des{$temp[0]}=$temp[6];
}
close IN;
my %protein_seqs;
my $in = Bio::SeqIO->new(-file => " $protein_faa", -format =>"fasta");
while(my $seq = $in->next_seq())  {
    my $id = $seq->id;
    $protein_seqs{$id} = $seq->seq;
}
my %hash_chr=('Chromosome'=>'Chr','Chromosome1'=>'Chr1',,'Chromosome2'=>'Chr2','Chromosome3'=>'Chr3','Chromosome4'=>'Chr4');
my %hash_pla=('Plasmid'=>'p','PlasmidA'=>'pA',,'PlasmidB'=>'pB','PlasmidC'=>'pC','PlasmidD'=>'pD','PlasmidE'=>'pE',,'PlasmidF'=>'pF','PlasmidG'=>'pG','PlasmidH'=>'pH');
my $in1 = Bio::SeqIO->new(-file => "$genome_faa", -format => "fasta");
while(my $seq = $in1->next_seq()) {
     my $id =$seq->id; 
      `mkdir -p $id`;
       my $dee;
      if($type_l eq "complete"){
         $dee ="[organism=.] [strain=.] [topology=circular]";
      }elsif($type_l eq "uncomplete"){
         $dee ="[organism=.] [strain=.] [topology=linear]";
      }
     my $des =$seq->seq;
     my $seq1 =Bio::Seq->new(-id =>$id,-seq =>$des,-description=>$dee);
     my $out=Bio::SeqIO->new(-file =>">$id/$id.fsa",-format =>"fasta");
     $out->write_seq($seq1);
      open(OUT, ">$id/$id.tbl") || die "could not open $id.tbl \n";
     print OUT ">Feature $id\n";
    open(IN, $tRNA_gff3)  || die $!;
    while(<IN>){
      chomp;
     next if(/^Gene/);
     my @temp=split /\t/;
     my @temp2 = split /_/,$temp[1];
     my $seq = $temp2[0];
     my $type = 'tRNA';
     if($type_l eq "complete"){
          # print "$hash_chr{$id}\n";
     if($seq eq $hash_chr{$id} || $seq eq $hash_pla{$id} || $seq eq $id){
    print  OUT "$temp[2]\t$temp[3]\tgene\n";
    print  OUT "\t\t\tlocus_tag\t$temp[0]\n";
    print  OUT "$temp[2]\t$temp[3]\ttRNA\n";
    print OUT  "\t\t\tproduct\ttRNA-$temp[4]\n";
     }
    }elsif($type_l eq "uncomplete"){
      if($seq eq $id){
     print  OUT "$temp[2]\t$temp[3]\tgene\n";
    print  OUT "\t\t\tlocus_tag\t$temp[0]\n";
    print  OUT "$temp[2]\t$temp[3]\ttRNA\n";
    print OUT  "\t\t\tproduct\ttRNA-$temp[4]\n";

}
}    
} 
    open(IN, $rRNA_gff3)  || die $!;
    while(<IN>){
      chomp;
     next if(/^Gene/);
     my @temp=split /\t/;
     my ($seq,$type)=split /_/,$temp[1];
    if($type_l eq "complete"){
       if($seq eq $hash_chr{$id} || $seq eq $hash_pla{$id} || $seq eq $id ){
           my @arry=split/;/,$temp[7];
     my ($de,$des)=split /=/,$arry[1];
    print OUT  "$temp[2]\t$temp[3]\tgene\n";
    print OUT  "\t\t\tlocus_tag\t$temp[0]\n";
    print OUT  "$temp[2]\t$temp[3]\trRNA\n";
    print OUT  "\t\t\tproduct\t$des\n";
     }
    }elsif($type_l eq "uncomplete"){
      if($seq eq $id){ 
     my @arry=split/;/,$temp[7];
     my ($de,$des)=split /=/,$arry[1];
    print OUT  "$temp[2]\t$temp[3]\tgene\n";
    print OUT  "\t\t\tlocus_tag\t$temp[0]\n";
    print OUT  "$temp[2]\t$temp[3]\trRNA\n";
    print OUT  "\t\t\tproduct\t$des\n";

}
}    
}
  open(IN, $gene_gff3)  || die $!;
    while(<IN>){
      chomp;
     next if(/^Gene/);
     my @temp=split /\t/;
     my ($seq,$type)=split /_/,$temp[1];
      if($type_l eq "complete"){
        if($seq eq $hash_chr{$id} || $seq eq $hash_pla{$id} || $seq eq $id){
           print OUT  "$temp[2]\t$temp[3]\tgene\n";
           if(exists $kegg{$temp[0]}){
           print OUT  "\t\t\tgene\t$kegg{$temp[0]}\n";
             }
           print OUT  "\t\t\tlocus_tag\t$temp[0]\n";
           print OUT  "$temp[2]\t$temp[3]\tCDS\n";
           if(exists $kegg{$temp[0]}){
           print OUT  "\t\t\tgene\t$kegg{$temp[0]}\n";
              }
           print OUT  "\t\t\tlocus_tag\t$temp[0]\n";
           print OUT  "\t\t\tcodon_start\t1\n";
           print OUT  "\t\t\ttransl_table\t11\n";
           if(exists $nr_des{$temp[0]}){
           print OUT  "\t\t\tproduct\t$nr_des{$temp[0]}\n";
                }
        }else{
            next;
          }
    }elsif($type_l eq "uncomplete"){
      if($seq eq $id){ 
         print OUT  "$temp[2]\t$temp[3]\tgene\n";
        if(exists $kegg{$temp[0]}){
           print OUT  "\t\t\tgene\t$kegg{$temp[0]}\n";
             }
        print OUT  "\t\t\tlocus_tag\t$temp[0]\n";
        $exon = $temp[12];
        my @matches = $exon =~ /(\d+..\d+)/g;
        my $exon_num = 1;
        foreach my $my_exon (@matches){
            my @exon_temp = split(/\.\./,$my_exon);
            my $exon_strat = $exon_temp[0];
            my $exon_end = $exon_temp[1];
            if ($exon_num eq 1){
                print OUT  "$exon_strat\t$exon_end\tCDS\n";
            }else{
                print OUT  "$exon_strat\t$exon_end\n";
            }
            $exon_num += 1;
        }
        if(exists $kegg{$temp[0]}){
        print OUT  "\t\t\tgene\t$kegg{$temp[0]}\n";
            }
        print OUT  "\t\t\tlocus_tag\t$temp[0]\n";
        print OUT  "\t\t\tcodon_start\t1\n";
        print OUT  "\t\t\ttransl_table\t11\n";
        if(exists $nr_des{$temp[0]}){
        print OUT  "\t\t\tproduct\t$nr_des{$temp[0]}\n";
             }
        $exon = $temp[12];
        my @matches2 = $exon =~ /(\d+..\d+)/g;
        my $exon_num2 = 1;
        foreach my $my_exon (@matches2){
            my @exon_temp = split(/\.\./,$my_exon);
            my $exon_strat2 = $exon_temp[0];
            my $exon_end2 = $exon_temp[1];
            if ($exon_num2 eq 1){
                print OUT  "$exon_strat2\t$exon_end2\tmRNA\n";
            }else{
                print OUT  "$exon_strat2\t$exon_end2\n";
            }
            $exon_num2 += 1;
        }
        if(exists $kegg{$temp[0]}){
        print OUT  "\t\t\tgene\t$kegg{$temp[0]}\n";
            }
        print OUT  "\t\t\tlocus_tag\t$temp[0]\n";
        print OUT  "\t\t\tcodon_start\t1\n";
        print OUT  "\t\t\ttransl_table\t11\n";
        if(exists $nr_des{$temp[0]}){
        print OUT  "\t\t\tproduct\t$nr_des{$temp[0]}\n";
             }

     }else{
         next;
          }

}
}
}
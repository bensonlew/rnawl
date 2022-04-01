#!perl -w
use strict;

if(@ARGV != 4 and @ARGV != 5)  {
    print "perl generation_gff.pl gene.gff3 tRNA.gff rRNA.gff  anno \n";
    exit;
}

my $remove_gene;
if (@ARGV == 4){
    my($gene_gff3, $trna_gff3, $rrna_gff3, $anno) = @ARGV;
    $remove_gene = "-";
}else{
    my($gene_gff3, $trna_gff3, $rrna_gff3, $anno, $remove_gene) = @ARGV;
}


my($gene_gff3, $trna_gff3, $rrna_gff3, $anno) = @ARGV;


my (%kegg,%besthit);

open(IN, $anno) || die "could not open $anno \n";
<IN>;
while(<IN>){
    chomp;
    my @temp = split /\t/;
    #if($temp[11] ne '-'){
    $kegg{$temp[0]}=$temp[11];  ##11 chagne 14
    #}
    $temp[6]=~s/^ +//;
    $besthit{$temp[0]}= $temp[6];
   }
close IN;

my $all = "all.gff3";
my $all2 = "all2.gff3";
`cat $trna_gff3 $rrna_gff3 $gene_gff3 > $all`;
`cat  $all|sort -k 1 > $all2`;

if ( not -e  "gff"){`mkdir -p gff`};


open(IN, "all2.gff3")  || die $!;

open(OUT, ">gff/out.gff");
print OUT "##gff-version 3\n";

my $type;

while(<IN>){
    chomp;
    next if (/^Gene ID/);
    my @temp = split /\t/;
    my ($seq_id, $gene) = split /_/, $temp[1];

    my $start = $temp[2];
    my $end = $temp[3];

    if(int($start) > int($end)){
        my $t = $start;
        $start = $end;
        $end = $t ;
    }
    my $id = $temp[0];
    if ($remove_gene eq "-"){
        $gene =~ s/\D+/GeFixPreAbcdefg/;
    }
    if($gene =~/ORF/){
        $type = 'CDS';
        my $strand = $temp[4];
        my $name = $kegg{$temp[0]};
        my $exon = $temp[12];
        my @matches = $exon =~ /(\d+..\d+)/g;
        if(exists $besthit{$temp[0]}){
            my $product = $besthit{$temp[0]};
            if($name eq '-' or $name eq ""){
                print OUT "$seq_id\t.\tgene\t$start\t$end\t.\t$strand\t.\tID=$id;Locus_tag=$id\n";
                my $mrna_id = "rna" . "_" . $id;
                print OUT "$seq_id\t.\tmRNA\t$start\t$end\t.\t$strand\t.\tID=$mrna_id;Parent=$id;product=$product;Locus_tag=$id\n";
                my $exon_number = 1;
                foreach my $my_exon (@matches){
                    my @exon_temp = split(/\.\./,$my_exon);
                    my $exon_strat = $exon_temp[0];
                    my $exon_end = $exon_temp[1];
                    my $new_id = "exon" . "_" . $id . "_" . $exon_number;
                    my $parent_id = "rna" . "_" . $id;
                    print OUT  "$seq_id\t.\texon\t$exon_strat\t$exon_end\t.\t$strand\t.\tID=$new_id;Parent=$parent_id;product=$product;Locus_tag=$id\n";
                    $exon_number += 1;
                }
                foreach my $my_exon (@matches){
                    my @exon_temp = split(/\.\./,$my_exon);
                    my $exon_strat = $exon_temp[0];
                    my $exon_end = $exon_temp[1];
                    my $cds_id = "cds" . "_" . $id;
                    my $parent_id = "rna" . "_" . $id;
                    print OUT  "$seq_id\t.\tCDS\t$exon_strat\t$exon_end\t.\t$strand\t.\tID=$cds_id;Parent=$parent_id;product=$product;Locus_tag=$id\n";
                }
            }else{
                print OUT "$seq_id\t.\tgene\t$start\t$end\t.\t$strand\t.\tID=$id;Name=$name;Locus_tag=$id\n";
                my $mrna_id = "rna" . "_" . $id;
                print OUT "$seq_id\t.\tmRNA\t$start\t$end\t.\t$strand\t.\tID=$mrna_id;Parent=$id;Name=$name;product=$product;Locus_tag=$id\n";
                my $exon_number = 1;
                foreach my $my_exon (@matches){
                    my @exon_temp = split(/\.\./,$my_exon);
                    my $exon_strat = $exon_temp[0];
                    my $exon_end = $exon_temp[1];
                    my $new_id = "exon" . "_" . $id . "_" . $exon_number;
                    my $parent_id = "rna" . "_" . $id;
                    print OUT  "$seq_id\t.\texon\t$exon_strat\t$exon_end\t.\t$strand\t.\tID=$new_id;Parent=$parent_id;Name=$name;product=$product;Locus_tag=$id\n";
                    $exon_number += 1;
                }
                foreach my $my_exon (@matches){
                    my @exon_temp = split(/\.\./,$my_exon);
                    my $exon_strat = $exon_temp[0];
                    my $exon_end = $exon_temp[1];
                    my $cds_id = "cds" . "_" . $id;
                    my $parent_id = "rna" . "_" . $id;
                    print OUT  "$seq_id\t.\tCDS\t$exon_strat\t$exon_end\t.\t$strand\t.\tID=$cds_id;Parent=$parent_id;Name=$name;product=$product;Locus_tag=$id\n";
                }

            }
        }else{
            if($name eq '-' or $name eq ""){
                print OUT "$seq_id\t.\tgene\t$start\t$end\t.\t$strand\t.\tID=$id;Locus_tag=$id\n";
                my $mrna_id = "rna" . "_" . $id;
                print OUT "$seq_id\t.\tmRNA\t$start\t$end\t.\t$strand\t.\tID=$mrna_id;Parent=$id;Locus_tag=$id\n";
                my $exon_number = 1;
                foreach my $my_exon (@matches){
                    my @exon_temp = split(/\.\./,$my_exon);
                    my $exon_strat = $exon_temp[0];
                    my $exon_end = $exon_temp[1];
                    my $new_id = "exon" . "_" . $id . "_" . $exon_number;
                    my $parent_id = "rna" . "_" . $id;
                    print OUT  "$seq_id\t.\texon\t$exon_strat\t$exon_end\t.\t$strand\t.\tID=$new_id;Parent=$parent_id;Locus_tag=$id\n";
                    $exon_number += 1;
                }
                foreach my $my_exon (@matches){
                    my @exon_temp = split(/\.\./,$my_exon);
                    my $exon_strat = $exon_temp[0];
                    my $exon_end = $exon_temp[1];
                    my $cds_id = "cds" . "_" . $id;
                    my $parent_id = "rna" . "_" . $id;
                    print OUT  "$seq_id\t.\tCDS\t$exon_strat\t$exon_end\t.\t$strand\t.\tID=$cds_id;Parent=$parent_id;Locus_tag=$id\n";
                }

            }else{
                print OUT "$seq_id\t.\tgene\t$start\t$end\t.\t$strand\t.\tID=$id;Name=$name;Locus_tag=$id\n";
                my $mrna_id = "rna" . "_" . $id;
                print OUT "$seq_id\t.\tmRNA\t$start\t$end\t.\t$strand\t.\tID=$mrna_id;Parent=$id;Name=$name;Locus_tag=$id\n";
                my $exon_number = 1;
                foreach my $my_exon (@matches){
                    my @exon_temp = split(/\.\./,$my_exon);
                    my $exon_strat = $exon_temp[0];
                    my $exon_end = $exon_temp[1];
                    my $new_id = "exon" . "_" . $id . "_" . $exon_number;
                    my $parent_id = "rna" . "_" . $id;
                    print OUT  "$seq_id\t.\texon\t$exon_strat\t$exon_end\t.\t$strand\t.\tID=$new_id;Parent=$parent_id;Name=$name;Locus_tag=$id\n";
                    $exon_number += 1;
                }
                foreach my $my_exon (@matches){
                    my @exon_temp = split(/\.\./,$my_exon);
                    my $exon_strat = $exon_temp[0];
                    my $exon_end = $exon_temp[1];
                    my $cds_id = "cds" . "_" . $id;
                    my $parent_id = "rna" . "_" . $id;
                    print OUT  "$seq_id\t.\tCDS\t$exon_strat\t$exon_end\t.\t$strand\t.\tID=$cds_id;Parent=$parent_id;Name=$name;Locus_tag=$id\n";
                }
            }
        };
    }elsif($gene =~/rRNA/){
        my $type = 'rRNA';
        my $strand = $temp[5];
        my $name_product = $temp[7];
        my $trna_id = "rna" . "_" . $id;
        my $exon_number = 1;
        my $new_id = "exon" . "_" . $id . "_" . $exon_number;
        print OUT "$seq_id\t.\tgene\t$start\t$end\t.\t$strand\t.\tID=$id;Locus_tag=$id\n";
        print OUT "$seq_id\t.\t$type\t$start\t$end\t.\t$strand\t.\tID=$trna_id;Parent=$id;$name_product;Locus_tag=$id\n";
        print OUT "$seq_id\t.\texon\t$start\t$end\t.\t$strand\t.\tID=$new_id;Parent=$trna_id;$name_product;Locus_tag=$id\n";
    }elsif($gene =~/tRNA/){
        my $type = 'tRNA';
        my $strand = '.';
        my $product = $temp[4];
        my $rrna_id = "rna" . "_" . $id;
        my $exon_number = 1;
        my $new_id = "exon" . "_" . $id . "_" . $exon_number;
        print OUT "$seq_id\t.\tgene\t$start\t$end\t.\t$strand\t.\tID=$id;Locus_tag=$id\n";
        print OUT "$seq_id\t.\t$type\t$start\t$end\t.\t$strand\t.\tID=$rrna_id;Parent=$id;product=tRNA-$product;Locus_tag=$id\n";
        print OUT "$seq_id\t.\texon\t$start\t$end\t.\t$strand\t.\tID=$new_id;Parent=$rrna_id;product=tRNA-$product;Locus_tag=$id\n";
    }

}

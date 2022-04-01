#!perl -w
use strict;
#zouguanqing


if(@ARGV != 4 and @ARGV != 5)  {
    print "perl generation_gff.pl gene.gff3 tRNA.gff rRNA.gff  anno \n";
    exit;
}
my $remove_gene;
my $gene_gff3 = $ARGV[0];
my $trna_gff3 = $ARGV[1];
my $rrna_gff3 = $ARGV[2];
my $anno = $ARGV[3];

if (@ARGV == 4){
    $remove_gene = "-";
}else{
    $remove_gene = $ARGV[4];
}

my (%kegg,%besthit);

open(IN, $anno) || die "could not open $anno \n";
<IN>;
while(<IN>){
    chomp;
    my @temp = split /\t/;
    #if($temp[11] ne '-'){
    $kegg{$temp[0]}=$temp[14];  #11 change 14
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

my %data = ('Chr'=>'Chromosome', 'Chr1'=>'Chromosome2', 'Chr3'=>'Chromosome3','pA'=>'PlasmidA', 'pB'=>'PlasmidB', 'pC'=>'PlasmidC','pD'=>'PlasmidD', 'pE'=>'PlasmidE', 'pF'=>'PlasmidF','pG'=>'PlasmidG', 'pH'=>'PlasmidH', 'pI'=>'PlasmidI','pJ'=>'PlasmidJ', 'pK'=>'PlasmidK', 'pL'=>'PlasmidL','pM'=>'PlasmidM','pN'=>'PlasmidN','pO'=>'PlasmidO','pP'=>'PlasmidP');
my $seq_pre_id;
my $type;
while(<IN>){
    chomp;
    next if (/^Gene ID/);
    my @temp = split /\t/;
    my ($seq_id, $gene) = split /_/, $temp[1];
    if(exists $data{$seq_id}){
        $seq_id = $data{$seq_id};
    }
    if($seq_id ne $seq_pre_id){
        open(OUT, ">gff/$seq_id".".gff");
        print OUT "##gff-version 3\n";
        $seq_pre_id = $seq_id
    }

    my $start = $temp[2];
    my $end = $temp[3];

    if(int($start) > int($end)){
        my $t = $start;
        $start = $end;
        $end = $t ;
    }
    my $id = $temp[0];
    if ($remove_gene eq "-"){
        $id =~ s/\D+/GeFixPreAbcdefg/;
    }

    if($gene =~/ORF/){
        $type = 'CDS';
        my $strand = $temp[4];
        my $name = $kegg{$temp[0]};
        #if($name eq '-' or $name eq ""){$name=$id};
        if(exists $besthit{$temp[0]}){
            my $product = $besthit{$temp[0]};
            if($name eq '-' or $name eq ""){
            print OUT "$seq_id\t.\tgene\t$start\t$end\t.\t$strand\t.\tID=$id;Locus_tag=$id\n";
            print OUT "$seq_id\t.\t$type\t$start\t$end\t.\t$strand\t.\tID=$id;Parent=$id;product=$product;Locus_tag=$id\n";
            }else{
            print OUT "$seq_id\t.\tgene\t$start\t$end\t.\t$strand\t.\tID=$id;Name=$name;Locus_tag=$id\n";
            print OUT "$seq_id\t.\t$type\t$start\t$end\t.\t$strand\t.\tID=$id;Parent=$id;Name=$name;product=$product;Locus_tag=$id\n";
            }

        }else{
            if($name eq '-' or $name eq ""){
            print OUT "$seq_id\t.\tgene\t$start\t$end\t.\t$strand\t.\tID=$id;Locus_tag=$id\n";
            print OUT "$seq_id\t.\t$type\t$start\t$end\t.\t$strand\t.\tID=$id;Parent=$id;Locus_tag=$id\n";
            }else{
            print OUT "$seq_id\t.\tgene\t$start\t$end\t.\t$strand\t.\tID=$id;Name=$name;Locus_tag=$id\n";
            print OUT "$seq_id\t.\t$type\t$start\t$end\t.\t$strand\t.\tID=$id;Parent=$id;Name=$name;Locus_tag=$id\n";

            }

        };
    }elsif($gene =~/rRNA/){
        my $type = 'rRNA';
        my $strand = $temp[5];
        my $name_product = $temp[7];
        print OUT "$seq_id\t.\tgene\t$start\t$end\t.\t$strand\t.\tID=$id;Locus_tag=$id\n";
        print OUT "$seq_id\t.\t$type\t$start\t$end\t.\t$strand\t.\tID=$id;Parent=$id;$name_product;Locus_tag=$id\n";
    }elsif($gene =~/tRNA/){
        my $type = 'tRNA';
        my $strand = '.';
        my $product = $temp[4];
        print OUT "$seq_id\t.\tgene\t$start\t$end\t.\t$strand\t.\tID=$id;Locus_tag=$id\n";
        print OUT "$seq_id\t.\t$type\t$start\t$end\t.\t$strand\t.\tID=$id;Parent=$id;product=tRNA-$product;Locus_tag=$id\n";
    }

}


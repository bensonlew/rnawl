#!perl -w
use strict;
#zouguanqing


if(@ARGV != 2 and @ARGV != 3)  {
    print "perl generation_ptt.pl gene.gff3  anno \n";
    exit;
}
my $remove_gene;
my $gene_gff3 = $ARGV[0];
my $anno = $ARGV[1];
if (@ARGV == 2){
    $remove_gene = "-";
}else{
    $remove_gene = @ARGV[2];
}

my (%cog_id,%cog_code,%besthit);

open(IN, $anno) || die "could not open $anno \n";
<IN>;
while(<IN>){
    chomp;
    my @temp = split /\t/;
    #if($temp[11] ne '-'){
    $cog_id{$temp[0]}=$temp[9];
    $cog_code{$temp[0]} = $temp[10];
    #}
    $temp[6]=~s/^ +//;
    $besthit{$temp[0]}= $temp[6];
   }
close IN;


if (not -e "ptt"){`mkdir -p ptt`};


open(IN, "$gene_gff3")  || die $!;
#open(OUT, ">ptt/out.ptt");

my %data = ('Chr'=>'Chromosome', 'Chr1'=>'Chromosome2', 'Chr3'=>'Chromosome3','pA'=>'PlasmidA', 'pB'=>'PlasmidB', 'pC'=>'PlasmidC','pD'=>'PlasmidD', 'pE'=>'PlasmidE', 'pF'=>'PlasmidF','pG'=>'PlasmidG', 'pH'=>'PlasmidH', 'pI'=>'PlasmidI','pJ'=>'PlasmidJ', 'pK'=>'PlasmidK', 'pL'=>'PlasmidL','pM'=>'PlasmidM','pN'=>'PlasmidN','pO'=>'PlasmidO','pP'=>'PlasmidP');
my @result;
my $pid = 0;
my %out;
my $seq_id;
while(<IN>){
    chomp;
    next if (/^Gene ID/);

    my @temp = split /\t/;
    my ($seq_id_current, $gene1) = split /_/, $temp[1];
    if(exists $data{$seq_id_current}){
        $seq_id_current = $data{$seq_id_current};
    }
    if($seq_id ne $seq_id_current){
        if ($pid ne 0){
            my $num = @result;
            print OUT "$seq_id\n";
            print OUT "$num proteins\n";
            print OUT "Location\tStrand\tLength\tPID\tGene\tSynonym\tCode\tCOG\tProduct\n";
            foreach my $i (@result){print OUT $i;}
            print OUT "\n";
        }
        open(OUT, ">ptt/$seq_id_current".".ptt");
        undef @result;
        $seq_id = $seq_id_current;
    }
    $pid +=1 ;
    my $start = $temp[2];
    my $end = $temp[3];
    my $length = $temp[6];

    if(int($start) > int($end)){
        my $t = $start;
        $start = $end;
        $end = $t ;
    }
    my $gene = $temp[0];
    if ($remove_gene eq "-"){
        $gene =~ s/\D+/GeFixPreAbcdefg/;
    }
    my $cog_code_tmp = $cog_code{$temp[0]};
    my $cog_id_tmp = $cog_id{$temp[0]};
    my $strand = $temp[4];
    my $product = $besthit{$temp[0]};

    push @result,  "$start..$end\t$strand\t$length\t$pid\t$gene\t-\t$cog_code_tmp\t$cog_id_tmp\t$product\n";

}
my $num = @result;
print OUT "$seq_id\n";
print OUT "$num proteins\n";
print OUT "Location\tStrand\tLength\tPID\tGene\tSynonym\tCode\tCOG\tProduct\n";
foreach my $i (@result){ print OUT $i;}


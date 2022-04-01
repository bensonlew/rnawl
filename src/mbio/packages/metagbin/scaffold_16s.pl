#! /usr/bin/perl -w 
use Bio::SeqIO;

if (@ARGV != 4){
   print "perl $0 bin_dir all.scaffold.fa rRNA.gff output\n";
   exit;
}

my ($dir,$in,$rrna,$out)=@ARGV;
open (IN2,$rrna) || die $!;
open (OUT,">$out") || die $!;
print OUT "Scaffold Id\tBin Id\t16s_start\t16s_end\tStrand\tSeq\n";
my %hash;
my @files = glob "$dir/*.fa";
for my $de (@files){
    my $in = Bio::SeqIO->new(-file => "$de", -format => "fasta");
    my $file= (split /\//, $de)[-1];
    my ($name)=$file =~ /(.*)\.fa/;
    while(my $seq = $in->next_seq()) {
        my $id=$seq->id;
        my $len=$seq->length;
        $hash{$id}=$name;
    }
}

while(<IN2>){
   chomp;
   next if (/^#/);
   if($_=~/Name=16S_rRNA/){
    my @temp=split /\t/,$_;
    if (exists $rr{$temp[0]}){
      $rr{$temp[0]} .=";" . $temp[3] . "\t" . $temp[4] . "\t" . $temp[6];
}else{
    $rr{$temp[0]}=$temp[3] . "\t" . $temp[4] . "\t" . $temp[6];
}
}
}
close IN2;
my $seqs = Bio::SeqIO->new(-file => "$in", -format => "fasta");
while(my $seq = $seqs->next_seq()) {
    my $id = $seq->id;
    my $length = $seq->length;
    my  $sqes =$seq->seq;
    my $total_atcg = ($seq->seq=~ tr/atcgATCG/atcgATCG/);
    my $gc = ($seq->seq=~ tr/cgCG/cgCG/);
    my $de = $gc/$total_atcg;
    if(exists $rr{$id}){
        if($rr{$id} =~/;/){
            my @arry=split /;/,$rr{$id};
              for my $ds (@arry){
               my @gg=split /\t/,$ds;
               if($gg[2] eq '+'){
                 my $str = $seq->subseq($gg[0], $gg[1]);
                             if (exists $hash{$id}){
                   print OUT "$id\t$hash{$id}\t$gg[0]\t$gg[1]\t$gg[2]\t$str\n";
                               }else{
                   print OUT "$id\t-\t$gg[0]\t$gg[1]\t$gg[2]\t$str\n";
                            }
               }elsif($gg[2] eq '-'){
                 my $st = $seq->subseq($gg[0], $gg[1]);
                 my $str=reverse $st;
                 $str=~tr/agctAGCT/TCGATCGA/;
                   if (exists $hash{$id}){
                   print OUT "$id\t$hash{$id}\t$gg[0]\t$gg[1]\t$gg[2]\t$str\n";
                               }else{
                   print OUT "$id\t-\t$gg[0]\t$gg[1]\t$gg[2]\t$str\n";
                            }
               }
           }
        }else{
          my @gg=split /\t/,$rr{$id};
            my $str;
            if($gg[2] eq '+'){ 
                 $str = $seq->subseq($gg[0], $gg[1]);
               }elsif($gg[2] eq '-'){
                 my $st = $seq->subseq($gg[0], $gg[1]);
       	       	 $str=reverse	$st;
       	       	 $str=~tr/agctAGCT/TCGATCGA/;
               }
               if (exists $hash{$id}){
                   print OUT "$id\t$hash{$id}\t$gg[0]\t$gg[1]\t$gg[2]\t$str\n";
                               }else{
                   print OUT "$id\t-\t$gg[0]\t$gg[1]\t$gg[2]\t$str\n";
                            }
           }
     }
}


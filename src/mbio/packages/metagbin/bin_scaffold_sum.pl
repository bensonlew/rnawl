#! /usr/bin/perl -w 
use Bio::SeqIO;

if (@ARGV != 4){
   print "perl $0 bin_dir all.scaffold.fa depth.txt output\n";
   exit;
}

my ($dir,$in,$depth,$out)=@ARGV;
open (OUT,">$out") || die $!;
open (IN,$depth) || die $!;
print OUT "saffold_id\tBin_name\tstat\tend\tGC(%)\tcovergae\n";
my (%dep,%hash);
while(<IN>){ 
   chomp;
   next if(/^contigName/);
   my @temp=split /\t/;
   $dep{$temp[0]}=$temp[2];   
}
close IN;

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

my $seqs = Bio::SeqIO->new(-file => "$in", -format => "fasta");
while(my $seq = $seqs->next_seq()) {
    my $id = $seq->id;
    my $length = $seq->length;
    my $total_atcg = ($seq->seq=~ tr/atcgATCG/atcgATCG/);
    my $gc = ($seq->seq=~ tr/cgCG/cgCG/);
    my $de = $gc/$total_atcg;
    if(exists $hash{$id}){
          print OUT "$id\t$hash{$id}\t1\t$length\t$de\t$dep{$id}\n";
        }else{
          print OUT "$id\t-\t1\t$length\t$de\t$dep{$id}\n";
       }
}





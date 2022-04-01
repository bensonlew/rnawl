#! perl -W

use Bio::SeqIO;

my ($list,$dir)=@ARGV;

my @files =glob "$dir/*.fasta";
my %hash;
my %seq;

for my $de (@files){
    my @name=split /\//,$de;
    my $name =$name[-1];
    my $in = Bio::SeqIO->new(-file => "$de", -format => "fasta");
     my $len;
     my $seqs;
     while(my $seq = $in->next_seq())  {
      $len =$seq->length;
      $seqs =$seq->seq;
}
  $hash{$name}=$len;
  $seq{$name}=$seqs;
}
open (IN,$list) || die $!;
my %hash2;
<IN>;
while(<IN>){
  chomp;
  s/[\r\n]//g;
  my @temp =split /\t/;
my $de =@temp;
if($de ==3){
   $hash2{$temp[1]}=$temp[2];
}
}
close IN;

my (@chr,@pla);
foreach my $key (sort {$hash{$b} <=> $hash{$a}} keys %hash){
      if(exists $hash2{$key} && ($hash2{$key} eq 'chromosome' || $hash2{$key} eq 'Chromosome')){
         $des =$key . "\t" . $hash2{$key};
         push @chr,$des;

}elsif(exists $hash2{$key} && ($hash2{$key} eq 'plasmid' || $hash2{$key} eq 'Plasmid')){
           $des =$key . "\t" . $hash2{$key};
          push @pla,$des;
}
}

my %chr=&get_chr(@chr);
my %pla=&get_pla(@pla);
open (OUT,">chromosome.fasta") || dir $!;
foreach (sort keys %chr){
   if(exists $seq{$_}){
    my @temp=split /\t/,$chr{$_};
     my $de=ucfirst($temp[0]);
   print OUT ">$de\n$seq{$_}\n"
}
}
close OUT;
open (OUT2,">plasmid.fasta") || dir $!;
open (OUT3,">plasmid.type.xls") || dir $!;
foreach (sort keys %pla){
   print OUT3 "$pla{$_}\n";
   if(exists $seq{$_}){
    my @temp=split /\t/,$pla{$_};
    my $de=ucfirst($temp[0]);
   print OUT2 ">$de\n$seq{$_}\n";
}
}
close OUT2;
sub get_chr{
   my (@arry)=@_;
   my %new_chr; 
   my $de =@arry;
   if($de == 1){
   my @temp =split /\t/,$arry[0];
   $new_chr{$temp[0]} = $temp[1];
    }else{
    for(my $i=0;$i<$de;$i++){
    my @temp =split /\t/,$arry[$i];
     my $num =$i+1; 
   $new_chr{$temp[0]} = $temp[1] . $num;
    }
   }
 return %new_chr;
}

sub get_pla{
   my (@arry)=@_;
   my %hash =(1=>'A',2 =>'B',3 =>'C',4 =>'D',5 =>'E',6 =>'F',7 =>'G',8 =>'H',9 =>'I',10 =>'J', 11 =>'K', 12 =>'L', 13 =>'M',14 =>'N', 15 =>'O');
   my %new_pla; 
   my $de =@arry;
   if($de == 1){
   my @temp =split /\t/,$arry[0];
   my $de = ucfirst($temp[1]);
   $new_pla{$temp[0]} = $de . "\tp_gene";
}elsif($de != 1){
  for(my $i=0;$i< $de;$i++){
   my @temp =split /\t/,$arry[$i];
   my $de = ucfirst($temp[1]);
   $new_pla{$temp[0]}=$de . $hash{$i+1} . "\t" . "p" . $hash{$i+1} . '_gene';
}
}
return %new_pla; 
}

`cat chromosome.fasta plasmid.fasta >all.fasta`;

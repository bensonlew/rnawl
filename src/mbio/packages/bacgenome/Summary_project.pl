#! /usr/bin/perl -W
if (@ARGV != 8){
 print "perl $0 type sample assemble gene rrna trna kegg cog \n";
 exit;
}

my ($type,$sample,$assemble,$gene,$rrna,$trna,$kegg,$cog)=@ARGV;

open (IN,$assemble) || die $!;
open (OUT,">$sample.project.summary") || die $!;
my $pla;
my $chr;
my $g_size;
my $scf_no;
my $gc;
my $cog_no =0;
my $cds;
my $rrna_no =0;
my $kegg_no =0;
my $trna_no =0;
while(<IN>){
  chomp;
  next if(/^Total/);
  next if (/^Chromosome/);
  my @temp =split /\t/;
  if($type eq 'uncomplete'){
   $g_size =$temp[1];
   $scf_no = $temp[0];
   $gc = $temp[7];
  }elsif($type eq 'complete'){
   $chr=$temp[0];
   $pla =$temp[1];
   $g_size = $temp[2];
   $gc =$temp[3];
}
}

open (IN2,$gene) || die $!;

while(<IN2>){ 
  chomp;
  next if (/^Sample/);
  my @temp =split /\t/;
  $cds =$temp[1];
}

open (IN3,$rrna) || die $!;

while(<IN3>){
  chomp;
next if (/^Gene/);
  $rrna_no +=1;
}


open (IN4,$trna) || die $!;

while(<IN4>){
  chomp;
next if	(/^Gene/);
  $trna_no +=1;
}


open (IN5,$kegg) || die $!;

while(<IN5>){
  chomp;
next if (/^Gene/);
  $kegg_no +=1;
}


open (IN6,$cog) || die $!;
while(<IN6>){
  chomp;
next if (/^Gene/);
  $cog_no +=1;
}


if ($type eq 'uncomplete'){
  print OUT "Sample\tGenome Size\tscaffold no\tGC Content(%)\tCDS No.\trRNA No.\ttRNA No.\tGenes of KEGG\tGenes of COG\n";
  print OUT "$sample\t$g_size\t$scf_no\t$gc\t$cds\t$rrna_no\t$trna_no\t$kegg_no\t$cog_no\n";
}elsif($type eq 'complete'){
  print	OUT "Sample\tGenome Size\tChrom No.\tPlas No.\tGC Content(%)\tCDS No.\trRNA No.\ttRNA No.\tGenes of KEGG\tGenes of COG\n";
  print	OUT "$sample\t$g_size\t$chr\t$pla\t$gc\t$cds\t$rrna_no\t$trna_no\t$kegg_no\t$cog_no\n";
}
close IN;
close IN2;
close IN3;
close IN4;
close IN5;
close IN6;
close OUT;
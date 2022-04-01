#! /usr/bin/perl -W

if (@ARGV !=3){
  print "perl $0 go1234level_statistics.xls gene_profile all\n";
  exit;
}
my ($file,$file2,$level)=@ARGV;
my (%go1,%go12,%go123,%go1234,%go_sel_level);

my (%abund, @name);
if ($level =~ "all"){
open (IN2,$file2) || die $!;
my $name =<IN2>;
@name=split /\t/,$name;
shift(@name);
my $dess =join("\t",@name);
while(<IN2>){
  chomp;
  my @temp =split /\t/;
  my $des=join("\t",@temp[1..$#temp]);
  $abund{$temp[0]}=$des;
}
close IN2;

open (IN,$file) || die $!;
<IN>;
while(<IN>){
  chomp;
 my @temp=split /\t/;
 my $go1=join("\t",@temp[0..2]);
 my $go2=join("\t",@temp[0..4]);
 my $go3=join("\t",@temp[0..6]);
 my $go =$temp[0];
 if(exists $go1{$go}){
    $go1{$go} .=";" . $temp[9];
 }else{
    $go1{$go}=$temp[9];
 }
 if (exists $go12{$go1}){
   $go12{$go1} .=";" . $temp[9];
 }else{
   $go12{$go1}=$temp[9];
 }
 if (exists $go123{$go2}){ 
   $go123{$go2} .=";" . $temp[9];
 }else{ 
   $go123{$go2}=$temp[9];
 }
 if (exists $go1234{$go3}){ 
   $go1234{$go3} .=";" . $temp[9];
 }else{ 
   $go1234{$go3}=$temp[9];
 }
}
close IN;

my $go1_str=&go_anno(\%go1,\%abund);
my $go12_str=&go_anno(\%go12,\%abund);
my $go123_str=&go_anno(\%go123,\%abund);
my $go1234_str=&go_anno(\%go1234,\%abund);
open (OUT,">all.go1.function.xls")|| die $!;
print OUT "GO (Lev1)\t$dess$go1_str";
open (OUT,">all.go12.function.xls")|| die $!;
print OUT "GO (Lev1)\tGO Term (Lev2)\tGO ID (Lev2)\t$dess$go12_str";
open (OUT,">all.go123.function.xls")|| die $!;
print OUT "GO (Lev1)\tGO Term (Lev2)\tGO ID (Lev2)\tGO Term (Lev3)\tGO ID (Lev3)\t$dess$go123_str";
open (OUT,">all.go1234.function.xls")|| die $!;
print OUT "GO (Lev1)\tGO Term (Lev2)\tGO ID (Lev2)\tGO Term (Lev3)\tGO ID (Lev3)\tGO Term (Lev4)\tGO ID (Lev4)\t$dess$go1234_str";
}else{
    open (IN2,$file2) || die $!;
    my $name =<IN2>;
    @name=split /\t/,$name;
    shift(@name);
    splice (@name, -1);
    my $dess =join("\t",@name);
    while(<IN2>){
      chomp;
      my @temp =split /\t/;
      splice (@temp, -1);
      my $des=join("\t",@temp[1..$#temp]);
      $abund{$temp[0]}=$des;
    }
    close IN2;

    open (IN,$file) || die $!;
    <IN>;
    while(<IN>){
    chomp;
    my @temp=split /\t/;
    my ($go_level,$level_abbr);
    if($level eq "GO (Lev1)"){
        $go_level= $temp[0];
        $level_abbr = "l1";
    }elsif($level eq "GO Term (Lev2)"){
        $go_level= $temp[1];
        $level_abbr = "l2";
    }elsif($level eq "GO Term (Lev3)"){
        $go_level= $temp[3];
        $level_abbr = "l3";
    }elsif($level eq "GO Term (Lev4)"){
        $go_level= $temp[5];
        $level_abbr = "l4";
    }else{
      print "level name is wrong\n";
      exit;
    }
    if(exists $go_sel_level{$go_level}){
        $go_sel_level{$go_level} .=";" . $temp[9];
    }else{
        $go_sel_level{$go_level}=$temp[9];
    }
    }
    close IN;
    my $go_level_str=&go_anno(\%go_sel_level,\%abund);
    open (OUT,">level.xls")|| die $!;
    print OUT "$level\t$dess\n$go_level_str";
    close OUT;
}


sub go_anno{
  my ($go,$abu)=@_;
  my $tr;
  foreach (sort keys %$go){
     $tr .=$_;
       my @arry;
     if($$go{$_} =~/;/){
      my @temp=split /;/,$$go{$_};
        for $de (@temp){
        if(exists $$abu{$de}){
           my @temp3=split /\t/,$$abu{$de};
         push @arry,[@temp3];
     }
    }
   }else{
    if (exists $$abu{$$go{$_}}){
     my @temp3=split /\t/,$$abu{$$go{$_}};
      push @arry,[@temp3];
}
}
for my $j (0..$#name){
    my $sum=0;
     for my $i (0..$#arry){
      $sum +=$arry[$i][$j];
} 
 $tr .="\t" . $sum;
}
 $tr .= "\n";
}
return $tr;
}

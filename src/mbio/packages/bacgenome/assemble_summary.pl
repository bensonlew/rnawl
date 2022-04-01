#! /usr/bin/perl

#summary the contig status of an assembly


$scaffold = $ARGV[0];
$large = $ARGV[1];
$output =$ARGV[2];
open (OUT,">$output") || die $!;
if($large eq "") {
    $large = 1000;
}

open (SCAF, "<$scaffold") or die;


while (<SCAF>) {
    if($_ =~ /^>(\S+)/) {
	$name = $1;
	$names[@names] = $name;
	$id = @names -1;
    } else {
	chomp $_;
	$seq{$name} .= $_;
    }
}

@names = sort {length($seq{$b}) <=> length($seq{$a})} @names;

$largest_scaf = $names[0];
$largest_length = length($seq{$names[0]});
$len_large_scaf = 0;

for(local $i=0; $i < @names; ++$i) {

    @contseq = split /N+/, $seq{$names[$i]};
    for(local $n=0; $n<@contseq; ++$n) {
	$seq2{"$names[$i]-$n"} = $contseq[$n];
	$names2[@names2] = "$names[$i]-$n";
    }

    local $at = 0; $gc = 0; $basen = 0;
    for(local $j=0; $j< length($seq{$names[$i]}); ++$j) {
	$base = substr($seq{$names[$i]}, $j, 1);
	if($base =~ /[ATat]/) {
	    $at ++;
	} elsif ($base =~ /[GCgc]/) {
	    $gc ++;
	} else {
	    $basen ++;
	}
    }

    if(length($seq{$names[$i]}) >= $large) {
	$no_large_scaf ++;
	$len_large_scaf += length($seq{$names[$i]});
	$at_large_scaf += $at;
	$gc_large_scaf += $gc;
	$basen_large_scaf += $basen;
    }
    $no_scaf ++;
    $len_scaf += length($seq{$names[$i]});
    $at_scaf += $at;
    $gc_scaf += $gc;
    $basen_scaf += $basen;    
}

for(local $i=0; $i<@names; ++$i) {
    $length += length($seq{$names[$i]});
    if($length >= 0.5 * $len_large_scaf && $N50_large eq "") {
	$N50_large = length($seq{$names[$i]});
	$N50_large_num = $i+1;
    }
    if($length >= 0.9 * $len_large_scaf && $N90_large eq "") {
	$N90_large = length($seq{$names[$i]});
	$N90_large_num = $i+1;
    }
    if($length >= 0.5 * $len_scaf && $N50 eq "") {
	$N50 = length($seq{$names[$i]});
	$N50_num = $i+1;
    }
    if($length >= 0.9 * $len_scaf && $N90 eq "") {
	$N90 = length($seq{$names[$i]});
	$N90_num = $i+1;
    }
}


@names2 = sort {length($seq2{$b}) <=> length($seq2{$a})} @names2;

$largest_contig = $names2[0];
$largest_con_length = length($seq2{$names2[0]});

$N50_con_large_num="";
$N50_con_large = "";
$N50_length = "";
$N50_con_num="";
$con_N50 = "";
$N90_con_large_num="";
$N90_con_large = "";
$N90_con_length="";
$N90_con_num="";
$con_N90 = "";
$length=0;
$at = 0;
$gc = 0;
$basen = 0;

for(local $i=0; $i < @names2; ++$i) {

    local $at = 0; $gc = 0; $basen = 0;
    for(local $j=0; $j< length($seq2{$names2[$i]}); ++$j) {
	$base = substr($seq2{$names2[$i]}, $j, 1);
	if($base =~ /[ATat]/) {
	    $at ++;
	} elsif ($base =~ /[GCgc]/) {
	    $gc ++;
	} else {
	    $basen ++;
	}
    }

	
	if(length($seq2{$names2[$i]}) >= $large) {
	$no_large_contig ++;
	$len_large_contig += length($seq2{$names2[$i]});
	$at_large_contig += $at;
	$gc_large_contig += $gc;
	$basen_large_contig += $basen;
	}
    $no_contig ++;
    $len_contig += length($seq2{$names2[$i]});
    $at_contig += $at;
    $gc_contig += $gc;
    $basen_contig += $basen;
}

for(local $i=0; $i<@names2; ++$i) {
    $length += length($seq2{$names2[$i]});
    if($length >= 0.5 * $len_large_contig && $N50_con_large eq "") {
	$N50_con_large = length($seq2{$names2[$i]});
	$N50_con_large_num = $i+1;
    }
    if($length >= 0.9 * $len_large_contig && $N90_con_large eq "") {
	$N90_con_large = length($seq2{$names2[$i]});
	$N90_con_large_num = $i+1;
    }
    if($length >= 0.5 * $len_contig && $con_N50 eq "") {
	$con_N50 = length($seq2{$names2[$i]});
	$N50_con_num = $i+1;
    }
    if($length >= 0.9 * $len_contig && $con_N90 eq "") {
	$con_N90 = length($seq2{$names2[$i]});
	$N90_con_num = $i+1;
    }
}
my $scf_gc=int(($gc_large_scaf*100000)/($gc_large_scaf+$at_large_scaf))/1000;
my $n_rate=int($basen_scaf*100000/$len_scaf)/1000;


print OUT "Total Scaf No.\tTotal Bases in Scaf (bp)\tLarge Scaf No.\tLarge Scaf Bases(bp)\tLargest Scaf Len(bp)\tScaf N50(bp)\tScaf N90(bp)\tG+C(%)\tN Rate(%)\tTotal Ctg No.\tTotal Bases in Ctg (bp)\tLarge Ctg No.\tLarge Ctg Bases(bp)\tLargest Ctg Len(bp)\tCtg N50(bp)\tCtg N90(bp)\n";
print OUT "$no_scaf\t$len_scaf\t$no_large_scaf\t$len_large_scaf\t$largest_length\t$N50_large\t$N90_large\t$scf_gc\t$n_rate\t$no_contig\t$len_contig\t$no_large_contig\t$len_large_contig\t$largest_con_length\t$N50_con_large\t$N90_con_large\n";


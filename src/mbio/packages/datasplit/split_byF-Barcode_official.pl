#!/usr/bin/perl -w
use strict;

my ($inputFq,$config,$linkPrimer_th,$outPrefix);
#my $reversePrimer_th;

my $opt;
while($opt = shift){
	if($opt eq "-i"){
		$inputFq = shift;
	}elsif($opt eq "-f"){
		$config = shift;
	}elsif($opt eq "-l"){
		$linkPrimer_th = shift;
#	}elsif($opt eq "-r"){
#		$reversePrimer_th = shift;
	}elsif($opt eq "-o"){
		$outPrefix = shift;
	}elsif($opt eq "-h"){
		&usage();
		exit;
	}
}

unless($inputFq and $config and $outPrefix){
	&usage();
	exit;
}

#$barcode_th = 0 unless(defined($barcode_th));
$linkPrimer_th = 2 unless(defined($linkPrimer_th));
#$reversePrimer_th = -1 unless($linkPrimer_th);

my %abbrev=(
    'A' => "AA",    'T' => "TT",
    'C' => "CC",    'G' => "GG",
    'M' => "AC",    'R' => "AG",
    'W' => "AT",    'S' => "CG",
    'Y' => "CT",    'K' => "GT",
    'V' => "ACG",    'H' => "ACT",
    'D' => "AGT",    'B' => "CGT",
    'X' => "ACGT",    'N' => "ACGT",
    'I' => "ACGT",
);

my (%length,%sample,%linkPrimer,%Rbarcode);
#my %reversePrimer;

my (%randBase,%ReversePrimer);

open INC,$config or die "can\'t open sample config: $config\n";
while(<INC>){
	next if(/^\#/);
	my @temp = split;
	die "ERROR: F barcode $temp[1] used for more samples\n" if(exists $sample{$temp[1]}{$temp[2]} and $sample{$temp[1]}{$temp[2]} ne $temp[0]);
	die "ERROR: F barcode $temp[1] has different R barcode\n" if(exists $Rbarcode{$temp[1]} and $Rbarcode{$temp[0]} ne $temp[3]);
	my $rand;
	my $primer;
	if($temp[2] =~ /^(N*)([^N]+)/){
		$rand=$1; $primer=$2;
}
	#$sample{$temp[1]}{$temp[2]} = $temp[0];
	$sample{$temp[1]}{$primer} = $temp[0];
        $randBase{$temp[1]} = length($rand);
	$length{$temp[1]} = length($temp[1]);
#	$linkPrimer{$temp[1]} = $temp[2];
	$Rbarcode{$temp[1]} = $temp[3];
#	$reversePrimer{$temp[1]} = $temp[3];
	$ReversePrimer{$primer} = $temp[4];
}
close INC;


my $seqNum = 0;
open OUT,"> $outPrefix.split.fq" or die "$!\n";
open OUTE,"> $outPrefix.removed.fq" or die "$!\n";
if($inputFq =~ /gz$/){
	open INR,"gzip -dc $inputFq |" or die "$!\n";
}else{
	open INR,$inputFq or die "$!\n";
}
while(<INR>){
	$seqNum++;
	chomp;
	my $flag = 0;
	my $id = $_;
	$id =~ s/^@//;
	chomp(my $seq = <INR>);
	chomp(my $direction = <INR>);
	chomp(my $quality = <INR>);
	foreach my $Fbarcode(sort {$length{$a} <=> $length{$b}} keys %length){
		#if($seq =~ /^$Fbarcode/i){
		if($seq =~ /^\w/i){
#			print "forward sequencing...\n";
			$flag = 1;
			my $reFlag = 0;
			foreach my $linkPrimer(sort keys %{$sample{$Fbarcode}}){
				my $primer_mismatch = 0;
				for(my $i = 0;$i < length($linkPrimer);$i++){
					#my $seq_base2 = substr($seq,length($Fbarcode) + $i,1);
					my $seq_base2 = substr($seq,$i+$randBase{$Fbarcode},1);
					my $linkPrimer_base = substr($linkPrimer,$i,1);
					$primer_mismatch++ unless($abbrev{$linkPrimer_base} =~ /$seq_base2/i);
				}
				if($primer_mismatch <= $linkPrimer_th){
					#my $temp_seq = substr($seq,length($Fbarcode) + length($linkPrimer),length($seq) - length($Fbarcode) - length($linkPrimer));
					my $temp_seq = substr($seq,length($linkPrimer)+$randBase{$Fbarcode},length($seq)-length($linkPrimer)-$randBase{$Fbarcode});
					#my $temp_quality = substr($quality,length($Fbarcode) + length($linkPrimer),length($quality) - length($Fbarcode) - length($linkPrimer));
					my $temp_quality = substr($quality,length($linkPrimer)+$randBase{$Fbarcode},length($quality)-length($linkPrimer)-$randBase{$Fbarcode});
					#print OUT "\@$sample{$Fbarcode}{$linkPrimer}_$seqNum\t$id\torig_bc=$Fbarcode\tnew_bc=$Fbarcode\tbc_diffs=0\n$temp_seq\n$direction\n$temp_quality\n";
					print OUT "\@$sample{$Fbarcode}{$linkPrimer}_$seqNum\t$id\n$temp_seq\n$direction\n$temp_quality\n";
					$reFlag = 1;
					last;
				}
			}
			print OUTE "\@$id\tprimer\n$seq\n$direction\n$quality\n" unless($reFlag);
			last;
		}else{
			my $reverse_seq = reverse $seq;
			$reverse_seq =~ tr/ATGC/TACG/;
			$reverse_seq =~ tr/atgc/tacg/;
			my $reverse_quality = reverse $quality;
			next unless($reverse_seq =~ /^$Fbarcode/i);
			my $reFlag = 0;
			foreach my $linkPrimer(sort keys %{$sample{$Fbarcode}}){
				my $primer_mismatch = 0;
				for(my $i = 0;$i < length($linkPrimer);$i++){
					#my $seq_base2 = substr($reverse_seq,length($Fbarcode) + $i,1);
					my $seq_base2 = substr($reverse_seq,$i,1);
					my $linkPrimer_base = substr($linkPrimer,$i,1);
					$primer_mismatch++ unless($abbrev{$linkPrimer_base} =~ /$seq_base2/i);
				}
				if($primer_mismatch <= $linkPrimer_th){
					#my $temp_seq = substr($reverse_seq,length($Fbarcode) + length($linkPrimer),length($reverse_seq) - length($Fbarcode) - length($linkPrimer));
					my $temp_seq = substr($reverse_seq,length($linkPrimer),length($reverse_seq)- length($linkPrimer) -length($ReversePrimer{$linkPrimer}));
					#my $temp_quality = substr($reverse_quality,length($Fbarcode) + length($linkPrimer),length($reverse_quality) - length($Fbarcode) - length($linkPrimer));
					my $temp_quality = substr($reverse_quality,length($linkPrimer),length($reverse_quality) -length($linkPrimer) -length($ReversePrimer{$linkPrimer}));
					#print OUT "\@$sample{$Fbarcode}{$linkPrimer}_$seqNum\t$id\torig_bc=$Fbarcode\tnew_bc=$Fbarcode\tbc_diffs=0\n$temp_seq\n$direction\n$temp_quality\n";
					print OUT "\@$sample{$Fbarcode}{$linkPrimer}_$seqNum\t$id\n$temp_seq\n$direction\n$temp_quality\n";
					$reFlag = 1;
					last;
				}
			}
			print OUTE "\@$id\tprimer\n$seq\n$direction\n$quality\n" unless($reFlag);
			last;
		}
	}
	print OUTE "\@$id\tF-barcode\n$seq\n$direction\n$quality\n" unless($flag);
}
close INR;
close OUTE;
close OUT;


sub usage{
print <<EOD
	
	Description: split sequences for each sample by barcode
	Version: V1.20140214
	Contact: hua.chen\@majorbio.com

	usage: perl $0 -i removeBox.fq -f sample.config -l linkPrimer.cutoff -o out.prefix
		-i	fastq file removed box sequence,required
		-f	sample config file,include sampleID,barcode,link primer. eg: /mnt/lustre/users/chenhua/program/diversity/V2.0/data/split_byBarcode.config, required
		-l	link primer cutoff,default 2
		-o	output prefix,required
EOD
}

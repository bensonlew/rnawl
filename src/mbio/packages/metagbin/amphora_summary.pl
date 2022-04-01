#!usr/bin/perl
use strict;
use warnings;

die "usage: perl $0 bindir outDir\n" unless (@ARGV == 2);
my ($bindir,$outDir) = @ARGV;

#my @level = ("Superkingdom","Phylum","Class","Order","Family","Genus","Species");
my @level = ("d","p","c","o","f","g","s");
my @bin=glob "$bindir/*.amphora_anno.xls";
open OT2,">$outDir/summary.anno.xls" or die "can't find rout $outDir $!\n";
print OT2 "Bin\tTaxonomy\n";
for my $eachBin(@bin){
        my @Bin =split /\//,$eachBin;
        print "$Bin[-1]\n";
        my ($bin_name)=$Bin[-1]=~/(.*).amphora_anno.xls/;
	#open IN,"$eachBin" or die "can't find file $eachBin $!\n";
        my $num =&get_size($eachBin);
        print "$num\n";
if ($num <=1){
   print OT2 "$bin_name\t-\n";
}else{
open IN,"$eachBin" or die "can't find file $eachBin $!\n";
Loop:	while(<IN>){
		next if(/Superkingdom/);
		chomp;
		my @tmp = split /\t/,$_;
		my $name;
		for my $i(0..$#tmp){
			$tmp[$i] =~ s/\)//;
			my @ar = split(/\(/,$tmp[$i]);
			if($ar[1] == 1){
				$name .= $level[$i]."_".$ar[0].";";
			}else{
				print OT2 "$bin_name\t$name\n";
				last Loop;
			}
		}
		print OT2 "$bin_name\t$name\n";
	}
	close IN;
}
}
close OT2;

sub get_size{
             my $num=0;
            my ($file)=@_;
          open (IN,$file) || die $!;
       my @arry;
         while(<IN>){
            chomp;
        push @arry,$_;
}
$num =@arry;
return $num;
close IN;
}

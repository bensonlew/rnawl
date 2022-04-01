#perl -W

if(@ARGV != 3){
  print "perl $0 motifs.gff sample output \n";
  exit;
}

my ($motifs_gff,$sample,$out)=@ARGV;

open (IN,$motifs_gff) || die $!;
open (OUT,"> $out") || die $!;
print OUT "Sample Name\tLocation\tStrand\tmodificationType\tStart\tEnd\tScore\tmotifString\tCoverage\tIPDRatio\n";

while(<IN>){
    chomp;
    next if(/^#/);
    my @tmp=split /\t/;
    my $loc = $tmp[0];
    my $motiftype = $tmp[2];
    my $start = $tmp[3];
    my $end = $tmp[4];
    my $score = $tmp[5];
    my $strand = $tmp[6];
    my $detail = &get_each_detail($tmp[8]);
    if ($detail ne "notfind"){
        print OUT "$sample\t$loc\t$strand\t$motiftype\t$start\t$end\t$score\t$detail\n";
    }
}
close IN;

sub get_each_detail{
    my ($detail) = @_;
    my ($motif,$coverage,$ipdratio,$each_detail);
    #print "$detail\n";
    if ($detail =~ /motif=/){
        #print "$detail\n";
        $motif = &get_value($detail,"motif");
        $coverage = &get_value($detail,"coverage");
        $ipdratio = &get_value($detail,"IPDRatio");
        $each_detail = "$motif\t$coverage\t$ipdratio";
    }else{
        $each_detail = "notfind";
    }
    return $each_detail;
}

sub get_value{
    my ($detail,$mytype) = @_;
    my $value;
    $splitstring = $mytype."=";
    my @tmp1 = split($splitstring, $detail);
    my @tmp2 = split(";", $tmp1[1]);
    $value = $tmp2[0];
    if ($value eq ""){
        $value = "-";
    }
    return $value;
}
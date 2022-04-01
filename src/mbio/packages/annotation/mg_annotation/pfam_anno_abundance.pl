#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
my %opts;
GetOptions(\%opts, "q=s", "p=s", "des=s", "o=s");

my $usage = <<"USAGE";
        Program : $0
        Contact : qingchen.zhang\@majorbio.com
        Description: Calcualte pfam abundance.
        Usage:perl $0 [options]
                -q      imput gene_pfam_anno detail file
                -p      input gene profile file
                -o      output Dir

USAGE
die $usage if (!($opts{q}&&$opts{p}&&$opts{o}));

`mkdir -p $opts{o}` unless -e ($opts{o});

my (%gene_profile, @sams, %profile_genes);
open INP, $opts{p} or die "can not open $opts{p}!\n";
while(<INP>){
    chomp;
    if (/^GeneID/){
    @sams = split /\t/;
    shift @sams;# É¾³ýÁËheadµÄGeneIDÁÐ
    }else{
    my @tmp = split /\t/;
    my $gene = $tmp[0];
    $profile_genes{$gene} = 1;
    shift @tmp;
    my $tmp_profile = join("\t", @tmp);
    $gene_profile{$gene} = $tmp_profile;
    }
}
close INP;

my(%type, %clan, %pfam, %desc, %pfam_id);
open INF, $opts{q} or die "can not open $opts{q}!\n";
open OUTS,"> $opts{o}/gene_pfam_anno_stat.xls" or die "can not open $opts{o}/gene_pfam_anno_stat.xls!\n";
print OUTS "#Pfam ID\tDomain\tDomain Description\tType\tClan ID\n";
while (<INF>){
    chomp;
    next if (/^#/);
    my @tmp = split /\t/;
    my $query = $tmp[0];
    my $desc = $tmp[6];
    my $pfam_id = $tmp[4];
    print OUTS "$tmp[1]\t$pfam_id\t$desc\t$tmp[5]\t$tmp[8]\n";
    if (exists $profile_genes{$query}){
        my $type = $tmp[5];
        my $clan = $tmp[8];
        my $pfam = $tmp[1];
        #$type{$pfam_acc} = $type;
        #$clan{$pfam_acc} = $clan;
        #$desc{$pfam_acc} = $desc;
        #$pfam_id{$pfam_acc} = $pfam_id;
        $type{$type} .= ",".$query;
        $clan{$clan} .= ",".$query;
        $pfam{$pfam} .= ",".$query;
    }
}
close INF;
close OUTS;


open OUTC, "> $opts{o}/pfam_acc_profile.xls" or die "can not open $opts{o}/pfam_acc_profile.xls!\n";
open OUTT, "> $opts{o}/pfam_type_profile.xls" or die "can not open $opts{o}/pfam_type_profile.xls!\n";
open OUTG, "> $opts{o}/pfam_clan_profile.xls" or die "can not open $opts{o}/pfam_clan_profile.xls!\n";
my $samples = join("\t", @sams);
print OUTC "#Pfam ID\t","$samples\n";
print OUTT "#Type\t","$samples\n";
print OUTG "#Clan ID\t","$samples\n";

foreach my $eachpfam (sort keys %pfam){
    print OUTC "$eachpfam";
    my $pfam_gene = $pfam{$eachpfam};
    my @genes = split(/,/, $pfam_gene);
    shift @genes;
    for(my $i = 0; $i < @sams; $i++){
        my $pfam_abu;
        foreach my $gene (@genes){
            my @sams_profile = split(/\t/, $gene_profile{$gene});
            my $sam_profile = $sams_profile[$i];
            $pfam_abu += $sam_profile;
        }
    print OUTC "\t", $pfam_abu;
    }
    print OUTC "\n";
}
close OUTC;

foreach my $eachtype (sort keys %type){
    print OUTT "$eachtype";
    my $type_gene = $type{$eachtype};
    my @genes = split(/,/, $type_gene);
    shift @genes;
    for (my $i = 0; $i < @sams; $i++){
        my $type_abu;
        foreach my $gene (@genes){
            my @sams_profile = split(/\t/, $gene_profile{$gene});
            my $sam_profile = $sams_profile[$i];
            $type_abu += $sam_profile;
        }
        print OUTT "\t", $type_abu;
    }
    print OUTT "\n"
}
close OUTT;

foreach my $eachclan (sort keys %clan){
    print OUTG "$eachclan";
    my $clan_gene = $clan{$eachclan};
    my @genes = split(/,/, $clan_gene);
    shift @genes;
    for (my $i = 0; $i < @sams; $i++){
        my $clan_abu;
        foreach my $gene (@genes){
            my @sams_profile = split(/\t/, $gene_profile{$gene});
            my $sam_profile = $sams_profile[$i];
            $clan_abu += $sam_profile;
        }
    print OUTG "\t", $clan_abu;
    }
    print OUTG "\n"
}
close OUTG;
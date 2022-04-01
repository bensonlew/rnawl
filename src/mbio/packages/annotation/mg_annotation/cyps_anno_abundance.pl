#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
my %opts;
GetOptions(\%opts, "q=s", "p=s", "o=s");

my $usage = << "USAGE";
        Program : $0
        contact : qingchen.zhang\@majorbio.com
        Description : p450_anno_statistic.
        Usage:perl $0 [options]
                -q      input gene_cyps_anno detail file
                -p      input gene profile file
                -o      output Dir

USAGE

die $usage if (!($opts{q}&&$opts{o}&&$opts{p}));

`mkdir -p $opts{o}` unless -e ($opts{o});

my(%gene_profile, @sams, %profile_genes);
open INP, $opts{p} or die "can not open $opts{p}!\n";
while(<INP>){
    chomp;
    if (/^GeneID/){
    @sams = split /\t/;
    ####去除head的GeneID列和最后的total列
    shift @sams;
    #pop @sams;
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
##Query\tIdentity\tEvalue\thfam_id\tsfam_id\tnr_hit\tspecies\thomo_family\tsuper_family\n"

my (%homo,%homo_genes,%super,%super_genes,%hfid,%sid,%sfid,%sp,%gi);
open INF,$opts{q} or die "can not open $opts{q}!\n";
open OUTC,"> $opts{o}/cyps_anno_stat.xls" or die "can not open $opts{o}/cyps_anno_stat.xls!\n";
print OUTC "#NR_Hit\tSfam_ID\tHfam_ID\tSpecies\tHomologous_Family\tSuperfamily\n";

while(<INF>){
    chomp;
    next if (/^#/);
    my @tmp = split /\t/;
    my $query = $tmp[0];
    my $hfid = $tmp[4];
    my $sfid = $tmp[5];
    #my $sid = $tmp[1];
    my $sp = $tmp[7];
    my $gi = $tmp[6];
    print OUTC "$gi\t$sfid\t$hfid\t$sp\t$tmp[8]\t$tmp[9]\n";
    if (exists $profile_genes{$query}){
        my $homo = $tmp[8];
        my $super = $tmp[9];
        my $sid = $tmp[1];
        #$hfid{$sid} = $hfid;
        #$sfid{$sid} = $sfid;
        #$sp{$sid} = $sp;
        #$gi{$sid} = $gi;
        #$homo{$sid} = $homo;
        #$super{$sid} = $super;
        $homo_genes{$homo} .= ",".$query;
        $super_genes{$super} .= ",".$query;
        $sid{$sid} .= ",".$query;
    }
}
close INF;
close OUTC;


open OUTA,"> $opts{o}/cyps_homo_family_profile.xls" or die "can not open $opts{o}/cyps_homo_family_profile.xls!\n";
open OUTB,"> $opts{o}/cyps_super_family_profile.xls" or die "can not open $opts{o}/cyps_homo_family_profile.xls!\n";
open OUTD,"> $opts{o}/cyps_sid_profile.xls" or die "can not open $opts{o}/cyps_sid_profile.xls!\n";
my $samples = join("\t", @sams);
#print "$samples\n";
print OUTA "#Homo_family\t","$samples\n";
print OUTB "#Super_family\t","$samples\n";
print OUTD "#Sid\t", "$samples\n";
foreach my $eachhomo (sort keys %homo_genes){
    #print "-----------$eachhomo\n";
    print OUTA "$eachhomo";
    my $homo_gene = $homo_genes{$eachhomo};
    my @genes = split(/,/, $homo_gene);
    shift @genes;
    for(my $i = 0; $i < @sams; $i++){
        my $homo_abu;
        foreach my $gene(@genes){
            my @sams_profile = split(/\t/, $gene_profile{$gene});
            my $sam_profile = $sams_profile[$i];
            $homo_abu += $sam_profile;
        }
    print OUTA "\t", $homo_abu;
    }
    print OUTA "\n"
}
close OUTA;

foreach my $eachsuper (sort keys %super_genes){
    #print "-----------$eachsuper\n";
    print OUTB "$eachsuper";
    my $super_gene = $super_genes{$eachsuper};
    my @genes = split(/,/, $super_gene);
    shift @genes;
    for (my $i = 0; $i < @sams; $i++){
        my $super_abu;
        foreach my $gene(@genes){
            my @sams_profile = split(/\t/, $gene_profile{$gene});
            my $sam_profile = $sams_profile[$i];
            $super_abu += $sam_profile;
        }
    print OUTB "\t", $super_abu;
    }
    print OUTB "\n";
}
close OUTB;

foreach my $eachsid (sort keys %sid){
    print OUTD "$eachsid";
    my $sid_gene = $sid{$eachsid};
    my @genes = split(/,/, $sid_gene);
    shift @genes;
    for (my $i = 0;$i < @sams; $i++){
        my $sid_abu;
        foreach my $gene(@genes){
            my @sams_profile = split(/\t/, $gene_profile{$gene});
            my $sam_profile = $sams_profile[$i];
            $sid_abu += $sam_profile;
        }
        print OUTD "\t", $sid_abu;
    }
    print OUTD "\n";
}
close OUTD;














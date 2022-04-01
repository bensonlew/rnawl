#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
my %opts;
GetOptions(\%opts,"q=s","p=s","o=s");

my $usage = <<"USAGE";
        Program : $0
        Contact : shaohua.yuan\@majorbio.com
        Discription: Calculate nr taxon abundance.
        Usage:perl $0 [options]
                -q     input gene_ardb_anno detail file
                -p     input gene profile file
                -o     output Dir

USAGE

die $usage  if (!($opts{q}&&$opts{o}&&$opts{p}));

`mkdir -p $opts{o}` unless -e ($opts{o});


my (%gene_profile,@sams,%profile_genes);
open INP,$opts{p} or die "can not open $opts{p}!\n";
while(<INP>){
    chomp;
    if (/^GeneID/){
    @sams = split /\t/;
    #### 筛除head的GeneID列和最后的total列
    shift @sams;
    #pop @sams;
    }else{
    my @tmp = split /\t/;
    my $gene = $tmp[0];
    $profile_genes{$gene} = 1;
    shift @tmp;
    my $tmp_profile = join("\t",@tmp);
    $gene_profile{$gene} = $tmp_profile;
    }
}
close INP;

my (%Domain,%Kingdom,%Phylum,%Class,%Order,%Family,%Genus,%Species);
open INF,$opts{q} or die "can not open $opts{q}!\n";
while(<INF>){
    chomp;
    next if(/^#/);
    my @tmp = split;
    my $query = $tmp[0];
    if(exists $profile_genes{$query}){
        $Domain{$tmp[2]} .= ",".$query;
        $Kingdom{join(";",@tmp[2..3])} .= ",".$query;
        $Phylum{join(";",@tmp[2..4])} .= ",".$query;
        $Class{join(";",@tmp[2..5])} .= ",".$query;
        $Order{join(";",@tmp[2..6])} .= ",".$query;
        $Family{join(";",@tmp[2..7])} .= ",".$query;
        $Genus{join(";",@tmp[2..8])} .= ",".$query;
        $Species{join(";",@tmp[2..9])} .= ",".$query;
     }
}
close INF;


my $d_out = "$opts{o}/tax_d.xls";
my $k_out = "$opts{o}/tax_k.xls";
my $p_out = "$opts{o}/tax_p.xls";
my $c_out = "$opts{o}/tax_c.xls";
my $o_out = "$opts{o}/tax_o.xls";
my $f_out = "$opts{o}/tax_f.xls";
my $g_out = "$opts{o}/tax_g.xls";
my $s_out = "$opts{o}/tax_s.xls";
my $samples = join("\t",@sams);
&calculate_profile(\%Domain, $d_out);
&calculate_profile(\%Kingdom, $k_out);
&calculate_profile(\%Phylum, $p_out);
&calculate_profile(\%Class, $c_out);
&calculate_profile(\%Order, $o_out);
&calculate_profile(\%Family, $f_out);
&calculate_profile(\%Genus, $g_out);
&calculate_profile(\%Species, $s_out);

sub calculate_profile{
    my ($name,$output)=@_;
    open OUTF, ">$output" or die "can not open $output!\n";
    print OUTF "#Taxonomy\t","$samples\n";
    foreach my $eachname (sort keys %$name){
        print OUTF "$eachname";
        my $name_abu;
        my $name_gene = $$name{$eachname};
        my @genes = split(/,/,$name_gene);
        shift @genes;
        for(my $i = 0;$i < @sams ;$i++){
            my $name_abu;
            foreach my $gene (@genes){
            my @sams_profile = split(/\t/,$gene_profile{$gene});
            my $sam_profile = $sams_profile[$i];
            $name_abu += $sam_profile;
       }
        print OUTF "\t",$name_abu;
    }
    print OUTF "\n";
    }
    close OUTF;
}

#! perl

my $afile = $ARGV[0];
my $pfile = $ARGV[1];
my $outp = $ARGV[2];
my %idx;

open PF, '<', $pfile or die $!;
my @line;
while(1) {
    my $offset = tell PF;
    my $line = <PF>;
    last unless $line;
    @line = split /\t/, $line;
    $idx{$line[0]} = $offset;
}
close PF;


my (%KO_gene, %gene_KO, %gene, %KO, %KO_des, %KO_link, %KO_info, %path_query);
open LE, '<', $afile or die $!;
<LE>;
while (<LE>) {
    chomp;
    my @line = split '\t';
    $KO{$line[2]}{$line[0]} = 1;
    $KO_info{$line[2]} = "$line[3]\t$line[7]";

    $gene{$line[1]}{$line[0]} = 1;
    $gene2KO{$line[1]} = $line[2];

    $query_KO{$line[0]} = $line[2];
    map { $path_query{$_}{$line[0]} = 1 unless $_ eq "-" } split /\;/, $line[4];
}
close LE;


my (%module, %module_des, %t_abu);
my $md_f = $ARGV[3];
open MD, '<', $md_f or die $!;
while (<MD>) {
    chomp;
    my @line = split '\t';
    my @m = split ';', $line[1];
    my @m_des = split ';', $line[2];
    map {$module{$m[$_]}{$line[0]} = 1} 0..$#m;
    map {$module_des{$m[$_]} = $m_des[$_]} 0..$#m;
}
close MD;


my @per_kos;

sub outabu {
    my ($level, $outfile, $des, $h, $d) = @_;
    print "# $outfile @_#\n";
    open PF, '<', $pfile or die $!;
    open OUT, '>', $outfile or die $!;
    open OUT2, '>', "$outp/kegg_pathway_eachmap.xls" or die;
    my $header = <PF>;
    chomp $header;
    my @header = split /\t/, $header;
    my $has_t = lc($header[-1]) eq "total" ? 1 : 0;
    print OUT "$h\t",join "\t", @header[1..$#header], "$d\n";
    print OUT2 "#Pathway\tDescription\tKO_list\tTotal_abundance\tQuery_count\t", join "_map\t", @header[1..$#header];
    print OUT2 "_map\n";
    foreach my $le (keys %{$level}) {
        next if $le eq '-';
        my %genes = %{$level -> {$le}};
        my @tmp;
        my @line;
        foreach my $g(keys %genes) {
            next unless $idx{$g};
            seek(PF, $idx{$g}, 0);
            my $line = <PF>;
            chomp $line;
            @line = split "\t", $line;
            @tmp = map { $tmp[$_ - 1] + $line[$_] } 1 .. $#line;
            map {$pre_kos{$le}{$_ - 1}{$query_KO{$line[0]}} = 1 if $line[$_] > 0} 1..$#line;
        }
        next unless @tmp;
        my $t_abu;
        if ($total) {
            $t_abu = $tmp[-1];
        } else {
            map { $t_abu += $_ } @tmp;
        }
        my $pa = "http://www.genome.jp/kegg-bin/show_pathway?$le+";
        my %kos;
        my @pre_kos = map { my @k = sort keys %{$pre_kos{$le}{$_ - 1}}; map {$kos{$_}=1} @k; $pa . join "+", @k;} 1..$#line;
        my @kos = sort keys %kos;

        print OUT $le, "\t", join "\t", @tmp , $des -> {$le}, $pa . join "+", @kos;
        print OUT "\n";

        print OUT2 "$le\t$des{$le}\t", join ";", @kos;
        print OUT2 "\t$t_abu\t", 0 + keys %genes, join "\t", @pre_kos;
        print OUT2 "\n";
        @tmp = ();
    }
}

outabu(\%module, "$outp/kegg_pathway_profile.xls", \%module_des, "#Pathway", "Description\tPathwaymap");




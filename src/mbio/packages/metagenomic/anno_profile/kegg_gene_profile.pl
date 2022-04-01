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

sub outabu {
    my ($level, $outfile, $des, $h, $d) = @_;
    print "# $outfile @_#\n";
    open PF, '<', $pfile or die $!;
    open OUT, '>', $outfile or die $!;
    my $header = <PF>;
    chomp $header;
    my @header = split /\t/, $header;
    print OUT "$h\t",join "\t", @header[1..$#header], "$d\n";
    foreach my $le (keys %{$level}) {
        next if $le eq '-';
        my %genes = %{$level -> {$le}};
        my @tmp;
        foreach my $g(keys %genes) {
            next unless $idx{$g};
            seek(PF, $idx{$g}, 0);
            my $line = <PF>;
            chomp $line;
            my @line = split "\t", $line;
            @tmp = map { $tmp[$_ - 1] + $line[$_] } 1 .. $#line;
        }
        next unless @tmp;
        print OUT $le, "\t", join "\t", @tmp , $des -> {$le};
        print OUT "\n" ;
        @tmp = ();
    }
}

outabu(\%gene, "$outp/kegg_gene_profile.xls", \%gene2KO, "#Gene", "KO");

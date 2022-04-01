#! perl

my $pfile = $ARGV[1];
my $afile = $ARGV[0];
my $outp = $ARGV[2];
my %idx;

open PF, '<', $pfile or die $!;
my @line;
while(1) {
    # print "$offset\n";
    my $offset = tell PF;
    my $line = <PF>;
    last unless $line;
    @line = split /\t/, $line;
    $idx{$line[0]} = $offset;
    # print $line[0], "\t", $offset, "\n";
}
close PF;


my (%lv1, %lv2, %lv3, %des);
open LE, '<', $afile or die $!;
<LE>;
while (<LE>) {
    chomp;
    my @line = split '\t';
    $des{$line[1]} //= $line[1];
    $des{$line[2]} //= "$line[1]\t$line[2]";
    $des{$line[3]} //= "$line[1]\t$line[2]\t$line[3]";
    $lv1{$line[1]}{$line[0]} = 1;
    $lv2{$line[2]}{$line[0]} = 1;
    $lv3{$line[3]}{$line[0]} = 1;
}
close LE;

sub outabu {
    my ($level, $outfile, $max_l) = @_;
    print "# $outfile @_#\n";
    open PF, '<', $pfile or die $!;
    open OUT, '>', $outfile or die $!;
    my $header = <PF>;
    my @header = split /\t/, $header;
    my @h = map { "Level$_" } 1 .. $max_l;
    print OUT "#", join "\t", @h;
    print OUT "\t", join "\t", @header[1..$#header];
    # my %l = %{$level};
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
        print OUT $des{$le}, "\t", join "\t", @tmp;
        print OUT "\n";
        @tmp = ();
    }
}

outabu(\%lv1, "$outp/kegg_level1_profile.xls", 1);
outabu(\%lv2, "$outp/kegg_level2_profile.xls", 2);
outabu(\%lv3, "$outp/kegg_level3_profile.xls", 3);

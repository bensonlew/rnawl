#! perl
use List::Util qw(max);

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


my (%cog,%fun,%cate,%cog_des,%fun_des);
open LE, '<', $afile or die $!;
<LE>;
while (<LE>) {
    chomp;
    my @line = split '\t';
    my @cogs = split ';', $line[1];
    my @cog_dess = split ';', $line[2];
    my @funcs = split ';', $line[3];
    my @func_dess = split ';', $line[4];
    my @cates = split ';', $line[5];
    my $max = max ($#cogs, $#funcs, $#cates);
    map {
        $cog_des{$cogs[$_]} //= $cog_dess[$_];
        $fun_des{$funcs[$_]} //= $func_dess[$_];
        $cog{$cogs[$_]}{$line[0]} = 1;
        $fun{$funcs[$_]}{$line[0]} = 1;
        $cate{$cates[$_]}{$line[0]} = 1;
    } 0 .. $max;
}
close LE;

sub outabu {
    my ($level, $outfile, $type, $h) = @_;
    print "# $outfile @_#\n";
    open PF, '<', $pfile or die $!;
    open OUT, '>', $outfile or die $!;
    my $header = <PF>;
    chomp $header;
    my @header = split /\t/, $header;
    print OUT $h, "\t", join "\t", @header[1..$#header];
    $type eq "" ? print OUT "\n" : print OUT "\tDescription\n";
    # my %l = %{$level};
    foreach my $le (keys %{$level}) {
        next if $le eq '-';
        next unless $le;
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
        print OUT $le, "\t", join "\t", @tmp;
        $type eq "" ? print OUT "\n" : print OUT "\t", $type->{$le}, "\n";
        @tmp = ();
    }
}

outabu(\%cate, "$outp/cog_category_profile.xls", "", "#Category");
outabu(\%fun, "$outp/cog_function_profile.xls", \%fun_des, "#Function");
outabu(\%cog, "$outp/cog_nog_profile.xls", \%cog_des, '#NOG');

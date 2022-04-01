#!perl

my (%t, %out, %al);
while(<>){
    next if /^#/;
    chomp;
    @line = split /\s+/;
    next if $done{$line[3]} == 1;
    $t{$line[3]} //= $line[0];
    next if $t{$line[3]} ne $line[0];
    $out{$line[3]} = "$line[0]\t$line[3]\t$line[2]\t$line[6]\t$line[7]\t";
    $al{$line[3]} //= '1' x $line[2];
    my $a = '0' x ($line[16]-$line[15] + 1);
    my $d = $line[16]-$line[15] + 1;
    substr($al{$line[3]}, $line[15] - 1, $d, $a);
    #print STDERR "$line[3]\t$line[2]\t$al{$line[3]}\n";
    $done{$line[3]} = 1 if $line[9] eq $line[10];
}

print "arghub_id\tgene_id\talign_length\tevalue\tscore\tq_coverage\n";
foreach my $k(keys %out){
    my $l = length($al{$k});
    $al{$k}=~s/1//g;
    next if length($al{$k}) / $l < 0.2;
    print $out{$k}, length($al{$k}) / $l, "\n";

}
#next if /^#/;next if $done{$line[3]} == 1; $a{$line[3]} //= $line[0]; print "$line[0]\t$line[2]\t$line[3]\t$line[6]\t$line[7]\t$line[15]\t$line[16]\t$line[17]\t$line[18]\t", ($line[16]-$line[15])/$line[2] if $a{$line[3]} eq $line[0]; $done{$line[3]} = 1 if $line[9] eq $line[10]

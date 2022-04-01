#!perl

%opts = @ARGV;
$input = $opts{'-i'};
$reset = $opts{'-r'};

# $reset == 0 对输入文件的第一列用行号代替, 并保存与原始的对应关系
# $reset == 1 根据对应关系还原
if ($reset == 0) {
    open O1, ">reindex.table" or die;
    open O2, ">reindex.map" or die;

    open IN, "<$input", or die;
    my $header = <IN>;
    print O1 $header;
    while(<IN>){
        my @a = split "\t";
        my $row_name = $a[0];
        $a[0] = "newname$.";
        my $newline = join "\t", @a;
        print O1 $newline;
        print O2 "newname$.\t$row_name\n";
    }
} else {
    open I1, "<$input" or die "reindex.table does not exists";
    open I2, "<reindex.map" or die "reindex.table does not exists";
    my %index_map;
    while(<I2>) {
        chomp;
        my @l = split /\t/;
        $index_map{$l[0]} = $l[1];
    }
    
    open O, ">$input.tmp" or die;
    my $h = <I1>;
    print O $h;
    while(<I1>){
        s/^(\S+)/$index_map{$1}/ if /^(\S+)/ and $index_map{$1};
        print O
    }
    `mv $input.tmp $input`

}

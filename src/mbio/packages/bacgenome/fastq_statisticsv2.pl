#! /usr/bin/perl -w
#The author: xiaoyue.wang\@majorbio.com
use strict;
use warnings;

die "Usage : perl $0 <sample> <Q20.stats> <Q30.stats> <trim1.fastq> <output file> \n" unless(@ARGV==5);

my ($sample_name, $q20_f, $q30_f, $trim1, $output)=@ARGV if(@ARGV==5);

my $trimpair_read=&get_read($trim1)."*2";
my ($clean_totalBase, $clean_q20)=&open_Qfile($q20_f);
my ($cleantemp, $clean_q30)=&open_Qfile($q30_f);

(-s $output) || system ("echo -e \"Sample_lib\\tpair reads(#)\\ttotal bases(bp)\\tQ20(%)\\tQ30(%)\" >> $output");


open OUT,">>$output" || die $!;


print OUT "$sample_name\t$trimpair_read\t$clean_totalBase\t$clean_q20\t$clean_q30\n";


sub get_read{
        my $fq = $_[0];
        my $reads_num=0;
        open(FQ,"<$fq") or die;
        while(<FQ>){
                <FQ>;
                $reads_num++;
                <FQ>;
                <FQ>;
        }
        close(FQ);
        return $reads_num;
}
sub open_Qfile{
	my $file=$_[0];
	open QU,"<$file" || die "don't open file: $file!\n";
	chomp (my $trimfq1=<QU>);
	chomp (my $trimfq2=<QU>);
	my @line3=split /\t/,$trimfq1;
	my @line4=split /\t/,$trimfq2;
	my $clean_totalbase=$line3[1]+$line4[1];
	my $clean_Qvalue=sprintf("%.2f",($line3[2]+$line4[2])*100/$clean_totalbase);
	return ($clean_totalbase, $clean_Qvalue);
}


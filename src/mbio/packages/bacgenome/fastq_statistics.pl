#! /usr/bin/perl -w
#The author: xiaoyue.wang\@majorbio.com
use strict;
use warnings;

die "Usage : perl $0 <sample> <insert length> <Q20.stats> <Q30.stats> <fq1> <trim1.fastq> <trim_unpaired.fastq> <output file> \n" unless(@ARGV==8);

my ($sample_name, $avg_ins, $q20_f, $q30_f, $fq1_f, $trim1, $trim_single, $output)=@ARGV if(@ARGV==8);

my $rawpair_read=&get_read($fq1_f)."*2";
my $rawss = &get_read($fq1_f);
my $trimpair_read=&get_read($trim1)."*2";
my $trimsingle_read=&get_read($trim_single);
my ($raw_totalBase, $raw_q20, $clean_totalBase, $clean_q20)=&open_Qfile($q20_f);
my $read =$raw_totalBase/($rawss*2)-1;
my ($rawtemp, $raw_q30, $cleantemp, $clean_q30)=&open_Qfile($q30_f);

(-s $output) || system ("echo -e \"Sample\\tInsert length(bp)\\tRead Len\\tRaw data pair reads(#)\\tRaw data total bases(bp)\\tRaw data Q20(%)\\tRaw data Q30(%)\\\tClean data pair reads(#)\\tClean data single reads(#)\\tClean data total bases(bp)\\tClean data Q20(%)\\tClean data Q30(%)\" >> $output");
open OUT,">>$output" || die $!;

print OUT "$sample_name\t$avg_ins\t$read\t$rawpair_read\t$raw_totalBase\t$raw_q20\t$raw_q30\t$trimpair_read\t$trimsingle_read\t$clean_totalBase\t$clean_q20\t$clean_q30\n";

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
	chomp (my $rawfq1=<QU>);
	chomp (my $rawfq2=<QU>);
	chomp (my $trimfq1=<QU>);
	chomp (my $trimfq2=<QU>);
	chomp (my $singlefq=<QU>);
	my @line1=split /\t/,$rawfq1;
	my @line2=split /\t/,$rawfq2;
	my @line3=split /\t/,$trimfq1;
	my @line4=split /\t/,$trimfq2;
	my @line5=split /\t/,$singlefq;
	my $raw_totalbase=$line1[1]+$line2[1];
	my $raw_Qvalue=sprintf("%.2f",($line1[2]+$line2[2])*100/$raw_totalbase);
	my $clean_totalbase=$line3[1]+$line4[1]+$line5[1];
	my $clean_Qvalue=sprintf("%.2f",($line3[2]+$line4[2]+$line5[2])*100/$clean_totalbase);
	return ($raw_totalbase, $raw_Qvalue, $clean_totalbase, $clean_Qvalue);
}


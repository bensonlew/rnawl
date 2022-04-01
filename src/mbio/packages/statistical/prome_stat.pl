#!/usr/bin/perl -w
use strict;
use warnings;

die "perl $0 <proresult><cut><sample><out>\n" unless  @ARGV==4;

my ($proresult,$cut,$sample,$out)=@ARGV;
my $nextline;
open IN,$proresult||die die "can not open proresult" ;
open OUT, "> $out" || die "can not open output files for write";
print OUT "Gene ID\tLocation\tSample Name\tUpstream Pos\tPromoter Len\tPromoter Seq\tLsp\tLspe\tDmax Pos\tDmax\tDave\n";
while(<IN>) {
  chomp;
  if ($_ =~ /^ID/) {
      my @line = split /\s+/,$_;
      my @id = split /[+-]/,$line[1];
      my @info = split /\|/,$line[1];
      $nextline = readline IN;
      if ($nextline && $nextline =~/^>/){
          chomp $nextline;
          $nextline =~s/^>//;
          my @line1 = split /\s+/,$nextline;
          my @se = split /\.\./, $line1[0];
          if($se[0]>=$cut){
          next;
          }else{
          my $up = $cut - int($se[0]);
          my $lsp = $cut - int($line1[3]);
          my $dp = $cut - int($line1[5]);
          print OUT "$id[0]\t$info[1]\t$sample\t$up\t$line1[1]\t$line1[2]\t$lsp\t$line1[4]\t$dp\t$line1[6]\t$line1[7]\n"
          }
      }
      }
  }
#!/usr/bin/perl -w
use strict;
use warnings;
use Carp;
use Bio::SearchIO;
use Getopt::Long;
use FindBin qw($RealBin);
use DBI qw(:sql_types);
my %opts;
my $VERSION="2.2";
GetOptions( \%opts,"i=s", "format|f=s","o=s","maxEvalue|e=f","parse_id!","minIdentity|d=i","minCoverage|c=i","rank=i","img=s","db=s","mblast!","h!");
my $usage = <<"USAGE";
       Program : $0
       Version : $VERSION
       Contact : quan.guo\@majorbio.com
       Lastest modify:2012-06-06
       Discription:parse blast out file which blast with string databse result and mapping to COG/KOG/NOG
       Usage :perl $0 [options]
			-i*	inputfile	input file,can use wildcard character in bash,like '*_blast.out',but must use '' when using  wildcard character
			--format|-f		blast	input file format
											     blast      BLAST (WUBLAST, NCBIBLAST,bl2seq)   defualt
 												 fasta      FASTA -m9 and -m0
  												 blasttable BLAST -m9 or -m8 output (both NCBI and WUBLAST tabular)
  												 megablast  MEGABLAST
  												 psl        UCSC PSL format
												 waba       WABA output
												 axt        AXT format
  												 sim4       Sim4
  												 hmmer      HMMER2 hmmpfam and hmmsearch or HMMER3 hmmscan and hmmsearch
												 exonerate  Exonerate CIGAR and VULGAR format
												 blastxml   NCBI BLAST XML
												 wise       Genewise -genesf format
			-o	outputdir	output dir:defualt:cog_out
			--maxEvalue|-e	1e-5	Max E value,Defualt:1E-6
			--minIdentity|-d	75	Min percent(positives) identity over the alignment of COG/KOG/NOG region(hit),Defualt:75
			--minCoverage|-c	50	Min Coverage of COG/KOG/NOG region,Defualt:50
			-rank	10	rank cutoff for valid hits from BLAST result, Defualt:10
			-parse_id		parse_id from query description not the Query ID,for blastxml sometimes
			-img	COG,KOG	Image type,COG or KOG or NOG,split by ",",defualt:COG
			-db	/mnt/ilustre/users/sanger-dev/app/database/COG/cog.db		database name for sqlite,defualt:	/state/partition1/COG/cog.db or /mnt/ilustre/app/rna/database/COG/cog.db
			-mblast				for mblast result only
			-h			display help message

        exmaple:perl $0 -i 'unfinish_*.out' -format blastxml  -minIdentity 80 -parse_id
USAGE

die $usage if((! $opts{i}) || $opts{h} );

$opts{format}=$opts{format}?$opts{format}:"blast";
$opts{maxEvalue}=$opts{maxEvalue}?$opts{maxEvalue}:"1e-6";
$opts{minIdentity}=$opts{minIdentity}?$opts{minIdentity}:75;
$opts{minCoverage}=$opts{minCoverage}?$opts{minCoverage}:50;
$opts{rank}=$opts{rank}?$opts{rank}:10;
$opts{img}=$opts{img}?$opts{img}:'COG';
$opts{o}=$opts{o}?$opts{o}:"cog_out";
unless($opts{db}){
	$opts{db}="/mnt/ilustre/users/sanger-dev/app/database/COG/cog.db";
	# $opts{db}="/state/partition1/COG/cog.db" unless (-e "/mnt/ilustre/app/rna/database/COG/cog.db");
}

my @images=split(",",$opts{img});
foreach my $img (@images){
	die "img type input error!\n" unless ($img eq "COG" || $img eq "KOG" ||$img eq "NOG");
}
my @file= glob $opts{i};
warn("Input blast result files:\n");
warn(join("\n",@file)."\n");
mkdir("$opts{o}","493") or die "Can't create dir at $opts{o}\n" unless( -e $opts{o});
my $dbh = DBI->connect("dbi:SQLite:dbname=$opts{db}","","",{AutoCommit => 1});


my @seqs;
my %cogs;
my %lists;
open(TABLE, "> $opts{o}/cog_table.xls") || die "Can't open $opts{o}/cog_table.xls\n";
print TABLE "#Qeury_name\tQeury_length\tHsp_start_of_query\tHsp_end_of_query\tHsp_strand_of_query\tHit_name\tHit_description\tHit_length\tHsp_start_of_hit\tHsp_end_of_hit\tCOG/NOG_group\tCOG/NOG_group_description\tCOG/NOG_group_categeries\tCOG/NOG_region_start\tCOG/NOG_region_end\tCoverage_of_COG/NOG_region\tIdentities_of_COG/NOG_region\tPositives_Identities_of_COG/NOG_region\n";
foreach my $f (@file){
	warn("Parsing blast result file $f ...\n");
	my $searchio= Bio::SearchIO->new(-format => $opts{format},
								 -file => $f,
								 -best => 1,
								);
		while(my $result = $searchio->next_result){
			my $algorithm=$result->algorithm();
			die "Only support blastp and blastx result!\n" unless($algorithm=~/blastx|blastp/i);
			my $query_name=$result->query_name;
			if($opts{parse_id}||$query_name=~/^Query_\d+$/){
				$query_name=$result->query_description;
				$query_name=~/^\s*(\S+)/;
				$query_name=$1;
			}else{
				$query_name=$result->query_name;
				$query_name=~/^\s*(\S+)/;
				$query_name=$1;
			}
			push(@seqs,$query_name);
			my $query_length=$result->query_length;
			while(my $hit = $result->next_hit){
				last if $hit->rank() > $opts{rank};
				my $hit_length=$hit->length();
				my $score=$hit->score();
				my $hit_name=$hit->name();
				my $hash=&getRegion($hit_name);
				while(my $hsp= $hit->next_hsp){#Bio::Search::HSP::HSPI
					if($opts{'maxEvalue'}){
						last if $hsp->evalue > $opts{'maxEvalue'};
					}
					foreach my $cog (keys %$hash){
						my $coverage=&checkCoverage($hsp,$hash->{$cog}{'start'},$hash->{$cog}{'end'});
						next if($coverage<$opts{minCoverage});
						my ($identities,$positives)=&cogRegion($hsp,$hash->{$cog}{'start'},$hash->{$cog}{'end'});
						next if($positives<$opts{minIdentity});
						if($cog=~/COG/i){
							$cogs{$query_name}{"COG"}{$cog}=1;
						}elsif($cog=~/KOG/i){
							$cogs{$query_name}{"KOG"}{$cog}=1;
						}elsif($cog=~/NOG/i){
							$cogs{$query_name}{"NOG"}{$cog}=1;
						}
						my $strand;
						if($algorithm=~/blastx/i){
							if($opts{mblast}){
								if($hsp->start('hit')>$hsp->end('hit')){
									$strand="-";
								}else{
									$strand="+";
								}
							}else{
								$strand=$hsp->strand('query')==1?"+":"-";
							}
						}else{
							$strand=" ";
						}
						my $catlog="";
						# foreach my $c (keys(%{$hash->{$cog}{'categeries'}})){
						# 	$catlog.="[$c] ";
						# }
                        $catlog = join(",", keys(%{$hash->{$cog}{'categeries'}}));  # edited by shijin
						if(defined $cog){$lists{$query_name}=1}
						my $des=defined($hash->{$cog}{'discription'})?$hash->{$cog}{'discription'}:"_";
						#my ($id,$cons)=$hsp->matches(-seq => 'hit',-start=>$hash->{$cog}{'start'},-stop=>$hash->{$cog}{'end'});
						print TABLE "$query_name\t$query_length\t".$hsp->start('query')."\t".$hsp->end('query')."\t$strand\t$hit_name\t$hash->{$cog}{'annotation'}\t$hit_length\t".$hsp->start('hit')."\t".$hsp->end('hit')."\t$cog\t$des\t$catlog\t$hash->{$cog}{'start'}\t$hash->{$cog}{'end'}\t$coverage%\t$identities%\t$positives%\n";
					}
				}
			}
		}
}

close TABLE;


warn("Outputing list file ...\n");
open(LIST, "> $opts{o}/cog_list.xls") || die "Can't open $opts{o}/cog_list.xls\n";
print LIST "Query_name\tCOG\tNOG\n";
foreach my $q (@seqs){
	if(! exists $cogs{$q}){
		next;
	}
	print LIST "$q";
	foreach my $t (qw/COG NOG/){
		print LIST "\t".join(";",keys(%{$cogs{$q}{$t}}));
	}
	print LIST "\n";
}
close LIST;

my %sumary;
my $general=&getGeneral();
warn("Outputing sumary file ...\n");
foreach my $q (keys(%cogs)){
	foreach my $t (qw/COG NOG/){
		foreach my $c (keys(%{$cogs{$q}{$t}})){
			my $cat=&getCat($c);
			if($cat){
				foreach my $aa (@$cat){
					$sumary{$aa}{$t}{$q}=1;
				}
			}
		}
	}
}
#{$general->{$aa}{type}}
#3 $general->{$aa}{type}
open(SUMARY, "> $opts{o}/cog_summary.xls") || die "Can't open $opts{o}/cog_summary.xls\n";
print SUMARY "#Total seqs with COG/KOG/NOGs:".scalar(keys(%lists))."\n";
print SUMARY "#Type\tfunctional_categories\tCOG\tNOG\n";
my @types=('INFORMATION STORAGE AND PROCESSING','CELLULAR PROCESSES AND SIGNALING','METABOLISM','POORLY CHARACTERIZED');
#foreach my $t (@types){
foreach my $c (sort keys %sumary){
    print SUMARY "$general->{$c}{type}\t[$c] $general->{$c}{'name'}";
    foreach my $z (qw/COG NOG/){
        print SUMARY "\t".scalar(keys %{$sumary{$c}{$z}});
		}
    foreach my $z (qw/COG NOG/){
        print SUMARY "\t".join(";", keys % {$sumary{$c}{$z}});
    }
    print SUMARY "\n";
}

close SUMARY;
# warn("Generating image file ...\n");
# chdir("./$opts{o}/");
# foreach my $img (@images){
# open(PIC, "> pic_$img") || die "Can't open pic_$img\n";
# foreach my $t (@types){
# 	foreach my $c (sort keys %{$sumary{$t}}){
# 		my $tile=$general->{$c}{'name'};
# 		$tile=~/^\s*(\S+.*\S+)\s*$/;
# 		print PIC "$1\t".scalar(keys(%{$sumary{$t}{$c}{$img}}))."\n";
# 	}
# }
# close PIC;
# system("$RealBin/cog-bar.pl -i pic_$img");
# }

sub getRegion(){
	my $name=shift;
	my $query=$dbh->prepare(<<SQL
select `start`,`end`,`group`,`annotation` from mappings where seqid=?;
SQL
			  );
	my $discription=$dbh->prepare(<<SQL
select `name` from description where id=?;
SQL
	);
	my $categeries=$dbh->prepare(<<SQL
select `type` from categeries where id=?;
SQL
	);
	$query->execute($name);
	my $l = $query->fetchall_hashref('group');
	foreach my $cog ( keys(%$l) ){
		$discription->execute($cog);
		my $res=$discription->fetch();
		$l->{$cog}{'discription'}=$res->[0];
		$categeries->execute($cog);
		while(my $res1=$categeries->fetch()){
			$l->{$cog}{'categeries'}{$res1->[0]}=1;
		}
	}
	return $l;
}

sub getCat(){
	my $cog=shift;
	my @catelog;
	my $categeries=$dbh->prepare(<<SQL
select `type` from categeries where id=?;
SQL
	);
	$categeries->execute($cog);
	while(my $res1=$categeries->fetch()){
		push(@catelog,$res1->[0]);
	}
	return \@catelog;
}

sub getGeneral(){
	my $query=$dbh->prepare(<<SQL
select `id`,`name`,`type` from `general`;
SQL
			  );
	$query->execute();
	my $l = $query->fetchall_hashref('id');
	return $l;
}

sub checkCoverage(){
	my $hsp=shift;
	my $start=shift;
	my $end=shift;
	my $hsp_start;
	my $hsp_end;
	my ($s,$e);
	if($hsp->start('hit')>$hsp->end('hit')){
		$s=$hsp->end('hit');
		$e=$hsp->start('hit');
	}else{
		$s=$hsp->start('hit');
		$e=$hsp->end('hit');
	}
	if($opts{mblast}){
		my $hit_seq=$hsp->seq_str('hit');
		my $gap=0;
		for(my $i=0;$i<length($hit_seq);$i++){
			my $char=substr($hit_seq,$i,1);
			$gap++ if($char eq "-");
		}
		$e=$e-$gap;
	}
	if($s>$end || $e<$start){
		return 0;
	}
	if($s>$start||$e<$end){
		if($s>$start){
			$hsp_start=$s;
		}else{
			$hsp_start=$start;
		}
		if($e<$end){
			$hsp_end=$e;
		}else{
			$hsp_end=$end;
		}
		my $cover=sprintf("%.2f",($hsp_end-$hsp_start+1)*100/($end-$start+1));
		return $cover;
	}
	return 100;
}

sub cogRegion(){
	my $hsp=shift;
	my $start=shift;
	my $end=shift;
	my $identities=0;
	my $positives=0;
	my ($s,$e);
	if($hsp->start('hit')>$hsp->end('hit')){
		$s=$hsp->end('hit');
		$e=$hsp->start('hit');
	}else{
		$s=$hsp->start('hit');
		$e=$hsp->end('hit');
	}

	my $hit_seq=$hsp->seq_str('hit');
	if($opts{mblast}){
		my $gap=0;
		for(my $i=0;$i<length($hit_seq);$i++){
			my $char=substr($hit_seq,$i,1);
			$gap++ if($char eq "-");
		}
		$e=$e-$gap;
	}

	if($s>$end || $e<$start){
		return (0,0);
	}
		if($s>$start){
			$start=$s
		}
		if($e<$end){
			$end=$e
		}
	my $match_str=$hsp->seq_str('match');
	my $x=$start-$s;
	my $y=$end-$s;
	if($s<$start){
		for(my $i=0;$i<$start-$s;$i++){
			my $s=substr($hit_seq,$i,1);
			$x++ if($s eq "-");
		}
	}

	for(my $i=$x;$i<=$y;$i++){

		my $char=substr($hit_seq,$i,1);
		#print "$i $char\n";
		$y++ if($char eq "-");
		my $char1=substr($match_str,$i,1);
		if($char1=~/[A-Z]/){
			$identities++;
			$positives++;
		}elsif($char1=~/\+/){
			$positives++;
		}
	}
	$identities=sprintf("%.2f",$identities*100/($end-$start+1));
	$positives=sprintf("%.2f",$positives*100/($end-$start+1));
	return ($identities,$positives);
}

$dbh->disconnect();

#!/mnt/ilustre/centos7users/dna/.env/bin/env perl
#chao
#20190215

use strict;
use warnings;
use POSIX qw(setsid);
use Getopt::Long;
use File::Basename;
use File::Spec;

my $script = basename($0);
my $INFO=<<LINE;
Discription:
	Use pb-assembly to assembly pacbio data.
Usage:
	$script --infile dir_contains_subreads.fasta --cpu 20 --mem 50
Options:
	--help   NA		display this information
	--infile  STR	input directory contains subreads.fasta format files [required]
	--gsize  INT	genome size,default is 5000000
	--name   STR	cfg prefix name, default is "falcon"
	--cpu    INT	cpu number used for falcon, default is 10
	--mem    STR	memory used for falcon, default is 30G
LINE

my %opts;
GetOptions (\%opts, 
	'help!',
	'infile=s',
	'gsize=i',
	'name=s',
	'cpu=i',
	'mem=s'
);

die $INFO if defined($opts{help});
die "Please specify infile!!\n\n$INFO" unless defined($opts{infile});
die "$opts{infile} does not exist!!\n\n$INFO" unless -e $opts{infile};
#$opts{outdir} ||= "./Out";
$opts{cpu}    ||= 10;
$opts{mem}    ||= "30G";
$opts{name}   ||= "falcon";
$opts{gsize}  ||= 5000000;
$opts{infile} = File::Spec->rel2abs($opts{infile});

#make fc_run.cfg
open CFG,">$opts{name}.cfg";
print CFG "[General]
input_fofn=input.fofn
input_type=raw
pa_DBdust_option=
pa_fasta_filter_option=pass
target=assembly
skip_checks=False
LA4Falcon_preload=false
#### Data Partitioning
pa_DBsplit_option=-x500 -s50
ovlp_DBsplit_option=-x500 -s50

#### Repeat Masking
pa_HPCTANmask_option=
pa_REPmask_code=1,100;2,80;3,60

####Pre-assembly
genome_size=$opts{gsize}
seed_coverage=20
length_cutoff=-1
pa_HPCdaligner_option=-v -B128 -M24
pa_daligner_option=-e.7 -l1000 -k18 -h80 -w8 -s100
falcon_sense_option=--output-multi --min-idt 0.70 --min-cov 2 --max-n-read 800
falcon_sense_greedy=False
####Pread overlapping
ovlp_daligner_option=-e.96 -l2000 -k24 -h1024 -w6 -s100
ovlp_HPCdaligner_option=-v -B128 -M24

####Final Assembly
overlap_filtering_setting=--max-diff 100 --max-cov 300 --min-cov 2
fc_ovlp_to_graph_option=
length_cutoff_pr=10000


[job.defaults]
job_type=localshell
pwatcher_type = blocking
submit = sh -c \"\${JOB_SCRIPT}\" > \"\${JOB_STDOUT}\"
[job.step.fc]
NPROC=$opts{cpu}
MB=4000
njobs=8
[job.step.pda]
NPROC=$opts{cpu}
MB=32768
njobs=2
[job.step.da]
NPROC=$opts{cpu}
MB=32768
njobs=2
[job.step.la]
NPROC=$opts{cpu}
MB=32768
njobs=2
[job.step.cns]
NPROC=$opts{cpu}
MB=32768
njobs=2
[job.step.pla]
NPROC=$opts{cpu}
MB=32768
njobs=2
[job.step.asm]
NPROC=$opts{cpu}
MB=32768
njobs=1

";
close CFG;
`ls $opts{infile} > input.fofn`;


#!/usr/bin/perl -w
use strict;

### START OF INPUT ############################################################################################################################################

# This file converts from reads aligned against the genome (the Bowtie file) to counts of each of the four possible nucleotides at each position along the pre-miRNA sequence, for all the pre-miRNAs.  
# Usage: $ perl Analyze_mutation_modified.pl spc      spc.filtered.hit spc.hairpin.rnafold spc.hairpin.u2t.hit spc.mature.fa spc.mutation.txt
# Loc Des                                    $ARGV[0] $ARGV[1]         $ARGV[2]            $ARGV[3]            $ARGV[4]      $ARGV[5]

# Previous steps:
# 1) process_reads.pl
# 2) bowtie -n 1 -e 50 -a -m 1 --best --strata --phred33-quals --trim3 2 BOWTIE_INDEXES_PATHWAY FASTQ_FILE > BOWTIE_OUTPUT_FILE

### main input and output
my $miRBase_string = $ARGV[0]; # type the miRBase string for the organism in question. For example, for human write "hsa" and for mouse write "mmu" (this data can be found in miRBase)
my $file_filtered_hit = $ARGV[1];
my $file_hairpin_rnafold = $ARGV[2];
my $file_hairpin_hit = $ARGV[3];
my $file_mature = $ARGV[4];
my $file_main_output = $ARGV[5];

### specie files
open FILE_RNAfold, "<", $file_hairpin_rnafold or die $!; # RNAfold output file
open FILE_miRNA_Bowtie, "<", $file_hairpin_hit or die $!; # pre_miRNA aligned to genome

### important parameters
my $minimum_quality_factor=30; # the minimum quality factor that allows a mismatch to be considered
my $ignore_N_bases_at_the_end=2; # trim N bases at the end of the read

### general input files
open FILE_Bowtie_genome, "<", $file_filtered_hit or die $!; # bowtie output file: aligning to genome using bowtie 
open FILE_mature, "<", $file_mature or die $!; # mature miRNA file

### general output files
open FILE_out,">", $file_main_output or die $!; # The main output file
open FILE_out_tab,">", "miRNA_expression.txt" or die $!; # expression output file

### files for de-bugging
open FILE_out_matlab_mature,">", "temp_file1.txt" or die $!; 
open FILE_out_matlab_star,">", "temp_file2.txt" or die $!;
open FILE_reads_name,">", "temp_file3.txt" or die $!; 
open FILE_all_reads_name,">", "temp_file4.txt" or die $!;

### fine-tuning parameters
my $algorithm_used=4; # 1 for N-in-a-row, 2 for probalistic naive model, 3 for probalistic bootstrap model and 4 for binomial method (performed afterwards by the file Binomial_analysis.pl )
my $mismatch_in_mature_star_non_or_all=3; # 1 if the mutation is in mature miRNA, 2 if the mutation is in star miRNA and 0 if neither, 3 if all 
my $mismatch_in_loop_helix_or_all=2; # 0 if the mutation is in loop, 1 if the mutation is in helix, 2 if all; 
my $mismatch_in_a_row=3; # used only if $algorithm_used==1 
my $min_read_depth=5; 
my $probability_cutoff=1e-6; # not used if $algorithm_used==1 or $algorithm_used==4 
my $maximum_percent_noisy=100; # not used if $algorithm_used==1
my $cutoff_for_finding_mature_miRNA=0.5; # part of the maximum read depth
my $length_cutoff_for_comparing_mature_miRNA=0.9; # part of the miRBase definition
my $max_miRBase_after_end=2; # miRBase definition after the experimental end
my $number_of_high_quality_bases=19; # used only if $algorithm_used==3
my $p_value_cutoff=0.001; # used only if $algorithm_used==3
my $number_of_rand_runs=1000; # used only if $algorithm_used==3
my $length_of_very_short_miRNA=18;

### fixed parameters
my $number_of_fieldes=8;
my $location_of_read_name_field=1;
my $location_of_miRNA_field=3;
my $location_of_ref_bais_field=4;
my $location_of_read_field=5;
my $location_of_quality_field=6;
my $ascii_to_quality_convertor=33;
my $number_of_barcodes=1;
my $typical_length_of_miRNA=21;
my @mismatch_types=qw(AG AC AU GA GC GU UG UC UA CG CA CU);
my @types_of_letters=qw(A U G C);

### END OF INPUT #############################################################################################################################################

# parameter declaration 
my $string;
my @fields;
my $i;
my $j;
my $name_of_miRNA;
my $ref_bais;
my $quality_factor;
my %RNAfold_hash;
my @split_array;
my %mature_miRNA_hash;
my %star_miRNA_hash;
my %miRNA_seq_hash;
my $location_of_mature_miRNA;
my $location_of_star_miRNA;
my $length_of_mature_miRNA;
my $length_of_star_miRNA;
my $length_of_read;
my $identity_of_mutation; # 1 if the mutation is in mature miRNA, 2 if the mutation is in star miRNA and 0 if neither
my $string_fold;
my $name_of_read;
my $altered_name_of_miRNA;
my $location_of_mature_miRNA2;
my $length_of_mature_miRNA2;
my $condition1;
my $condition2;
my $condition3;
my $counting_mismatch;
my $counting_mismatch_total;
my $temp_index;
my $cut_name_of_miRNA;
my %global_hash_of_hashes;
my $read_seq;
my @all_the_miRNA;
my @all_the_miRNA_keys;
my @number_of_reads_per_barcode;
my @reads_per_miRNA_for_each_barcode;
my @reads_per_star_miRNA_for_each_barcode;
my $total_reads_per_miRNA;
my $temp_number;
my $number_of_read;
my $temp_read_number;
my $total_num_of_mismatch;
my $total_length_aligned;
my $error_rate;
my %num_of_letters_hash;
my $current_base;
my @other_types;
my $max_is;
my $different_base;
my %all_mutations_hash;
my $temp_quality;
my $test_flag;
my $altered_name;
my $string_mis;
my $string_org;
my $temp_string_mis;
my $temp_string_org;
my $N_in_a_row_location;
my $current_name;
my @depth_per_position;
my %hash_model_for_mismatch;
my $peak_location;
my $temp_depth_location;
my $possible_start_of_mature;
my $possible_end_of_mature;
my $cum_sum;
my $max_read;
my $no_problem_total;
my $problem_total;
my $no_reads_total;
my $flag_problem;
my $temp_location;
my $count_inside_mature_miRNA;
my %hash_mismatch_vs_position;
my @general_keys;
my %final_hash_mismatch_vs_position;
my %final_hash_model_for_mismatch;
my $counter;
my $odds_mismatch;
my @counts_array;
my @expected_array;
my $chi_sqr;
my $p_value;
my $exp_max;
my $pseudo_counts=1;
my $start_mature;
my $end_mature;
my $start_star;
my $end_star;
my $stop_looking;
my $mature_miRNA_representive_location;
my $star_miRNA_representive_location;
my $at_least_some_expression;
my %hash_of_surronding_mismatch;
my %hash_of_raw_mismatch;
my @all_keys;
my $barcode_number;
my $ignore_index;
my $cum_sum_after_end;
my $flag_after_end;
my $strand;
my $chr;
my $genomic_start;
my $genomic_end;
my $tmp_location;
my $tmp_position;
my %hash_genomic_locations_of_miRNA;
my $was;
my $now;
my $end_in_miRNA;
my $start_in_miRNA;
my @fields_now;
my $running_counter;
my @all_numbers;
my @all_numbers1;
my @all_numbers2;
my $temp_string;
            
### read the mature miRNA file into hash ####################################################################  
while ($string=<FILE_mature>)
{
    if ($string =~ /^>$miRBase_string-/)
    {
        $string =~ />(\S+)/;
        # $name_of_miRNA=$1;
        # $string=<FILE_mature>;
        # @split_array=split(/\s/,$string);
        # $mature_miRNA_hash{lc($name_of_miRNA)}=$split_array[0];
        if (substr($1,-1) eq "*")
        {
            $name_of_miRNA=substr($1,0,-1);
            $string=<FILE_mature>;
            @split_array=split(/\s/,$string);
            $star_miRNA_hash{lc($name_of_miRNA)}=$split_array[0];
        }
        else
        {
            $name_of_miRNA=$1;
            $string=<FILE_mature>;
            @split_array=split(/\s/,$string);
            $mature_miRNA_hash{lc($name_of_miRNA)}=$split_array[0];
        }
    }
}
close(FILE_mature);

### read the RNAfold file and locate the mature and star for each miRNA #####################################   
while ($string=<FILE_RNAfold>)
{
    if ($string =~ />(\S+)/)
    {
        $name_of_miRNA=$1;
        $string=<FILE_RNAfold>;
        if ($string =~ /(\S+)/)
        {
            # read the pre-miRNA sequence
            $miRNA_seq_hash{$name_of_miRNA}=$1;
            for ($i=0;$i<length($1);$i++)
            {
                $global_hash_of_hashes{$name_of_miRNA}{"miRNA_seq"}[$i]=substr($1,$i,1);
            }
            
            # find the location of mature miRNA on the pre-miRNA
            $location_of_mature_miRNA=0;
            $length_of_mature_miRNA=0;            
            $location_of_mature_miRNA2=0;
            $length_of_mature_miRNA2=0;
            if (defined $mature_miRNA_hash{$name_of_miRNA})
            {
                if (index($miRNA_seq_hash{$name_of_miRNA},$mature_miRNA_hash{$name_of_miRNA})!=(-1))
                {
                    $location_of_mature_miRNA=index($miRNA_seq_hash{$name_of_miRNA},$mature_miRNA_hash{$name_of_miRNA});
                    $length_of_mature_miRNA=length($mature_miRNA_hash{$name_of_miRNA});                  
                }
                else
                {
                    print "problem in mature miRNA: ",$name_of_miRNA,"\n";
                }
            }
            elsif (defined $mature_miRNA_hash{$name_of_miRNA."-5p"})
            {
                $altered_name_of_miRNA=$name_of_miRNA."-5p";
                if (index($miRNA_seq_hash{$name_of_miRNA},$mature_miRNA_hash{$altered_name_of_miRNA})!=(-1))
                {
                    $location_of_mature_miRNA=index($miRNA_seq_hash{$name_of_miRNA},$mature_miRNA_hash{$altered_name_of_miRNA});
                    $length_of_mature_miRNA=length($mature_miRNA_hash{$altered_name_of_miRNA});
                }
                else
                {
                    print "problem in mature miRNA: ",$name_of_miRNA,"\n";
                }                
                $altered_name_of_miRNA=$name_of_miRNA."-3p";
                if (defined $mature_miRNA_hash{$altered_name_of_miRNA})
                {
                    if (index($miRNA_seq_hash{$name_of_miRNA},$mature_miRNA_hash{$altered_name_of_miRNA})!=(-1))
                    {
                        $location_of_mature_miRNA2=index($miRNA_seq_hash{$name_of_miRNA},$mature_miRNA_hash{$altered_name_of_miRNA});
                        $length_of_mature_miRNA2=length($mature_miRNA_hash{$altered_name_of_miRNA});
                    }
                    else
                    {
                        print "problem in mature miRNA: ",$name_of_miRNA,"\n";
                    }
                }
            }            
            else
            {
                $temp_index='';
                $cut_name_of_miRNA='';
                $temp_index=rindex($name_of_miRNA,"-");
                $cut_name_of_miRNA=substr($name_of_miRNA,0,$temp_index);
                if (defined $mature_miRNA_hash{$cut_name_of_miRNA})
                {
                    if (index($miRNA_seq_hash{$name_of_miRNA},$mature_miRNA_hash{$cut_name_of_miRNA})!=(-1))
                    {
                        $location_of_mature_miRNA=index($miRNA_seq_hash{$name_of_miRNA},$mature_miRNA_hash{$cut_name_of_miRNA});
                        $length_of_mature_miRNA=length($mature_miRNA_hash{$cut_name_of_miRNA});
                    }
                    else
                    {
                        print "problem in mature miRNA: ",$name_of_miRNA,"\n";
                    }                    
                }
                elsif (defined $mature_miRNA_hash{$cut_name_of_miRNA."-5p"})
                {
                    $altered_name_of_miRNA=$cut_name_of_miRNA."-5p";
                    if (index($miRNA_seq_hash{$name_of_miRNA},$mature_miRNA_hash{$altered_name_of_miRNA})!=(-1))
                    {
                        $location_of_mature_miRNA=index($miRNA_seq_hash{$name_of_miRNA},$mature_miRNA_hash{$altered_name_of_miRNA});
                        $length_of_mature_miRNA=length($mature_miRNA_hash{$altered_name_of_miRNA});
                    }
                    else
                    {
                        print "problem in mature miRNA: ",$name_of_miRNA,"\n";
                    }
                    $altered_name_of_miRNA=$cut_name_of_miRNA."-3p";
                    if (defined $mature_miRNA_hash{$altered_name_of_miRNA})
                    {
                        if (index($miRNA_seq_hash{$name_of_miRNA},$mature_miRNA_hash{$altered_name_of_miRNA})!=(-1))
                        {
                            $location_of_mature_miRNA2=index($miRNA_seq_hash{$name_of_miRNA},$mature_miRNA_hash{$altered_name_of_miRNA});
                            $length_of_mature_miRNA2=length($mature_miRNA_hash{$altered_name_of_miRNA});
                        }
                        else
                        {
                            print "problem in mature miRNA: ",$name_of_miRNA,"\n";
                        }
                    }
                }
                else
                {
                    print "can't find the mature miRNA: ",$name_of_miRNA,"\n";
                }
            }
            
            ### find the location of star miRNA on the pre-miRNA ###
            $location_of_star_miRNA=0;
            $length_of_star_miRNA=0;
            if (defined $star_miRNA_hash{$name_of_miRNA})
            {
                if (index($miRNA_seq_hash{$name_of_miRNA},$star_miRNA_hash{$name_of_miRNA})!=(-1))
                {
                    $location_of_star_miRNA=index($miRNA_seq_hash{$name_of_miRNA},$star_miRNA_hash{$name_of_miRNA});
                    $length_of_star_miRNA=length($star_miRNA_hash{$name_of_miRNA});
                }
                else
                {
                    print "problem in star miRNA: ",$name_of_miRNA,"\n";
                }                
            }            
            else
            {
                $temp_index='';
                $cut_name_of_miRNA='';
                $temp_index=rindex($name_of_miRNA,"-");
                $cut_name_of_miRNA=substr($name_of_miRNA,0,$temp_index);
                if (defined $star_miRNA_hash{$cut_name_of_miRNA})
                {
                    if (index($miRNA_seq_hash{$name_of_miRNA},$star_miRNA_hash{$cut_name_of_miRNA})!=(-1))
                    {
                        $location_of_star_miRNA=index($miRNA_seq_hash{$name_of_miRNA},$star_miRNA_hash{$cut_name_of_miRNA});
                        $length_of_star_miRNA=length($star_miRNA_hash{$cut_name_of_miRNA});
                    }                   
                }
            }

            ### create none(0)/mature(1)/star(2) (NMS) assignment ###
            for ($i=0;$i<length($1);$i++)
            {
                $condition1=0;
                $condition2=0;
                $condition3=0;
                $condition1=($i<($length_of_mature_miRNA+$location_of_mature_miRNA))&&($i>=$location_of_mature_miRNA);
                $condition2=($i<($length_of_mature_miRNA2+$location_of_mature_miRNA2))&&($i>=$location_of_mature_miRNA2);
                $condition3=($i<($length_of_star_miRNA+$location_of_star_miRNA))&&($i>=$location_of_star_miRNA);
                if (($condition1)||($condition2))
                {
                    $global_hash_of_hashes{$name_of_miRNA}{"miRNA_NMS"}[$i]=1;
                }
                elsif ($condition3)
                {
                    $global_hash_of_hashes{$name_of_miRNA}{"miRNA_NMS"}[$i]=2;
                }
                else
                {
                    $global_hash_of_hashes{$name_of_miRNA}{"miRNA_NMS"}[$i]=0;
                }
            }
        }
        else
        {
            print "can't read the miRNA sequence!\n";
        }
        ### read the pre-miRNA second structure (SS) ###
        $string=<FILE_RNAfold>;
        @split_array=split(/\s/,$string);
        for ($i=0;$i<length($split_array[0]);$i++)
        {
            $global_hash_of_hashes{$name_of_miRNA}{"miRNA_SS"}[$i]=substr($split_array[0],$i,1);
        }
    }
}
close(FILE_RNAfold);

for ($j=0;$j<$number_of_barcodes;$j++)
{
    $number_of_reads_per_barcode[$j]=0;
}

print "end of pre-process part\n";

### read the miRNA bowtie file into hash ##########################################################################

while ($string=<FILE_miRNA_Bowtie>)
{
    # reading the data from the alignment file 
    @fields=();
    @fields=split(/\t/,$string);
    $name_of_miRNA=$fields[0];
    $strand=$fields[1];
    $chr=$fields[2];
    $genomic_start=$fields[3];
    $read_seq=$fields[4];
    if ($strand eq "+")
    {
        for ($i=0;$i<length($read_seq);$i++)
        {
            $tmp_location=$genomic_start+$i;
            $hash_genomic_locations_of_miRNA{$strand."\t".$chr."\t".$tmp_location}="+\t".$name_of_miRNA."\t".$i;    
        }
    }
    else
    {
        for ($i=0;$i<length($read_seq);$i++)
        {
            $tmp_location=$genomic_start+$i;
            $tmp_position=length($read_seq)-1-$i;
            $hash_genomic_locations_of_miRNA{$strand."\t".$chr."\t".$tmp_location}="+\t".$name_of_miRNA."\t".$tmp_position;    
        }        
    }
}

### convert the bowtie file  ##########################################################################
open FILE_Bowtie,">","temp_aligned_to_miRNA.txt" or die $!; # tab output file
while ($string=<FILE_Bowtie_genome>)
{
    # reading the data from the alignment file 
    @fields=();
    @fields=split(/\t/,$string);
    $strand=$fields[1];
    $chr=$fields[2];
    $genomic_start=$fields[3];
    $read_seq=$fields[4];
    $genomic_end=$genomic_start+length($read_seq)-1;
    if ((exists  $hash_genomic_locations_of_miRNA{$strand."\t".$chr."\t".$genomic_start})&&(exists  $hash_genomic_locations_of_miRNA{$strand."\t".$chr."\t".$genomic_end}))
    {
        if ($strand eq "+")
        {
            $was=$strand."\t".$chr."\t".$genomic_start;
            $now=$hash_genomic_locations_of_miRNA{$strand."\t".$chr."\t".$genomic_start};
            # print "$was\n";
            # print "$now\n";
            $string =~ s/\Q$was\E/$now/;
            my @tmp_fields = split(/\t/,$string);
            if ($tmp_fields[3] < 0)
            {
                next;
            }
            print FILE_Bowtie $string;
        }
        else
        {
            $was=$strand."\t".$chr."\t".$genomic_start;
            $now=$hash_genomic_locations_of_miRNA{$strand."\t".$chr."\t".$genomic_start};
            @fields_now=split(/\t/,$now);
            $end_in_miRNA=$fields_now[2];
            $start_in_miRNA=$end_in_miRNA-(length($read_seq)-1);
            $now =$fields_now[0]."\t".$fields_now[1]."\t".$start_in_miRNA;
            $string =~ s/\Q$was\E/$now/;
            $was=$fields[4];
            $now=reverse $was;
            $now=~tr/ACTG/TGAC/;
            $string =~ s/\Q$was\E/$now/;
            $was=$fields[5];
            $now=reverse $was;
            $string =~ s/\Q$was\E/$now/;
            if (exists $fields[7])
            {
                $was=$fields[7];
                $now=$was;
                $now=~tr/ACTG/TGAC/;
                $string =~ s/\Q$was\E/$now/;
            }
            my @tmp_fields = split(/\t/,$string);
            if ($tmp_fields[3] < 0)
            {
                next;
            }
            print FILE_Bowtie $string;
        }
    }
}
close FILE_Bowtie;
print "done with temp_aligned_to_miRNA.txt","\n";

### read the bowtie file into hash ##########################################################################
open FILE_Bowtie,"<","temp_aligned_to_miRNA.txt" or die $!; # tab output file
$number_of_read=0;
$total_num_of_mismatch=0;
$total_length_aligned=0;
while ($string=<FILE_Bowtie>)
{
    $number_of_read=$number_of_read+1;
    $temp_read_number=$number_of_read/100000;
    if ($temp_read_number !~ /\D/)
    {
        print "read number=",$number_of_read,"\n";
    }
    # reading the data from the alignment file 
    @fields=();
    $name_of_miRNA='';
    $ref_bais='';
    $name_of_read='';
    @fields=split(/\t/,$string);
    $name_of_miRNA=$fields[$location_of_miRNA_field-1];
    $name_of_read=$fields[$location_of_read_name_field-1];
    $ref_bais=$fields[$location_of_ref_bais_field-1]; # $fields[3] 17
    $read_seq=$fields[$location_of_read_field-1]; # $fields[4] "TAGCTTATCAGACTGATGTT"
    $read_seq =~ s/T/U/g;
    $length_of_read=length($read_seq); # 20
    $total_length_aligned=$total_length_aligned+$length_of_read;
    $quality_factor=$fields[$location_of_quality_field-1];
    # counting the number of reads per barcode
    $barcode_number=0;
    if (exists $number_of_reads_per_barcode[$barcode_number])
    {
        $number_of_reads_per_barcode[$barcode_number]=$number_of_reads_per_barcode[$barcode_number]+1;
    }
    else
    {
        $number_of_reads_per_barcode[$barcode_number]=1;
    }
    # assigning the data into the hash of hashes
    for ($i=$ref_bais;$i<($ref_bais+$length_of_read);$i++) # for ($i=17;$i<37;$++)
    {
        # the quality cutoff for mismatches       
        if (substr($read_seq,($i-$ref_bais),1) ne $global_hash_of_hashes{$name_of_miRNA}{"miRNA_seq"}[$i])
        {
            $temp_quality=ord(substr($quality_factor,($i-$ref_bais),1))-$ascii_to_quality_convertor;
            # create %hash_of_raw_mismatch regardless of quality factor
            if (defined $hash_of_raw_mismatch{$global_hash_of_hashes{$name_of_miRNA}{"miRNA_seq"}[$i].substr($read_seq,($i-$ref_bais),1)})
            {
                $hash_of_raw_mismatch{$global_hash_of_hashes{$name_of_miRNA}{"miRNA_seq"}[$i].substr($read_seq,($i-$ref_bais),1)}=$hash_of_raw_mismatch{$global_hash_of_hashes{$name_of_miRNA}{"miRNA_seq"}[$i].substr($read_seq,($i-$ref_bais),1)}+1;
            }
            else
            {
                $hash_of_raw_mismatch{$global_hash_of_hashes{$name_of_miRNA}{"miRNA_seq"}[$i].substr($read_seq,($i-$ref_bais),1)}=1;
            }
            # create %hash_of_surronding_mismatch regardless of quality factor
            if (defined $hash_of_surronding_mismatch{substr($read_seq,($i-$ref_bais-1),3)})
            {
                $hash_of_surronding_mismatch{substr($read_seq,($i-$ref_bais-1),3)}=$hash_of_surronding_mismatch{substr($read_seq,($i-$ref_bais-1),3)}+1;
            }
            else
            {
                $hash_of_surronding_mismatch{substr($read_seq,($i-$ref_bais-1),3)}=1;
            }
            if (($temp_quality>=$minimum_quality_factor)&&(($length_of_read-($i-$ref_bais+1))>=$ignore_N_bases_at_the_end))
            {
                $total_num_of_mismatch=$total_num_of_mismatch+1;
                $global_hash_of_hashes{$name_of_miRNA}{$name_of_read."_0_quality"}[$i]=$temp_quality;
                $global_hash_of_hashes{$name_of_miRNA}{$name_of_read."_0_seq"}[$i]=substr($read_seq,($i-$ref_bais),1);
            }           
        }
        else
        {
            $global_hash_of_hashes{$name_of_miRNA}{$name_of_read."_0_seq"}[$i]=substr($read_seq,($i-$ref_bais),1);
        }
    }            
}
close FILE_Bowtie;

$error_rate=$total_num_of_mismatch/$total_length_aligned;
print FILE_out "total number of mismatches:\t",$total_num_of_mismatch,"\n";
print FILE_out "total length aligned:\t",$total_length_aligned,"\n";
print FILE_out "error rate:\t",$error_rate,"\n";
print FILE_out "\n";
print FILE_out "raw mismatches distribution:\n";
@all_keys=keys %hash_of_raw_mismatch;
foreach (@all_keys)
{
    print FILE_out $_,"\t",$hash_of_raw_mismatch{$_},"\n";
}

#print FILE_out "raw surronding mismatches distribution:\n";
@all_keys=keys %hash_of_surronding_mismatch;
foreach (@all_keys)
{
    #print FILE_out $_,"\t",$hash_of_surronding_mismatch{$_},"\n";
}

### create matlab output files ##############################################################################
@all_the_miRNA=keys %global_hash_of_hashes;
print FILE_out_matlab_mature "matrix_of_reads_mature=[ \n";
print FILE_out_matlab_star "matrix_of_reads_star=[ \n";
for ($j=0;$j<=$#all_the_miRNA;$j++)
{
    @all_the_miRNA_keys=();
    $total_reads_per_miRNA=0;
    @all_the_miRNA_keys=keys %{$global_hash_of_hashes{$all_the_miRNA[$j]}};
    for ($i=0;$i<$number_of_barcodes;$i++)
    {
        $reads_per_miRNA_for_each_barcode[$i]=0;
        $reads_per_star_miRNA_for_each_barcode[$i]=0;
    }
    
    # find the start and end of each mature miRNA 
    $stop_looking=0;
    $start_mature=0;
    $end_mature=0;
    for ($i=0;$i<=$#{$global_hash_of_hashes{$all_the_miRNA[$j]}{"miRNA_seq"}};$i++)
    {
        if (($global_hash_of_hashes{$all_the_miRNA[$j]}{"miRNA_NMS"}[$i]==1)&&($stop_looking==0))
        {
            $start_mature=$i;
            $stop_looking=1;
        }
        if (($global_hash_of_hashes{$all_the_miRNA[$j]}{"miRNA_NMS"}[$i]==0)&&($stop_looking==1))
        {
            $end_mature=$i-1;
            $stop_looking=0;
        }
    }
    
    # find the start and end of each star miRNA 
    $stop_looking=0;
    $start_star=0;
    $end_star=0;
    for ($i=0;$i<=$#{$global_hash_of_hashes{$all_the_miRNA[$j]}{"miRNA_seq"}};$i++)
    {
        if (($global_hash_of_hashes{$all_the_miRNA[$j]}{"miRNA_NMS"}[$i]==2)&&($stop_looking==0))
        {
            $start_star=$i;
            $stop_looking=1;
        }
        if (($global_hash_of_hashes{$all_the_miRNA[$j]}{"miRNA_NMS"}[$i]==0)&&($stop_looking==1))
        {
            $end_star=$i-1;
            $stop_looking=0;
        }
    }

    # calc depth per position for each mature miRNA and return the maximum read position
    @depth_per_position=();
    $at_least_some_expression=0;
    for ($i=$start_mature;$i<=$end_mature;$i++)
    {
        $depth_per_position[$i]=0;
        foreach (@all_the_miRNA_keys)
        {
            if ($_ =~ /\_([0-9]{1,2})\_seq/)
            {               
                if (defined $global_hash_of_hashes{$all_the_miRNA[$j]}{$_}[$i])
                {
                    $depth_per_position[$i]=$depth_per_position[$i]+1;
                    $at_least_some_expression=1;
                }
            }
        }
    }
    $mature_miRNA_representive_location=0;
    if ($at_least_some_expression)
    {
        $mature_miRNA_representive_location=array_member_location(\@depth_per_position,max(@depth_per_position));
    }
    
    # calc depth per position for each star miRNA and return the maximum read position
    $star_miRNA_representive_location=0;
    if ($start_star!=0)
    {
        @depth_per_position=();
        $at_least_some_expression=0;
        for ($i=$start_star;$i<=$end_star;$i++)
        {
            $depth_per_position[$i]=0;
            foreach (@all_the_miRNA_keys)
            {
                if ($_ =~ /\_([0-9]{1,2})\_seq/)
                {               
                    if (defined $global_hash_of_hashes{$all_the_miRNA[$j]}{$_}[$i])
                    {
                        $depth_per_position[$i]=$depth_per_position[$i]+1;
                        $at_least_some_expression=1;
                    }
                }
            }
        }
        if ($at_least_some_expression)
        {
            $star_miRNA_representive_location=array_member_location(\@depth_per_position,max(@depth_per_position));
        }
    }
    
    # print FILE_out_matlab_mature $all_the_miRNA[$j],", mature location: ",$mature_miRNA_representive_location," locations: ",$start_mature,":",$end_mature,"\n";
    # print FILE_out_matlab_mature $all_the_miRNA[$j],", star location: ",$star_miRNA_representive_location," locations: ",$start_star,":",$end_star,"\n";

    # count the number of reads in mature miRNA 
    foreach (@all_the_miRNA_keys)
    {
        if ($_ =~ /(.+?\_([0-9]{1,2}))\_seq/)
        {
            print FILE_all_reads_name $1,"\t",$all_the_miRNA[$j],"\n";  
            if (defined $global_hash_of_hashes{$all_the_miRNA[$j]}{$_}[$mature_miRNA_representive_location])
            {		
		print FILE_reads_name $1,"\t",$all_the_miRNA[$j],"\n";                
		$reads_per_miRNA_for_each_barcode[$2]=$reads_per_miRNA_for_each_barcode[$2]+1;
                $total_reads_per_miRNA=$total_reads_per_miRNA+1;
            }
        }
    }

    # count the number of reads in star miRNA 
    foreach (@all_the_miRNA_keys)
    {
        if ($_ =~ /\_([0-9]{1,2})\_seq/)
        {
            if ($star_miRNA_representive_location!=0)
            {
                if (defined $global_hash_of_hashes{$all_the_miRNA[$j]}{$_}[$star_miRNA_representive_location])
                {
                    $reads_per_star_miRNA_for_each_barcode[$1]=$reads_per_star_miRNA_for_each_barcode[$1]+1;
                }
            }
        }
    }    
    
    print FILE_out_tab $all_the_miRNA[$j];
    for ($i=0;$i<$number_of_barcodes;$i++)
    {
        $temp_number=$reads_per_miRNA_for_each_barcode[$i];
        print FILE_out_matlab_mature $temp_number," ";
        print FILE_out_tab "\t",$temp_number;       
    }
    print FILE_out_tab "\n",$all_the_miRNA[$j]."*";
    for ($i=0;$i<$number_of_barcodes;$i++)
    {
        $temp_number=$reads_per_star_miRNA_for_each_barcode[$i];
        print FILE_out_matlab_star $temp_number," ";
        print FILE_out_tab "\t",$temp_number; 
    }    
    print FILE_out_matlab_mature "\n";
    print FILE_out_matlab_star "\n";
    print FILE_out_tab "\n";
}
print FILE_out_matlab_mature "]; \n";
print FILE_out_matlab_star "]; \n";

print FILE_out_matlab_mature "reads_per_barcode=[\n";
for ($j=0;$j<$number_of_barcodes;$j++)
{
    print FILE_out_matlab_mature $number_of_reads_per_barcode[$j],"\n";
}
print FILE_out_matlab_mature "]; \n";

print FILE_out_matlab_mature "names_of_miRNA={ \n";
for ($j=0;$j<=$#all_the_miRNA;$j++)
{
    print FILE_out_matlab_mature "'",$all_the_miRNA[$j],"' \n";
}
print FILE_out_matlab_mature "}; \n";


close(FILE_out_matlab_mature);
close(FILE_out_matlab_star);
close(FILE_out_tab);

### create model for mismatch - only for the mature miRNA ###################################################

$no_problem_total=0;
$problem_total=0;
$no_reads_total=0;
$running_counter=0;

for ($j=0;$j<=$#all_the_miRNA;$j++)
{
    @all_the_miRNA_keys=();
    $total_reads_per_miRNA=0;
    @all_the_miRNA_keys=keys %{$global_hash_of_hashes{$all_the_miRNA[$j]}};
    @depth_per_position=();
    
    # calc depth per position for each miRNA
    for ($i=0;$i<=$#{$global_hash_of_hashes{$all_the_miRNA[$j]}{"miRNA_seq"}};$i++)
    {
        $depth_per_position[$i]=0;
        foreach (@all_the_miRNA_keys)
        {
            if ($_ =~ /\_([0-9]{1,2})\_seq/)
            {               
                if (defined $global_hash_of_hashes{$all_the_miRNA[$j]}{$_}[$i])
                {
                    $depth_per_position[$i]=$depth_per_position[$i]+1;
                }
            }
        }
    }
    $max_read=max(@depth_per_position);
    if ($max_read>=$min_read_depth)
    {
        $peak_location=array_member_location(\@depth_per_position,$max_read);
        if ($peak_location eq 'problem')
        {
            print "problem in peak location!\n";
        }
        
        $temp_depth_location=$peak_location;
        while (($depth_per_position[$temp_depth_location]>=($cutoff_for_finding_mature_miRNA*$max_read))&&($temp_depth_location>0))
        {
            $temp_depth_location=$temp_depth_location-1;
        }
        $possible_start_of_mature=$temp_depth_location+1;
        
        $temp_depth_location=$peak_location;
        while (($depth_per_position[$temp_depth_location]>=($cutoff_for_finding_mature_miRNA*max(@depth_per_position)))&&($temp_depth_location<$#depth_per_position))
        {
            $temp_depth_location=$temp_depth_location+1;
        }
        $possible_end_of_mature=$temp_depth_location-1;
        
        $cum_sum=0;
        for ($i=$possible_start_of_mature;$i<=$possible_end_of_mature;$i++)
        {
            if ($global_hash_of_hashes{$all_the_miRNA[$j]}{"miRNA_NMS"}[$i]==1)
            {
                $cum_sum=$cum_sum+1;
            }
        }
        $cum_sum_after_end=0;
        $flag_after_end=1;
        while ($flag_after_end)
        {
            $i++;
            if (defined($global_hash_of_hashes{$all_the_miRNA[$j]}{"miRNA_NMS"}[$i]))
            {
                if ($global_hash_of_hashes{$all_the_miRNA[$j]}{"miRNA_NMS"}[$i]==1)
                {
                    $cum_sum_after_end++;
                }
                else
                {
                    $flag_after_end=0;
                }
            }
            else
            {
                $flag_after_end=0;                
            }
        }        
        if ((($cum_sum/($possible_end_of_mature-$possible_start_of_mature+1))<$length_cutoff_for_comparing_mature_miRNA)||($cum_sum_after_end>$max_miRBase_after_end))
        {
            #print FILE_out "problem in the miRbase mature definition in miRNA: ",$all_the_miRNA[$j],"\t",$possible_start_of_mature,":",$possible_end_of_mature,"\n";
            $problem_total=$problem_total+1;
            $flag_problem=1;
        }
        else
        {
            #print FILE_out "no problem in miRNA: ",$all_the_miRNA[$j],"\n";
            $no_problem_total=$no_problem_total+1;
            $flag_problem=0;
        }
        if (($possible_end_of_mature-$possible_start_of_mature+1)<$length_of_very_short_miRNA)
        {
            #print FILE_out "The expressed part of this miRNA is very short: ",$all_the_miRNA[$j],"\t",$possible_start_of_mature,":",$possible_end_of_mature,"\n";
        }
    }
    else
    {
        #print FILE_out "no reads for miRNA: ",$all_the_miRNA[$j],"\n";
        $no_reads_total=$no_reads_total+1;
    }

 
    # calc mismatch rate per position for all miRNA
    $count_inside_mature_miRNA=0;
    if ($max_read>=$min_read_depth)
    {
        for ($i=0;$i<=$#{$global_hash_of_hashes{$all_the_miRNA[$j]}{"miRNA_seq"}};$i++)
        {
            if ($global_hash_of_hashes{$all_the_miRNA[$j]}{"miRNA_NMS"}[$i]==1)
            {
                if (($depth_per_position[$i]>=$min_read_depth)&&($count_inside_mature_miRNA<=$number_of_high_quality_bases))
                {
                    $string_mis='';
                    $string_org=$global_hash_of_hashes{$all_the_miRNA[$j]}{"miRNA_seq"}[$i];
                    foreach (@all_the_miRNA_keys)
                    {    
                        if ($_ =~ /\_([0-9]{1,2})\_quality/)
                        {
                            $current_name=$_;
                            if (defined $global_hash_of_hashes{$all_the_miRNA[$j]}{$current_name}[$i])
                            {
                                $altered_name=$current_name;
                                if ($altered_name =~ s/_quality/_seq/)
                                {
                                    $string_mis=$string_mis.$global_hash_of_hashes{$all_the_miRNA[$j]}{$altered_name}[$i];
                                }
                            }
                        }
                    }
                    if ((length($string_mis)/$depth_per_position[$i]*100)<=$maximum_percent_noisy)
                    {
                        $num_of_letters_hash{"A"}=0;
                        $num_of_letters_hash{"U"}=0;
                        $num_of_letters_hash{"G"}=0;
                        $num_of_letters_hash{"C"}=0;
                        delete($num_of_letters_hash{$string_org});
                        @general_keys=keys %num_of_letters_hash;
                        foreach(@general_keys)
                        {
                            while ($string_mis=~/$_/g)
                            {
                               $num_of_letters_hash{$_}=$num_of_letters_hash{$_}+1;
                            }
                            if (defined($hash_model_for_mismatch{$string_org.$_}))
                            {
                                $temp_location=$#{$hash_model_for_mismatch{$string_org.$_}};
                                $hash_model_for_mismatch{$string_org.$_}[$temp_location+1]=$num_of_letters_hash{$_}/$depth_per_position[$i];
                            }
                            else
                            {
                                $hash_model_for_mismatch{$string_org.$_}[0]=$num_of_letters_hash{$_}/$depth_per_position[$i];
                            }
                        }
                        $cum_sum=0;
                        foreach(@general_keys)
                        {
                            $cum_sum=$cum_sum+$num_of_letters_hash{$_};
                        }
                        foreach(@general_keys)
                        {
                            $temp_string=$depth_per_position[$i]-$cum_sum+$num_of_letters_hash{$_};
                            # print FILE_out $running_counter,"\t",$count_inside_mature_miRNA,"\t",$string_org.$_,"\t",$num_of_letters_hash{$_},"\t",$temp_string,"\n";                            
                            $hash_mismatch_vs_position{$count_inside_mature_miRNA}{$string_org.$_}[$running_counter]=$num_of_letters_hash{$_}."_".$temp_string;
                        }          
                    }
                }
                $count_inside_mature_miRNA=$count_inside_mature_miRNA+1;
            }
            else
            {
                $count_inside_mature_miRNA=0;
            }
        }
    }
    $running_counter++;
}

print FILE_out "\n";
print FILE_out "Number of miRNAs with reads and mature position which fits the miRBase definition:\t",$no_problem_total,"\n";
print FILE_out "Number of miRNAs with reads and mature position which doesn't fit the miRBase definition:\t",$problem_total,"\n";
print FILE_out "Number of miRNAs without reads:\t",$no_reads_total,"\n";
print FILE_out "\n";

print "start analyze mismatches\n";

print FILE_out "rate of mismatch vs position inside mature miRNA:\n";
print FILE_out "Position\tMismatch_type\trate\n";
@general_keys=();
@general_keys=keys %hash_mismatch_vs_position;
@general_keys = sort {$a <=> $b} @general_keys;
for ($i=0;$i<=($typical_length_of_miRNA-$ignore_N_bases_at_the_end);$i++)
{
    foreach (@mismatch_types)
    {
        @all_numbers1=();
        @all_numbers2=();
        for ($j=0;$j<=$running_counter;$j++)
        {
            #print $i,"\t",$_,"\t",$j,"\n";
            if (exists $hash_mismatch_vs_position{$i}{$_}[$j])
            {
                if ($hash_mismatch_vs_position{$i}{$_}[$j] =~ /(\w+)\_(\w+)/)
                {
                    push(@all_numbers1,$1);
                    push(@all_numbers2,$2);
                }
                else
                {
                    print "something is wrong!\n";
                }
            }
        }
        if (($#all_numbers1>=0)&&($#all_numbers2>=0))
        {
            $final_hash_mismatch_vs_position{$i}{$_}=sum(@all_numbers1)/sum(@all_numbers2);
            print FILE_out $i,"\t",$_,"\t",$final_hash_mismatch_vs_position{$i}{$_},"\n";
        }
        else
        {
            $final_hash_mismatch_vs_position{$i}{$_}=0;
            print FILE_out $i,"\t",$_,"\t",$final_hash_mismatch_vs_position{$i}{$_},"\n";            
        }
    }
}
print "printing output file\n";
#print FILE_out "model for mismatch:\n";
@general_keys=();
@general_keys=keys %hash_model_for_mismatch;
foreach (@general_keys)
{
    $final_hash_model_for_mismatch{$_}=mean(@{$hash_model_for_mismatch{$_}});
#    print FILE_out $_,"\t",$final_hash_model_for_mismatch{$_},"\n";
}

### analyze mismatchs part ##################################################################################################################################

foreach(@mismatch_types)
{
    $all_mutations_hash{$_}=0;
}

### run for each known miRNA
print FILE_out "\n";
print FILE_out "miRNA_number\tmiRNA_name\n";
print FILE_out "Format:\t[1 2 0]=mature, star or neither\t[().]=secondary structure\t[A G C U]=genomic sequence\t[A G C U]=the read in this position\n";
print FILE_out "\n";
for ($j=0;$j<=$#all_the_miRNA;$j++)
#for ($j=340;$j<341;$j++)
{
    @all_the_miRNA_keys=();
    $total_reads_per_miRNA=0;
    @all_the_miRNA_keys=keys %{$global_hash_of_hashes{$all_the_miRNA[$j]}};
    print FILE_out $j,"\t",$all_the_miRNA[$j],"\n";
    
    # calc depth per position for each miRNA
    @depth_per_position=();
    for ($i=0;$i<=$#{$global_hash_of_hashes{$all_the_miRNA[$j]}{"miRNA_seq"}};$i++)
    {
        $depth_per_position[$i]=0;
        foreach (@all_the_miRNA_keys)
        {
            if ($_ =~ /\_([0-9]{1,2})\_seq/)
            {               
                if (defined $global_hash_of_hashes{$all_the_miRNA[$j]}{$_}[$i])
                {
                    $depth_per_position[$i]=$depth_per_position[$i]+1;
                }
            }
        }
    }
    
    ### for algorithm: N in a row scheme ####################################################################
    if (($mismatch_in_a_row>1)&&($algorithm_used==1))
    {
        foreach (@all_the_miRNA_keys)
        {
            if ($_ =~ /_quality/)
            {
                $current_name=$_;
                $string_mis='';
                $string_org='';
                for ($i=0;$i<=$#{$global_hash_of_hashes{$all_the_miRNA[$j]}{"miRNA_seq"}};$i++)
                {
                    if (defined $global_hash_of_hashes{$all_the_miRNA[$j]}{$current_name}[$i])
                    {
                        $altered_name=$current_name;
                        if ($altered_name =~ s/_quality/_seq/)
                        {
                            $string_mis=$string_mis.$global_hash_of_hashes{$all_the_miRNA[$j]}{$altered_name}[$i];
                            $string_org=$string_org.$global_hash_of_hashes{$all_the_miRNA[$j]}{"miRNA_seq"}[$i];
                        }
                        
                        # check other conditions
                        $test_flag=1;
                        if ($mismatch_in_mature_star_non_or_all!=3)
                        {
                            if ($mismatch_in_mature_star_non_or_all!=$global_hash_of_hashes{$all_the_miRNA[$j]}{"miRNA_NMS"}[$i])
                            {
                                $test_flag=0;
                            }
                        }
                        if ($mismatch_in_loop_helix_or_all!=2)
                        {
                            if (($mismatch_in_loop_helix_or_all==0)&&($global_hash_of_hashes{$all_the_miRNA[$j]}{"miRNA_SS"}[$i] ne "."))
                            {
                                $test_flag=0;
                            }
                            elsif (($mismatch_in_loop_helix_or_all==1)&&($global_hash_of_hashes{$all_the_miRNA[$j]}{"miRNA_SS"}[$i] eq "."))
                            {
                                $test_flag=0;
                            }
                        }
                        if ($depth_per_position[$i]<$min_read_depth)
                        {
                            $test_flag=0;    
                        }
                        for ($ignore_index=1;$ignore_index<=$ignore_N_bases_at_the_end;$ignore_index++)
                        {
                            if (defined $global_hash_of_hashes{$all_the_miRNA[$j]}{"miRNA_NMS"}[$i+$ignore_index])
                            {
                                if ($mismatch_in_mature_star_non_or_all!=$global_hash_of_hashes{$all_the_miRNA[$j]}{"miRNA_NMS"}[$i+$ignore_index])
                                {
                                    $test_flag=0;
                                }
                            }
                            else
                            {
                                $test_flag=0;
                            }
                        }
                    }
                }
                if (($string_mis =~ /((A{2,})|(U{2,})|(G{2,})|(C{2,}))/)&&($test_flag==1))
                {
                    if (length($1)>=$mismatch_in_a_row)
                    {
                        $N_in_a_row_location=index($string_mis,$1);
                        $temp_string_mis=substr($string_mis,$N_in_a_row_location,length($1));
                        $temp_string_org=substr($string_org,$N_in_a_row_location,length($1));
                        if ($temp_string_org =~ /^(A{2,}|U{2,}|G{2,}|C{2,})$/)
                        {
                            $current_base=substr($temp_string_org,0,1);
                            $different_base=substr($temp_string_mis,0,1);
                            $all_mutations_hash{$current_base.$different_base}=$all_mutations_hash{$current_base.$different_base}+1;
                            print FILE_out $string_org,"\t",$string_mis,"\n"; 
                        }
                    }
                }
            }
        }
    }
    
    ### for algorithm: statistical naive scheme #############################################
    elsif (($algorithm_used==2)||($algorithm_used==4))
    {
        for ($i=0;$i<=$#{$global_hash_of_hashes{$all_the_miRNA[$j]}{"miRNA_seq"}};$i++)
        {   
            $test_flag=1;
            if ($mismatch_in_mature_star_non_or_all!=3)
            {
                if ($mismatch_in_mature_star_non_or_all!=$global_hash_of_hashes{$all_the_miRNA[$j]}{"miRNA_NMS"}[$i])
                {
                    $test_flag=0;
                }
            }
            if ($mismatch_in_loop_helix_or_all!=2)
            {
               if (($mismatch_in_loop_helix_or_all==0)&&($global_hash_of_hashes{$all_the_miRNA[$j]}{"miRNA_SS"}[$i] ne "."))
               {
                   $test_flag=0;
               }
               elsif (($mismatch_in_loop_helix_or_all==1)&&($global_hash_of_hashes{$all_the_miRNA[$j]}{"miRNA_SS"}[$i] eq "."))
               {
                   $test_flag=0;
               }
            }   
            $string=$global_hash_of_hashes{$all_the_miRNA[$j]}{"miRNA_NMS"}[$i];
            $string=$string.$global_hash_of_hashes{$all_the_miRNA[$j]}{"miRNA_SS"}[$i];
            $string=$string.$global_hash_of_hashes{$all_the_miRNA[$j]}{"miRNA_seq"}[$i];
            foreach (@all_the_miRNA_keys)
            {
                if ($_ =~ /\_([0-9]{1,2})\_seq/)
                {
                    if (defined $global_hash_of_hashes{$all_the_miRNA[$j]}{$_}[$i])
                    {
                        $string=$string.$global_hash_of_hashes{$all_the_miRNA[$j]}{$_}[$i];
                    }
                }
            }
            print FILE_out substr($string,0,50)," length: ",length($string),"\n";
            if ((length($string)>=($min_read_depth+3))&&($test_flag==1))
            {
                $string=substr($string,3);
                $current_base=$global_hash_of_hashes{$all_the_miRNA[$j]}{"miRNA_seq"}[$i];
                if ($string !~ /^$current_base+$/)
                {
                    $num_of_letters_hash{"A"}=0;
                    $num_of_letters_hash{"U"}=0;
                    $num_of_letters_hash{"G"}=0;
                    $num_of_letters_hash{"C"}=0;          
                    foreach(@types_of_letters)
                    {
                        while ($string=~/$_/g)
                        {
                           $num_of_letters_hash{$_}=$num_of_letters_hash{$_}+1;
                        }               
                    }
                    print FILE_out $num_of_letters_hash{"A"},"\t",$num_of_letters_hash{"U"},"\t",$num_of_letters_hash{"G"},"\t",$num_of_letters_hash{"C"},"\n"; 
                    if (((length($string)-$num_of_letters_hash{$current_base})/length($string)*100)<=$maximum_percent_noisy)
                    {
                        delete($num_of_letters_hash{$current_base});
                        @other_types=keys(%num_of_letters_hash);
                        $max_is=0;
                        foreach(@other_types)
                        {
                            if ($num_of_letters_hash{$_}>$max_is)
                            {
                                $max_is=$num_of_letters_hash{$_};
                                $different_base=$_;
                            }                   
                        }
                        if ($algorithm_used==2)
                        {
                            if ((length($string)*($error_rate**($num_of_letters_hash{$different_base})))<$probability_cutoff)
                            {
                                $all_mutations_hash{$current_base.$different_base}=$all_mutations_hash{$current_base.$different_base}+1;
                                print FILE_out $current_base.$different_base," : ",length($string),",",$num_of_letters_hash{$different_base},"\n";
                                print FILE_out (length($string)*($error_rate**($num_of_letters_hash{$different_base}))),"\n";
                            }
                        }
                        if ($algorithm_used==4)
                        {
                            if ($num_of_letters_hash{$different_base}/length($string)>$error_rate)
                            {
                                $all_mutations_hash{$current_base.$different_base}=$all_mutations_hash{$current_base.$different_base}+1;
                                print FILE_out $current_base.$different_base," : ",length($string),",",$num_of_letters_hash{$different_base},"\n";
                            }
                        }                        
                    }
                }
            }
        }
    }
    
    #### for algorithm: statistical not-naive scheme #############################################
    else
    {
        # check only in mature miRNA
        $count_inside_mature_miRNA=-1;
        for ($i=0;$i<=$#{$global_hash_of_hashes{$all_the_miRNA[$j]}{"miRNA_seq"}};$i++)
        {   
            $test_flag=0;
            if (($global_hash_of_hashes{$all_the_miRNA[$j]}{"miRNA_NMS"}[$i]==1)&&($count_inside_mature_miRNA<$number_of_high_quality_bases))
            {
                $test_flag=1;
                $count_inside_mature_miRNA=$count_inside_mature_miRNA+1;
            }
            else
            {
                $count_inside_mature_miRNA=-1;
            }
            if ($mismatch_in_loop_helix_or_all!=2)
            {
               if (($mismatch_in_loop_helix_or_all==0)&&($global_hash_of_hashes{$all_the_miRNA[$j]}{"miRNA_SS"}[$i] ne "."))
               {
                   $test_flag=0;
               }
               elsif (($mismatch_in_loop_helix_or_all==1)&&($global_hash_of_hashes{$all_the_miRNA[$j]}{"miRNA_SS"}[$i] eq "."))
               {
                   $test_flag=0;
               }
            }   
            $string=$global_hash_of_hashes{$all_the_miRNA[$j]}{"miRNA_NMS"}[$i];
            $string=$string.$global_hash_of_hashes{$all_the_miRNA[$j]}{"miRNA_SS"}[$i];
            $string=$string.$global_hash_of_hashes{$all_the_miRNA[$j]}{"miRNA_seq"}[$i];
            foreach (@all_the_miRNA_keys)
            {
                if ($_ =~ /\_([0-9]{1,2})\_seq/)
                {
                    if (defined $global_hash_of_hashes{$all_the_miRNA[$j]}{$_}[$i])
                    {
                        $string=$string.$global_hash_of_hashes{$all_the_miRNA[$j]}{$_}[$i];
                    }
                }
            }
            print FILE_out substr($string,0,100)," length: ",length($string),"\n";
            if ((length($string)>=($min_read_depth+3))&&($test_flag==1))
            {
                $string=substr($string,3);
                $current_base=$global_hash_of_hashes{$all_the_miRNA[$j]}{"miRNA_seq"}[$i];
                if ($string !~ /^$current_base+$/)
                {
                    $num_of_letters_hash{"A"}=0;
                    $num_of_letters_hash{"U"}=0;
                    $num_of_letters_hash{"G"}=0;
                    $num_of_letters_hash{"C"}=0;          
                    foreach(@types_of_letters)
                    {
                        while ($string=~/$_/g)
                        {
                           $num_of_letters_hash{$_}=$num_of_letters_hash{$_}+1;
                        }               
                    }
                    print FILE_out $num_of_letters_hash{"A"},"\t",$num_of_letters_hash{"U"},"\t",$num_of_letters_hash{"G"},"\t",$num_of_letters_hash{"C"},"\n"; 
                    if (((length($string)-$num_of_letters_hash{$current_base})/length($string)*100)<=$maximum_percent_noisy)
                    {
                        delete($num_of_letters_hash{$current_base});
                        @other_types=keys(%num_of_letters_hash);
                        $max_is=0;
                        $odds_mismatch=();
                        @counts_array=();
                        foreach(@other_types)
                        {
                            if ($num_of_letters_hash{$_}>$max_is)
                            {
                                $max_is=$num_of_letters_hash{$_};
                                $different_base=$_;
                                $odds_mismatch=$final_hash_model_for_mismatch{$current_base.$_};
                                $counts_array[0]=$num_of_letters_hash{$_}+$pseudo_counts;
                            }                  
                        }
                        $counts_array[1]=length($string)-$counts_array[0]+$pseudo_counts;
                                                
                        @expected_array=();
                        $expected_array[0]=round(length($string)*$odds_mismatch)+$pseudo_counts;
                        $expected_array[1]=length($string)-$expected_array[0]+$pseudo_counts;
                        
                        #print FILE_out "odds: ",$odds_mismatch,"\n";
                        
                        for ($counter=0;$counter<=$#counts_array;$counter++)
                        {
                            #print FILE_out "counts: ",$counts_array[$counter],"\texpected: ",$expected_array[$counter],"\n";
                        }
                                           
                        if ($counts_array[0]>$expected_array[0])
                        {
                            $chi_sqr=0;
                            for ($counter=0;$counter<=$#counts_array;$counter++)
                            {
                                $chi_sqr=$chi_sqr+(($counts_array[$counter]-$expected_array[$counter])**2)/$expected_array[$counter];
                            }
                            #print FILE_out "chi= ",$chi_sqr,"\n";
                            
                            $p_value=calc_prob($odds_mismatch,length($string),$chi_sqr,$number_of_rand_runs);
                            #print FILE_out "P-value= ",$p_value,"\n";
                        
                            if (($p_value<=$p_value_cutoff)&&((length($string)*($error_rate**($num_of_letters_hash{$different_base})))<$probability_cutoff))
                            {
                                $all_mutations_hash{$current_base.$different_base}=$all_mutations_hash{$current_base.$different_base}+1;
                                print FILE_out $current_base.$different_base," : ",length($string),",",$num_of_letters_hash{$different_base},"\n";
                                if (($current_base.$different_base) eq "AG")
                                {
                                    #print FILE_out "editing ?!\n";
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    
}

### display results #########################################################################################
print FILE_out "\n";
print FILE_out "summary of mutations with higher rate compared to the expected sequencing error rate:\n";
foreach(@mismatch_types)
{
    print FILE_out $_," : ",$all_mutations_hash{$_},"\n";
}

close(FILE_out);

#############################################################################################################################################
sub round {
    my($number) = shift;
    return int($number + .5 * ($number <=> 0));
}

#############################################################################################################################################
sub max {
    my @numbers = @_;
    my $max;
    my $running=-1;
    for($running=0;$running<=$#numbers;$running++)
    {
        if (defined($numbers[$running]))
        {
            $max = $numbers[$running];
            last;
        }
    } 
    
    for($running=0;$running<=$#numbers;$running++)
    {
        if (defined($numbers[$running]))
        {
            if ($numbers[$running] > $max)
            {
                $max = $numbers[$running];
            }
        }
    }
    return $max;
}

#############################################################################################################################################
sub min {
    my @numbers = @_;
    my $min;
    my $running=-1;
    for($running=0;$running<=$#numbers;$running++)
    {
        if (defined($numbers[$running]))
        {
            $min = $numbers[$running];
            last;
        }
    }
    
    for($running=0;$running<=$#numbers;$running++)
    {
        if (defined($numbers[$running]))
        {
            if ($numbers[$running] < $min)
            {
                $min = $numbers[$running];
            }
        }
    }
    return $min;
}

#############################################################################################################################################
# input array_member_location(\@some_array,$some_member)
# this function returns the first occurence of the member
# 'problem' is returned if the member is not in the array
sub array_member_location {
    my @array=@{$_[0]};
    my $member=$_[1];
    my $sub_i;
    my $location='problem';
    for ($sub_i=0;$sub_i<=$#array;$sub_i++)
    {
        if (defined($array[$sub_i]))
        {
            if ($array[$sub_i] == $member)
            {
                $location = $sub_i;
                last
            }
        }
    }
    return $location;
}

#############################################################################################################################################
sub mean {
    my @numbers = @_;
    my $mean=0;
    foreach (@numbers)
    {
        $mean=$mean+$_;
    }
    $mean=$mean/($#numbers+1);
    return $mean;
}

#############################################################################################################################################
sub sum {
    my @numbers = @_;
    my $sum=0;
    foreach (@numbers)
    {
        $sum=$sum+$_;
    }
    return $sum;
}

#############################################################################################################################################
# input sub_odds, sub_number_of_reads and sub_chi ($sub_odds,$sub_number_of_reads,$sub_chi,$sub_number_of_repeats)
# this function returns the probabilty of getting chi higher or equal by chance
sub calc_prob {
    my $sub_odds=$_[0];
    my $sub_number_of_reads=$_[1];
    my $sub_chi=$_[2];
    my $sub_number_of_repeats=$_[3];
    my $sub_i;
    my $sub_j;
    my $sub_rand;
    my @sub_rand_count;
    my $sub_probabily;    
    my @sub_expected;    
    my $sub_temp_sum=0;
    my $sub_chi_rand;
    my $sub_pseudo_counts=1;
    
    $sub_expected[0]=round($sub_number_of_reads*$sub_odds)+$sub_pseudo_counts;
    $sub_expected[1]=$sub_number_of_reads-$sub_expected[0]+$sub_pseudo_counts;
    
    $sub_probabily=0;
    for ($sub_j=0;$sub_j<$sub_number_of_repeats;$sub_j++)
    {     
        $sub_rand_count[0]=$sub_pseudo_counts;
        $sub_rand_count[1]=$sub_pseudo_counts;       
        for ($sub_i=0;$sub_i<$sub_number_of_reads;$sub_i++)
        {
            $sub_rand=rand();
            if ($sub_rand<=$sub_odds)
            {
                $sub_rand_count[0]=$sub_rand_count[0]+1;
            }
            else
            {
                $sub_rand_count[1]=$sub_rand_count[1]+1;
            }
        }
        #print "rand counts: ",$sub_rand_count[0],"\t",$sub_rand_count[1],"\t",$sub_rand_count[2],"\t",$sub_rand_count[3],"\n";
        $sub_chi_rand=0;
        for ($sub_i=0;$sub_i<=$#sub_rand_count;$sub_i++)
        {
            $sub_chi_rand=$sub_chi_rand+(($sub_rand_count[$sub_i]-$sub_expected[$sub_i])**2)/$sub_expected[$sub_i];
        }
        #print "repeat: ",$sub_j,"\tchi_rand= ",$sub_chi_rand,"\n";
        if ($sub_chi_rand>=$sub_chi)
        {
            $sub_probabily=$sub_probabily+1; 
        }
    }
    
    $sub_probabily=$sub_probabily/$sub_number_of_repeats;
    return $sub_probabily;
}

#############################################################################################################################################

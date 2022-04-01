#!/usr/bin/perl -w
use strict;
use Math::CDF qw(:all);
my $prob = pnorm(1.96);

### START OF INPUT ############################################################################################################################################
# In this step binomial statistics is performed on the output file (or files) of analyze_mutation_modified.pl to separate sequencing errors from statistically significant modifications.  
# Usage: $ perl binomial_analysis_modified.pl INPUT_file OUTPUT_file
# Loc Des                            $ARGV[0] $ARGV[1]   $ARGV[2]

my $file_mutation_txt = $ARGV[0]; # the output file of analyze_mutation_modified.pl
my $file_hairpin_hit = $ARGV[1]; # the output file of pre-miRNA.fa against to genome by bowtie
my $file_binomial_xls = $ARGV[2]; # the output file of binomial_analysis_modified.pl

open FILE_miRNA_Bowtie, "<", $file_hairpin_hit or die $!; # miRNA aligned to genome

my $p_mismatch=0.001; # the sequencing error probability
my $max_position_to_check=19; # ignore last two bases
my $BH_or_Bonferroni="Bonferroni"; # multiple testing correction: a choice between "BH" (Benjamini-Hochberg correction) or "Bonferroni" (Bonferroni correction) 
my $P_value=1.00; # the P-value required after multiple testing correction (this is the FDR if Benjamini-Hochberg correction is chosen)

my $minimal_number_of_reads=5; # the bionomial test will be performed in all locations with at least this number of reads
my $maximal_percent_change=100; # maximum percent change allowed, setting this value to 100 indicate that even locations in which all the sequencing reads are different from the genomic sequence are considered  

my $first_N_bases=2; # If the mismatch is in the first N bases of the miRNA (as defined by $first_N_bases) AND the mismatch is in isomir with expression level which is significantly below the 'main' isomir, the mismatch is ignored (as we observed 5' adenylation and uridylation in low-abundance isomirs)    
my $fold_change_between_isomirs=10; # If the mismatch is in the first N bases of the miRNA AND the mismatch is in isomir with expression level which is significantly below the 'main' isomir (as defined by $fold_change_between_isomirs), the mismatch is ignored 

my $file_test="de-bugging.txt"; # de-bugging output file

### END OF INPUT #############################################################################################################################################

my $total_number_of_files=$#ARGV+1;
my $string;
my $global_flag;
my $number_of_expressed_miRNA;
my $number_of_multiple_tests;
my $counter;
my $current_p_value;
my @p_values;
my @sorted_p_values;
my $length_p_values;
my $i;
my $number_of_reads;
my %hash_all_counts;
my $miR_name;
my $current_base;
my @all_miRs;
my @all_locations;
my %hash_location_pvalue;
my $key;
my @array_expression;
my %hash_max_expression;
my @array_values;
my $flag_start_of_read;
my $flag_order_of_magnitude_difference;
my $global_running;
my $global_counter;
my %hash_convert_pre_miRNA_location_to_mature_location;
my $temp_miR_name;
my $previous_global_counter;
my %hash_5p_or_3p;
my %hash_all_letters;
my %temp_hash_all_letters;
my $temp_edited_levels;
my $temp_edited_letter;
my $total_expression_levels;
my $edited_expression_levels;
my %processed_p_values;
my @all_processed_p_values;
my %BH_adjusted_p_values;
my @part_of_all_processed_p_values;
my $multiple_tests_corrected_p_value;
my @fields;
my $name_of_miRNA;
my $strand;
my $chr;
my $genomic_start;
my $read_seq;
my $tmp_location;
my %hash_convert_position_in_miRNA_to_genomic_position;
my $tmp_position;

### read the miRNA bowtie file into hash 
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
            $hash_convert_position_in_miRNA_to_genomic_position{$name_of_miRNA."\t".($i+1)}=$chr."\t".$tmp_location."\t".($tmp_location+1)."\t".$strand;  
        }
    }
    else
    {
        for ($i=0;$i<length($read_seq);$i++)
        {
            $tmp_location=$genomic_start+$i;
            $tmp_position=length($read_seq)-1-$i;
            $hash_convert_position_in_miRNA_to_genomic_position{$name_of_miRNA."\t".($tmp_position+1)}=$chr."\t".$tmp_location."\t".($tmp_location+1)."\t".$strand;    
        }        
    }
}
close FILE_miRNA_Bowtie;

### hash of the RNA letters
$hash_all_letters{"A"}=1;
$hash_all_letters{"U"}=1;
$hash_all_letters{"G"}=1;
$hash_all_letters{"C"}=1;

### read the input files into hash
open FILE_test, ">", $file_test or die $!;

# for ($global_running=0;$global_running<$total_number_of_files;$global_running++)
# {
    # open FILE_R, "<", $file_names[$global_running] or die $!; 

open FILE_R, "<", $file_mutation_txt or die $!; 

$counter=0;
$global_counter=0;
$global_flag=1;
$string=<FILE_R>;

while ($global_flag==1)
{
    if ($string =~ /^[1,2].+?length:\s(\w+)/)
    {
        $counter++;
        $global_counter++;
        $hash_convert_pre_miRNA_location_to_mature_location{$miR_name}{$global_counter}=$counter;
        $number_of_reads=$1-3;
        $current_base=substr($string,2,1);
        if ($counter<=$max_position_to_check)
        {
            if (!($string=<FILE_R>))
            {
                $global_flag=0;
            }
            if ($string !~ /length:/)
            {
                if ($string =~ /([0-9]+)\t([0-9]+)\t([0-9]+)\t([0-9]+)/)
                {
                    if (exists $hash_all_counts{$miR_name}{$global_counter}{"RNA"})
                    {
                        $hash_all_counts{$miR_name}{$global_counter}{"A"}=$hash_all_counts{$miR_name}{$global_counter}{"A"}+$1;
                        $hash_all_counts{$miR_name}{$global_counter}{"U"}=$hash_all_counts{$miR_name}{$global_counter}{"U"}+$2;
                        $hash_all_counts{$miR_name}{$global_counter}{"G"}=$hash_all_counts{$miR_name}{$global_counter}{"G"}+$3;
                        $hash_all_counts{$miR_name}{$global_counter}{"C"}=$hash_all_counts{$miR_name}{$global_counter}{"C"}+$4; 
                    }
                    else
                    {
                        $hash_all_counts{$miR_name}{$global_counter}{"RNA"}=$current_base;
                        $hash_all_counts{$miR_name}{$global_counter}{"A"}=$1;
                        $hash_all_counts{$miR_name}{$global_counter}{"U"}=$2;
                        $hash_all_counts{$miR_name}{$global_counter}{"G"}=$3;
                        $hash_all_counts{$miR_name}{$global_counter}{"C"}=$4;
                    }
                }
                if (!($string=<FILE_R>))
                {
                    $global_flag=0;
                }
                if ($string =~ /^([A|U|G|C]{2})\s:/)
                {
                    if (!($string=<FILE_R>))
                    {
                        $global_flag=0;
                    }                    
                }
            }
            else
            {
                %temp_hash_all_letters=%hash_all_letters;
                delete($temp_hash_all_letters{$current_base});
                if (exists $hash_all_counts{$miR_name}{$global_counter}{"RNA"})
                {
                    $hash_all_counts{$miR_name}{$global_counter}{$current_base}=$hash_all_counts{$miR_name}{$global_counter}{$current_base}+$number_of_reads;  
                }
                else
                {
                    $hash_all_counts{$miR_name}{$global_counter}{"RNA"}=$current_base;
                    $hash_all_counts{$miR_name}{$global_counter}{$current_base}=$number_of_reads;
                    foreach (keys %temp_hash_all_letters)
                    {
                        $hash_all_counts{$miR_name}{$global_counter}{$_}=0;
                    }
                }
            }
        }
        else
        {
            if (!($string=<FILE_R>))
            {
                $global_flag=0;
            } 
        }
    }
    elsif ($string =~ /length:/)
    {
        $global_counter++;
        $counter=0;
        if (!($string=<FILE_R>))
        {
            $global_flag=0;
        }
    }
    else
    {
        if ($string =~ /[0-9]+\t([a-zA-Z]{3}.+?)\n/)
        {
            $global_counter=0;
            $counter=0;
            $miR_name=$1;
        }
        if (!($string=<FILE_R>))
        {
            $global_flag=0;
        }
    }
}
close FILE_R;
# }

### asign P-value to each location in the mature/star sequence 
@all_miRs=keys %hash_all_counts;
@array_expression=();
foreach(@all_miRs)
{
    $previous_global_counter=0;
    $miR_name=$_;
    $temp_miR_name=$miR_name;
    @all_locations=keys %{$hash_all_counts{$miR_name}};
    @all_locations = sort { $a <=> $b } @all_locations;
    foreach(@all_locations)
    {
        $global_counter=$_;
        if (exists $hash_all_counts{$miR_name}{$global_counter}{"RNA"})
        {
            $current_base=$hash_all_counts{$miR_name}{$global_counter}{"RNA"};
            %temp_hash_all_letters=%hash_all_letters;
            delete($temp_hash_all_letters{$current_base});
            $temp_edited_levels=0;
            $temp_edited_letter="A";
            foreach (keys %temp_hash_all_letters)
            {
                if ($hash_all_counts{$miR_name}{$global_counter}{$_}>=$temp_edited_levels)
                {
                    $temp_edited_levels=$hash_all_counts{$miR_name}{$global_counter}{$_};
                    $temp_edited_letter=$_;
                }
            } 
            $total_expression_levels=$hash_all_counts{$miR_name}{$global_counter}{"A"}+$hash_all_counts{$miR_name}{$global_counter}{"U"}+$hash_all_counts{$miR_name}{$global_counter}{"G"}+$hash_all_counts{$miR_name}{$global_counter}{"C"};
            $edited_expression_levels=$temp_edited_levels;
            print FILE_test $miR_name,"\t",$global_counter,"\t",$current_base.$temp_edited_letter,"\t",$edited_expression_levels,"\t",$total_expression_levels,"\n";
            
            if (($edited_expression_levels>0)&&($total_expression_levels>=$minimal_number_of_reads))
            {
                $current_p_value=pbinom(($total_expression_levels-$edited_expression_levels),$total_expression_levels,(1-$p_mismatch));
                print FILE_test $miR_name,"\t",$global_counter,"\tP-val=",$current_p_value,"\n";
            }
            else
            {
                $current_p_value=1;
            }
            push(@p_values,$current_p_value);
            $hash_location_pvalue{$miR_name."\t".$global_counter."\t".$current_base.$temp_edited_letter."\t".$edited_expression_levels."\t".$total_expression_levels}=$current_p_value;
            
            if ((($global_counter-$previous_global_counter)!=1)&&($previous_global_counter!=0))
            {
                $hash_max_expression{$temp_miR_name}=max(@array_expression);
                @array_expression=();
                $temp_miR_name=$miR_name."-3p";
                push(@array_expression,$total_expression_levels);
            }
            else
            {
                push(@array_expression,$total_expression_levels);
            }
            $previous_global_counter=$global_counter;
            $hash_5p_or_3p{$miR_name}{$global_counter}=$temp_miR_name;
        }
    }
    $hash_max_expression{$temp_miR_name}=max(@array_expression);
    @array_expression=();  
}
close FILE_test;

### calc BH-corrected P-values
@sorted_p_values = sort {$a <=> $b} @p_values;
$length_p_values=$#sorted_p_values+1;
#print "number of tests:\t",$length_p_values,"\n";

$i=0;
foreach $key (sort { $hash_location_pvalue{$a} <=> $hash_location_pvalue{$b}} keys %hash_location_pvalue)
{
    $i=$i+1;
    $processed_p_values{$key}=$hash_location_pvalue{$key}*$length_p_values/$i;
}
@all_processed_p_values=values %processed_p_values;
@all_processed_p_values = sort { $a <=> $b } @all_processed_p_values;
$i=-1;
foreach $key (sort { $processed_p_values{$a} <=> $processed_p_values{$b}} keys %processed_p_values)
{
    $i=$i+1;
    @part_of_all_processed_p_values=@all_processed_p_values[$i..$#all_processed_p_values];
    $BH_adjusted_p_values{$key}=min(@part_of_all_processed_p_values);
    if ($BH_adjusted_p_values{$key}>$P_value)
    {
        last;
    }
}
 
### display the P-values after multiple testing correction
### export result to a file, instead of stdout

open FILE_output, ">", $file_binomial_xls or die $!;

print FILE_output "#CHROM\tSTART\tEND\tSTRAND\tmiRNA_name\tLocation_inside_pre_miRNA\tMismatch_type\tNumber_of_reads_with_the_mismatch\tTotal_number_of_reads_in_this_position\tRaw_P_value\tBonferroni_P_value\tBH_P_value\n";
foreach $key (sort { $hash_location_pvalue{$a} <=> $hash_location_pvalue{$b}} keys %hash_location_pvalue)
{
    if ($BH_or_Bonferroni eq "BH")
    {
        $multiple_tests_corrected_p_value=$BH_adjusted_p_values{$key};
    }
    else
    {
        $multiple_tests_corrected_p_value=$hash_location_pvalue{$key}*$length_p_values;
    }
    if ($multiple_tests_corrected_p_value<=$P_value)
    {
        @array_values=split("\t",$key);           
        if ($hash_max_expression{$hash_5p_or_3p{$array_values[0]}{$array_values[1]}}>($array_values[4]*$fold_change_between_isomirs))
        {
            $flag_order_of_magnitude_difference=1;
        }
        else
        {
            $flag_order_of_magnitude_difference=0;
        }
        if ($hash_convert_pre_miRNA_location_to_mature_location{$array_values[0]}{$array_values[1]}<=$first_N_bases)
        {
            $flag_start_of_read=1;
        }
        else
        {
            $flag_start_of_read=0;
        }
        if (($flag_start_of_read==0)||($flag_order_of_magnitude_difference==0))
        {
            print FILE_output $hash_convert_position_in_miRNA_to_genomic_position{$array_values[0]."\t".$array_values[1]},"\t",$key,"\t",$hash_location_pvalue{$key},"\t",$hash_location_pvalue{$key}*$length_p_values,"\t",$BH_adjusted_p_values{$key},"\n";
        }
    }
    else
    {
        last;
    }
}

close FILE_output;

#############################################################################################################################################
sub max {
    my @numbers = @_;
    my $max;
    my $running;
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
    my $running;
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

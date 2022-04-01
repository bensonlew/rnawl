#!/usr/bin/perl
#contact: ming.lei@majorbio.com
#create:2015/8/10
#upadte:2015/8/10

use strict;
use warnings;
use Getopt::Long;
use Parallel::ForkManager;
use File::Basename;
use Cwd;
use Pod::Usage;

my $help          = q{};
my $man           = q{};
my $scaffold_file = q{};
my $proteins_file = q{};
my $proteins_dir  = q{};
my $genblast_dir = q{};
my $genewise_dir  = q{};
my %scaffolds;

my $genewise_bin = q{}; #guanqing.zou 20180530
my $genblast_bin = q{}; #guanqing.zou 20180530



my $current_path = cwd();
my $arg_num      = scalar @ARGV;

&GetOptions(
    'help|h!'         => \$help,
    'man!'            => \$man,
    'scaffold_file=s' => \$scaffold_file,
    'proteins_file=s' => \$proteins_file,
    'proteins_dir=s'  => \$proteins_dir,
    'genblast_dir=s' => \$genblast_dir,
    'genewise_dir=s'  => \$genewise_dir,

    'genewise_bin=s' => \$genewise_bin,  #guanqing.zou 20180530
    'genblast_bin=s' => \$genblast_bin   #guanqing.zou 20180530
);

if ($help) {
    pod2usage(1);
}

if ($man) {
    pod2usage( -verbose => 2 );
}

if ( $arg_num == 0 ) {
    pod2usage(1);
    exit;
}

############################Main program###########################
#
main: {
    
    mkdir $proteins_dir;
    mkdir $genblast_dir;

    open INONE, "<", $proteins_file;
    $/ = ">";
    my @proteins = <INONE>;
    $/ = "\n";
    close INONE;

    foreach my $protein (@proteins) {
        if ( $protein =~ /^((\S+).*?)\n(.*)/msg ) {
            my $header   = $1;
            my $prefix   = $2;
            my $sequence = $3;
            $sequence =~ s/>//;
            $header =~ s/[\.|\|]/_/g;
            $prefix =~ s/[\.|\|]/_/g;
            open OUT, ">", "$proteins_dir/$prefix.fasta";
            print OUT ">", $header, "\n", $sequence;
            close OUT;
        }
    }
    my @protein_files = glob "$proteins_dir/*fasta";

    `cp $genblast_bin/alignscore.txt ./`;  #zouguanqing 20180530
    #`ln -s /mnt/ilustre/users/sanger-dev/app/bioinfo/annotation/my_interproscan/interproscan-5.27-66.0/bin/blast/2.2.24/bin/formatdb ./`;
    `ln -s $genblast_bin/formatdb ./`;
    `ln -s $genblast_bin/blastall ./`;
    `$genblast_bin/formatdb -i $scaffold_file`;  #zouguanqing 20180530
    &parallel_genblastg( \@protein_files, $scaffold_file, $genblast_dir );

    mkdir $genewise_dir;

    open INTWO, "<", $scaffold_file;
    while ( my $line = <INTWO> ) {
        if ( $line =~ /^(\S+)/ ) {
            my $key = $1;
            $key =~ s/>//;
            $/ = ">";
            my $seq = <INTWO>;
            $/ = "\n";
            $seq =~ s/>|\s//g;
            $scaffolds{$key} = $seq;
        }
    }
    close INTWO;

    my @genblast_files = glob "$genblast_dir/*0";
	my %start_end;
    &parallel_genewise( \@genblast_files, $genewise_dir, \%scaffolds, $proteins_dir, \%start_end);

    open OUT, ">", "pseudogene_information.xls";
    print OUT "Gene ID\tScaffold ID\tBegin\tEnd\tFrameshifts\tPremature stop codons\n";
    close OUT;

    my @genewise_out_files = glob "$genewise_dir/*genewise.out";
    foreach my $genewise_out_file (@genewise_out_files) {
        &deal_genewise_out( $genewise_out_file, $proteins_dir ,\%start_end);
    }
}

#############################subroutine##############################
sub parallel_genblastg {
    my ( $protein_files_ref, $scaffold_file, $genblast_dir ) = @_;
    my $pm = new Parallel::ForkManager(30);
    foreach my $protein_file (@$protein_files_ref) {
        $pm->start and next;
        &genblastg( $protein_file, $scaffold_file, $genblast_dir );
        $pm->finish;
    }
    $pm->wait_all_children;
}

sub genblastg {
    my ( $protein_file, $scaffold_file, $genblast_dir ) = @_;
    my $output_prefix = $1, if ( $protein_file =~ /\/(\S+)\.fasta/ );
    my $cmd =
"$genblast_bin/genblast_v138_linux_x86_64  -p  genblastg  -q  $protein_file   -t  $scaffold_file  -e  1e-5 -g  T  -f  F  -a 0.5 -d 100000  -r  10  -c  0.5  -s  0  -i  15  -x  20 -n 20 -v 2 -h 0 -j 3  -norepair -gff -cdna -pro -o  $genblast_dir/$output_prefix";
    print $cmd, "\n";
    system $cmd;
}

sub parallel_genewise {
    my (
        $genblast_files_array_ref, $genewise_dir,
        $scaffolds_hash_ref,        $proteins_dir,
		$start_end_hash_ref,
    ) = @_;

    my $pm = new Parallel::ForkManager(30);
    foreach my $genblast_file (@$genblast_files_array_ref) {
        $pm->start and next;
        &genewise(
            $genblast_file,     $genewise_dir,
            $scaffolds_hash_ref, $proteins_dir,
			$start_end_hash_ref,
        );
        $pm->finish;
    }
    $pm->wait_all_children;
}

sub genewise {
    my ( $genblast_file, $genewise_dir, $scaffold_seq_hash_ref, $proteins_dir ,$start_end_hash_ref)
      = @_;
    my $gene_name = $1, if ( $genblast_file =~ /\/(.*?)\_1.1c_2.3/ );
    my %ranks;

	my $ex_start;
	my $ex_end;

    open IN, "<", $genblast_file;
    while ( my $line = <IN> ) {
        if ( $line =~ /^.*?(rank:1)$/ ) {
            chomp $line;
            $ranks{$1} = $line;
        }
    }
    close IN;

    foreach my $key ( keys %ranks ) {
        if ( $ranks{$key} =~
            /^(.*?)\|((.*?)\:(\d+)\.\.(\d+)\|(\S)\|.*?rank:(\d+))/ )
        {
            my $gene_id      = $1;
            my $header       = $2;
            my $scaffold_id  = $3;
            my $start        = $4;
            my $end          = $5;
            my $orientation  = $6;
            my $rank         = $7;
            my $scaffold_seq = $$scaffold_seq_hash_ref{$scaffold_id};
            my $scaffold_len = length $scaffold_seq;
            my $dna_seq;
          
            $header =~ s/\|/_/g;
            open OUT, ">", "$genewise_dir/$gene_name\_$rank.fa";
            if ( $start < 20000 ) {
                if ( $scaffold_len - $end <= 20000 ) {
                    $dna_seq = $scaffold_seq;
					$ex_start = 1;
					$ex_end = $scaffold_len;
                }
                else {
                    $dna_seq = substr( $scaffold_seq, 0, $end + 20000 );
					$ex_start = 1;
					$ex_end = $end+20000;
                }
            }
            else {
                if ( $scaffold_len - $end <= 20000 ) {
                    $dna_seq = substr( $scaffold_seq, $start - 20001 );
					$ex_start = $start - 20001;
					$ex_end = $scaffold_len;
                }
                else {
                    $dna_seq = substr(
                        $scaffold_seq,
                        $start - 20001,
                        $end - $start + 40000
                    );
					$ex_start = $start - 20001;
					$ex_end = $start+20000;
                }
            }

            if ( $orientation =~ /\+/ ) {
                print OUT ">", $scaffold_id, "\t", $header, "\n";
                print OUT $dna_seq, "\n";
            }
            elsif ( $orientation =~ /\-/ ) {
                $dna_seq = reverse($dna_seq);
                $dna_seq =~ tr/atcgATCG/TAGCtagc/;
                print OUT ">", $scaffold_id, "\t", $header, "\n";
                print OUT $dna_seq;
				$ex_start = $ex_end;
            }
		    $$start_end_hash_ref{$gene_name} = $ex_start;
`$genewise_bin/genewise  -both  $proteins_dir/$gene_name.fasta  $genewise_dir/$gene_name\_$rank.fa  -alg 333  -silent   -pseudo  -divide "//" -trans  -cdna   >>  $genewise_dir/$gene_name\_$rank\_genewise.out`;
        }
    }
}

sub deal_genewise_out {
    my ( $genewise_out_file, $proteins_dir ,$start_end_hash_ref) = @_;
    print $genewise_out_file, "\n";
    my $prefix  = $1, if ( $genewise_out_file =~ /\/(.*?)\_\d+\_genewise.out/ );
    my $pro_seq = `grep -v ">"  $proteins_dir/$prefix.fasta`;
    my $flag    = 1;
    my %genewise_outs;
    my $target_seq = q{};
    chomp $pro_seq;

    open IN, "<", $genewise_out_file;
    $/ = "genewise output";
    while ( my $lines = <IN> ) {
        if ( $lines =~ /Score\s+(\d+\.\d+)\s+bits/ ) {
            $genewise_outs{$lines} = $1;
        }
        elsif ( $lines =~ /Target\s+Sequence\s+(.*?)\n/msg ) {
            $target_seq = $1;
        }
    }
    $/ = "\n";
    close IN;

    my @keys =
      sort { $genewise_outs{$b} <=> $genewise_outs{$a} } keys %genewise_outs;
    my $best_genewise_out = $keys[0];
    print $genewise_outs{$keys[0]},"\n";
    print "$genewise_out_file.pseudogene.fasta\n";
    open OUT,    ">",  "$genewise_out_file.pseudogene.fasta";
    open OUTONE, ">>", "pseudogene_information.xls";
    print OUT ">", $prefix, "\n", $pro_seq, "\n";

    my $frameshift = "NO";
    my $stopsite   = "NO";
    my $translation;

    my @blocks                 = split( "//", $best_genewise_out );
    my $align_information      = $blocks[0];
    my $pep_seq_information    = $blocks[1];
    my $coding_seq_information = $blocks[2];

    if ( $align_information =~ /\!/ ) {
        $frameshift = "YES";
    }

    if ( $align_information =~ /\*/ ) {
        $stopsite = "YES";
    }

    $coding_seq_information =~ s/\s//g;
    my $header = q{};
    if ( $coding_seq_information =~ />(.*?.)sp(.*?)$/msg ) {
        $header = $1;
        my $dna_seq = $2;
        $dna_seq =~ s/\s|\///g;
        $translation = &dnapeptide($dna_seq);
    }

    my $scaffold_id = q{};
    my $start       = q{};
    my $end         = q{};

	my $ex_start = $1, if($$start_end_hash_ref{$prefix} =~ /(\S+)/);
    if ( $header =~ /(.*?)\.\[(\d+)\:(\d+)\]/ ) {
        $scaffold_id = $1;
        $start       = $ex_start+$2;
        $end         = $ex_start+$3;
    }

    if ( $pep_seq_information =~ /Gene 0 is a pseudo gene/ ) {
        print OUT ">", $header, "\t", length($translation), "\n";
        print OUTONE $prefix, "\t", $scaffold_id, "\t", $start, "\t", $end,
          "\t", $frameshift, "\t", $stopsite, "\n";
        print OUT $translation, "\n";
    }
    else {
        $flag = 0;
        `rm   $genewise_out_file.pseudogene.fasta`;
    }
    close OUT;
    close OUTONE;
}

sub dnapeptide {
    my ($dna) = @_;
    use strict;
    use warnings;
    my $protein = '';
    for ( my $i = 0 ; $i < ( length($dna) - 2 ) ; $i += 3 ) {
        $protein .= &codonaa( substr( $dna, $i, 3 ) );
    }
    return $protein;
}

sub codonaa {
    my ($codon) = @_;
    if    ( $codon =~ /TCA/i ) { return 'S' }    # Serine
    elsif ( $codon =~ /TCC/i ) { return 'S' }    # Serine
    elsif ( $codon =~ /TCG/i ) { return 'S' }    # Serine
    elsif ( $codon =~ /TCT/i ) { return 'S' }    # Serine
    elsif ( $codon =~ /TCN/i ) { return 'S' }    # Serine
    elsif ( $codon =~ /TTC/i ) { return 'F' }    # Phenylalanine
    elsif ( $codon =~ /TTT/i ) { return 'F' }    # Phenylalanine
    elsif ( $codon =~ /TTA/i ) { return 'L' }    # Leucine
    elsif ( $codon =~ /TTG/i ) { return 'L' }    # Leucine
    elsif ( $codon =~ /TTN/i ) { return 'X' }    # Leucine
    elsif ( $codon =~ /TAC/i ) { return 'Y' }    # Tyrosine
    elsif ( $codon =~ /TAT/i ) { return 'Y' }    # Tyrosine
    elsif ( $codon =~ /TAN/i ) { return 'X' }    # Tyrosine
    elsif ( $codon =~ /TAA/i ) { return '*' }    # Stop
    elsif ( $codon =~ /TAG/i ) { return '*' }    # Stop
    elsif ( $codon =~ /TGC/i ) { return 'C' }    # Cysteine
    elsif ( $codon =~ /TGT/i ) { return 'C' }    # Cysteine
    elsif ( $codon =~ /TGN/i ) { return 'X' }    # Cysteine
    elsif ( $codon =~ /TGA/i ) { return '*' }    # Stop
    elsif ( $codon =~ /TGG/i ) { return 'W' }    # Tryptophan
    elsif ( $codon =~ /CTA/i ) { return 'L' }    # Leucine
    elsif ( $codon =~ /CTC/i ) { return 'L' }    # Leucine
    elsif ( $codon =~ /CTG/i ) { return 'L' }    # Leucine
    elsif ( $codon =~ /CTT/i ) { return 'L' }    # Leucine
    elsif ( $codon =~ /CTN/i ) { return 'L' }    # Leucine
    elsif ( $codon =~ /CCA/i ) { return 'P' }    # Proline
    elsif ( $codon =~ /CCC/i ) { return 'P' }    # Proline
    elsif ( $codon =~ /CCG/i ) { return 'P' }    # Proline
    elsif ( $codon =~ /CCT/i ) { return 'P' }    # Proline
    elsif ( $codon =~ /CCN/i ) { return 'P' }    # Proline
    elsif ( $codon =~ /CAC/i ) { return 'H' }    # Histidine
    elsif ( $codon =~ /CAT/i ) { return 'H' }    # Histidine
    elsif ( $codon =~ /CAA/i ) { return 'Q' }    # Glutamine
    elsif ( $codon =~ /CAG/i ) { return 'Q' }    # Glutamine
    elsif ( $codon =~ /CAN/i ) { return 'X' }    # Glutamine
    elsif ( $codon =~ /CGA/i ) { return 'R' }    # Arginine
    elsif ( $codon =~ /CGC/i ) { return 'R' }    # Arginine
    elsif ( $codon =~ /CGG/i ) { return 'R' }    # Arginine
    elsif ( $codon =~ /CGT/i ) { return 'R' }    # Arginine
    elsif ( $codon =~ /CGN/i ) { return 'R' }    # Arginine
    elsif ( $codon =~ /ATA/i ) { return 'I' }    # Isoleucine
    elsif ( $codon =~ /ATC/i ) { return 'I' }    # Isoleucine
    elsif ( $codon =~ /ATT/i ) { return 'I' }    # Isoleucine
    elsif ( $codon =~ /ATN/i ) { return 'X' }    # Isoleucine
    elsif ( $codon =~ /ATG/i ) { return 'M' }    # Methionine
    elsif ( $codon =~ /ACA/i ) { return 'T' }    # Threonine
    elsif ( $codon =~ /ACC/i ) { return 'T' }    # Threonine
    elsif ( $codon =~ /ACG/i ) { return 'T' }    # Threonine
    elsif ( $codon =~ /ACT/i ) { return 'T' }    # Threonine
    elsif ( $codon =~ /ACN/i ) { return 'T' }    # Threonine
    elsif ( $codon =~ /AAC/i ) { return 'N' }    # Asparagine
    elsif ( $codon =~ /AAT/i ) { return 'N' }    # Asparagine
    elsif ( $codon =~ /AAA/i ) { return 'K' }    # Lysine
    elsif ( $codon =~ /AAG/i ) { return 'K' }    # Lysine
    elsif ( $codon =~ /AAN/i ) { return 'X' }
    elsif ( $codon =~ /AGC/i ) { return 'S' }    # Serine
    elsif ( $codon =~ /AGT/i ) { return 'S' }    # Serine
    elsif ( $codon =~ /AGA/i ) { return 'R' }    # Arginine
    elsif ( $codon =~ /AGG/i ) { return 'R' }    # Arginine
    elsif ( $codon =~ /AGN/i ) { return 'X' }    # Arginine
    elsif ( $codon =~ /GTA/i ) { return 'V' }    # Valine
    elsif ( $codon =~ /GTC/i ) { return 'V' }    # Valine
    elsif ( $codon =~ /GTG/i ) { return 'V' }    # Valine
    elsif ( $codon =~ /GTT/i ) { return 'V' }    # Valine
    elsif ( $codon =~ /GTN/i ) { return 'V' }    # Valine
    elsif ( $codon =~ /GCA/i ) { return 'A' }    # Alanine
    elsif ( $codon =~ /GCC/i ) { return 'A' }    # Alanine
    elsif ( $codon =~ /GCG/i ) { return 'A' }    # Alanine
    elsif ( $codon =~ /GCT/i ) { return 'A' }    # Alanine
    elsif ( $codon =~ /GCN/i ) { return 'A' }    # Alanine
    elsif ( $codon =~ /GAC/i ) { return 'D' }    # Aspartic Acid
    elsif ( $codon =~ /GAT/i ) { return 'D' }    # Aspartic Acid
    elsif ( $codon =~ /GAA/i ) { return 'E' }    # Glutamic Acid
    elsif ( $codon =~ /GAG/i ) { return 'E' }    # Glutamic Acid
    elsif ( $codon =~ /GAN/i ) { return 'X' }
    elsif ( $codon =~ /GGA/i ) { return 'G' }    # Glycine
    elsif ( $codon =~ /GGC/i ) { return 'G' }    # Glycine
    elsif ( $codon =~ /GGG/i ) { return 'G' }    # Glycine
    elsif ( $codon =~ /GGT/i ) { return 'G' }    # Glycine
    elsif ( $codon =~ /GGN/i ) { return 'G' }    # Glycine
    elsif ( $codon =~ /.N./i ) { return 'X' }
    elsif ( $codon =~ /N../i ) { return 'X' }
    elsif ( $codon =~ /.NN/i ) { return 'X' }
    elsif ( $codon =~ /NN./i ) { return 'X' }
    elsif ( $codon =~ /N.N/i ) { return 'X' }
    elsif ( $codon =~ /NNN/i ) { return 'X' }
    else {
        print STDERR "Bad codon \"$codon\"!!\n";
        exit;
    }
}
#`rm alignscore.txt formatdb blastall`;

__END__

=head1 NAME

Runall_for_pseudogenes_identification.pl - identify pseudogenes of newly sequenced genome

=head1 SYNOPSIS

Runall_for_pseudogenes_identification.pl [options]

  Options:
     -help		*Help command
     -man		*Explanation command
     -scaffold_file	*Genome sequences file  
     -proteins_file	*Protein sequences file
     -proteins_dir	*Directory containing protein files 	
     -genblast_dir     *Directory containing results of genblast analyses 
     -genewise_dir      *Directory containing results of genewise analyses
     -genblast_bin      *genblast bin path
     -genewise_bin      *genewise bin path

=head1 DESCRIPTION

B<identification of pseudogenes>

The main purpose of this program is to identify psedogenes using genblast and genewise tools based on protein sequences

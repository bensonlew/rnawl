# -*- coding: utf-8 -*-
# __author__ = 'haidong.gu'
# __last_modify__ = '2019/1/8'
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
from mbio.packages.metagenomic.common import check_command


class BlasrAgent(Agent):
    """
    version 1.0
    """

    def __init__(self, parent):
        super(BlasrAgent, self).__init__(parent)
        options = [
            # Input Files
            {"name": "query", "type": "infile", "format": "sequence.fasta"},
            # a multi-fasta file of reads.preferable to use plx.h5 or bax.h5 files because they
            # contain more rich quality value information
            {"name": "reference", "type": "string"},
            {"name": "safile", "type": "string"},
            # suffixArrayFile, use the suffix array 'sa' for detecting matches between the reads
            # and the reference. The suffix array has been prepared by the sawriter program
            {"name": "ctab", "type": "string"},
            # a table of tuple counts used to estimate match significance.This is by the program
            # 'printTupleCountTable'. While it is quick to generate on the fly, if there are many
            # invocations of blasr, it is useful to precompute the ctab
            {"name": "region", "type": "string"},
            # Read in a read-region table in HDF format for masking portions of reads.This may be
            # a single table if there is just one input file, or a fofn. When a region table is
            # specified, any region table inside the reads.plx.h5 or reads.bax.h5 files are ignored
            #-----Options for modifying reads-----#
            # There is ancilliary information about substrings of reads that is stored in a 'region
            # table' for each read file. Because HDF is used, the region table may be part of the
            # .bax.h5 or .plx.h5 file, or a seperate file. A contiguously read substring from the
            # template is a subread, and any read may contain multiple subreads. The boundaries of
            # the subreads may be inferred from the region table either directly or by definition
            # of adapter boundaries. Typically region tables also contain information for the
            # location of the high and low quality regions of reads. Reads produced by spurrious
            # reads from empty ZMWs have a high quality start coordinate equal to high quality end,
            # making no usable read.
            {"name": "noSplitSubreads", "type": "string", "default": "false"},
            # Do not split subreads at adapters. This is typically only useful when the genome in
            # an unrolled version of a known template, and contains
            # template-adapter-reverse_template sequence
            {"name": "polymerase", "type": "string"},
            # Instead of reconstituting ZMW reads, this option reconstitutes polymerase reads,
            # omitting LQ regions. Polymerase reads are aligne, if at least one subread is present.
            {"name": "useccs", "type": "string"},
            # Align the circular consensus sequence(ccs), then report alignments of the ccs subreads
            # to the window that the ccs was mapped to. Only alignment of the subreads are reported
            {"name": "useccsall", "type": "string"},
            # Similar to -useccs, except all subreads are aligned, rather than just the subreads
            # used to call the ccs. This will include reads that only cover part of the template
            {"name": "useccsdenovo", "type": "string"},
            # Align the circular consensus, and report only the alignment of the ccs sequence
            {"name": "ignoreRegions", "type": "string", "default": "false"},
            # Ignore any information in the region table
            {"name": "ignoreHQRegions", "type": "string", "default": "false"},
            # Ignore any hq regions in the region table
            #-----Alignments To Report ------#
            {"name": "num_alignment", "type": "int", "default": 10},  # bestn
            # Report the top 'n' alignments
            {"name": "hitPolicy", "type": "string"},
            # (all) Specify a policy to treat multiple hits from [all, allbest, random, randombest, leftmost]
            # all       : report all alignments.
            # allbest   : report all equally top scoring alignments.
            # random    : report a random alignment.
            # randombest: report a random alignment from multiple equally top scoring alignments.
            # leftmost  : report an alignment which has the best alignmentscore and has the smallest mapping coordinate in any reference.
            {"name": "placeRepeatsRandomly", "type": "string"},
            # DEPRECATED! If true, equivalent to --hitPolicy randombest.
            {"name": "placeGapConsistently", "type": "string"},
            # Place gaps consistently in alignments of a read as alignments of its reverse complementary sequence.
            {"name": "randomSeed", "type": "int"},
            # Seed for random number generator. By default (0), use current time as seed.
             {"name": "noSortRefinedAlignments", "type": "string", "default": "false"},
            # Once candidate alignments are generated and scored via sparse dynamic programming, they
            # are rescored using local alignment that accounts for differrent error profiles. Resorting
            # based on the local aligment may change the order the hits are returned
            {"name": "allowAdjacentIndels", "type": "string"},
            # When specified, adjacent insertion or deletions are allowed. Otherwise, adjacent insertion
            # and deletions are merged into one operation. Using quality values to guide pairwise
            # alignments may dictate that the higher probability alignment contains adjacent insertions
            # or deletions. Current tools such as GATK do not permit this and so they are not reported
            # by default
            # output formats and files
            {"name": "clipping", "type": "string", "default": "none"},
            # Use none/hard/subread/soft clipping for SAM output
            # {"name": "out", "type": "string"},
            # Write out to 'out'
            # {"name": "unaligned", "type": "string"},
            # output reads that are not aligned to 'file'
            {"name": "outfmt", "type": "int", "default": 5},
            # if not printing SAM, modify the output of the alignment.
            # 0 Print blast like output with |'s connecting matched nucleotides
            # 1 Print only a summary: score and pos
            # 2 Print in Compare.xml format
            # 3 Print in vulgar format (deprecated)
            # 4 Print a longer tabular version of the alignment
            # 5 Print in a machine-parsable format that is read by compareSequences.py
            {"name": "header", "type": "string"},
            # Print a header as the first line of the output file describing the contents of each column
            {"name": "titleTable", "type": "string"},
            # Construct a table of reference sequence titles. The reference sequences are enumerated by
            # row, 0, 1,... The reference index is printed in alignment results rather than the full
            # reference name. This makes output concise, particularly when very verbose titles exists in
            # reference names
            {"name": "minPctIdentity", "type": "float", "default": 0},
            # Only report alignments if they are greater than p percent identity
            {"name": "holeNumbers", "type": "string"},
            # When specified, only align reads whose ZMW hole numbers are in LIST.
            # LIST is a comma-delimited string of ranges, such as '1,2,3,10-13'.
            # This option only works when reads are in base or pulse h5 format.

            {"name": "printSAMQV", "type": "string", "default": "false"},
            # Print quality values to sam files
            #-----Options for anchoring alignment regions. This will have the greatest effect on speed and sensitivity.
            {"name": "minMatch", "type": "int", "default": "25"},
            # Minimum seed length. Higher minMatch will speed up alignment, but decrease sensitivity
            {"name": "maxMatch", "type": "int"},
            # Stop mapping a read to the genome when the lcp length reaches l. This is useful when the query
            # is part of the reference, for example when constructing pairwise alignments for de novo assembly
            {"name": "maxLCPLength", "type": "string"},
            # The same as maxMatch
            {"name": "maxAnchorsPerPosition", "type": "int"},
            # Do not add anchors from a position if it matches to more than 'm' locations in the target
            {"name": "advanceExactMatches", "type": "string"},
            # Another trick for speeding up alignments with match - E fewer anchors. Rather than finding
            # anchors between the read and the genome at every position in the read, when an anchor is
            # found at position i in a read of length L, the next position in a read to find an anchor
            # is at i+L-E. Use this when alignining already assembled contigs
            {"name": "ncandidates", "type": "int", "default": 70},
            # Keep up to 'n' candidates for the best alignment. A large value of n will slow mapping
            # because the slower dynamic programming steps are applied to more clusters of anchors
            # which can be a rate limiting step when reads are very long
            {"name": "concordant", "type": "string"},
            # Map all subreads of a zmw (hole) to where the longest full pass subread of the zmw aligned
            # to. This requires to use the region table and hq regions. This option only works when reads
            # are in base or pulse h5 format
            {"name": "fastMaxInterval", "type": "string"},
            # Fast search maxium increasing intervals as alignment candidates. The search is not as exhaustive
            # as the default, but is much faster
            {"name": "aggressiveIntervalCut", "type": "string"},
            # Agreesively filter out non-promising alignment candidates, if there exists at least one promising
            # candidate. If this option is turned on, Blasr is likely to ignore short alignments of ALU elements
            {"name": "fastSDP", "type": "string"},
            # Use a fast heuristic algorithm to speed up sparse dynamic programming
            #-----Options for Refining Hits.------
            {"name": "sdpTupleSize", "type": "string"},
            # Use matches of length K to speed dynamic programming alignments. This controls accuracy of assigning
            # gaps in pairwise alignments once a mapping has been found, rather than mapping sensitivity itself
            {"name": "scoreMatrix", "type": "string"},
            # Specify an alternative score matrix for scoring fasta reads. The matrix is in the format
            #   ACGTN
            # A abcde
            # C fghij
            # G klmno
            # T pqrst
            # N uvwxy    The value a....y should be input as a quoted space separated string: "a b c ... y".
            # lower scores are better, so matches should be less than mismatches e.g. a,g,m,s = -5(match), mismatch = 6
            {"name": "affineOpen", "type": "int"},
            # Set the penalty for opening an affine alignment
            {"name": "affineExtend", "type": "int"},
            # Change affine (extension) gap penalty. Lower value allows more gaps
            #-----Options for overlap/dynamic programming alignments and pairwise overlap for de novo assembly-----
            {"name": "useQuality", "type": "string"},
            # use substitution/insertion/deletion/merge quality values to score gap and mismatch penalties in pairwise
            # alignments. Because the insertion and deletion rates are much higher than subsitution, this will make
            # many alignments favor an insertion/deletion over a substitution. Naive consensus calling methods will
            # then often miss substitution polymorphisms. This option should be used when calling consensus using the
            # Quiver method. Furthermore, when not using quality values to score alignments, there will be a lower
            # consensus accuracy in homolymer regions
            {"name": "affineAlign", "type": "string"},
            # Refine alignment using affine guided align
            #-----Options for filtering reads-----
            {"name": "minReadLength", "type": "string"},
            # Skip reads that have a full length less than l. Subreads may be shorter
            {"name": "minSubreadLength", "type": "string"},
            # Do not align subreads of length less than l.
            {"name": "maxScore", "type": "string"},
            # Maximum score to output (high is bad, negative good)
            #----- Options for parallel alignment-----
            {"name": "num_threads", "type": "int", "default": 10},  # nproc
            # Align using N processes. All large data structures such as the suffix array and tuple count table are shared
            {"name": "start", "type": "string"},
            # Index of the first read to begin aligning. This is useful when multiple instances are running on the same
            # data, for example when on a multi-rack cluster
            {"name": "stride", "type": "string"},
            # Align one read every 'S' reads
            #-----Options for subsampling reads-----
            {"name": "subsample", "type": "string"},
            # Proportion of reads to randomly subsample (expressed as a decimal) and align
            {"name": "memory", "type": "int", "default": 20},
            # {"name": "outxml", "type": "outfile", "format": "align.blast.blast_xml"},
            {"name": "outtable", "type": "outfile", "format": "sequence.profile_table"}
        ]
        self.add_option(options)
        self.queue = 'BLAST'

    def check_options(self):
        """
        检查参数是否正确
        """
        if not self.option("query").is_set:
            raise OptionError("必须设置参数query")
        if not self.option("reference"):
            raise OptionError("必须设置参数reference")
        if not 0 <= self.option("outfmt") < 6:
            raise OptionError("outfmt.必须在0-6之间")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = self.option('num_threads')
        self._memory = str(self.option('memory')) + "G"


class BlasrTool(Tool):
    def __init__(self, config):
        super(BlasrTool, self).__init__(config)
        self.db = os.path.join(self.config.SOFTWARE_DIR, "database/taxon_db")
        if self.option("reference") == "silva128":
            self.reference = os.path.join(self.db, self.option("reference"), "silva.16s.fasta")
        elif self.option("reference") == "silva132":
            self.reference = os.path.join(self.db, self.option("reference"), "silva.16s.fasta")
        elif self.option("reference") == "greengene":
            self.reference = os.path.join(self.db, "Greengenes135/greengenes.16s.fasta")
        elif self.option("reference") == "nt":
            self.reference = os.path.join(self.db, self.option("reference"), "nt.fasta")
        self.cmd_path = "bioinfo/Genomic/Sofware/smrttools/smrtcmds/bin/blasr"

    def run_blasr(self):
        """
        description
        :return:
        """
        option_map = {
            "num_alignment": " --bestn",
            "outfmt": " -m",
            "minMatch": " --minMatch",
            "ncandidates": " --nCandidates",
            "num_threads": " --nproc"
        }
        query_name = os.path.splitext(os.path.basename(self.option("query").prop['path']))[0]
        db_name = os.path.splitext(os.path.basename(self.option("reference")))[0]
        outputfile = os.path.join(self.output_dir, query_name + "_vs_" + db_name)
        cmd = "%s %s %s --out %s" % (self.cmd_path, self.option("query").path, self.reference, outputfile)
        for option in option_map.keys():
            if self.option(option):
                cmd += " %s %s" % (option_map[option], self.option(option))
        command = self.add_command("blasr", cmd, ignore_error=True).run()
        self.wait(command)
        def success():
            self.logger.info("blasr运行完成")
        def fail():
            self.set_error("blasr运行出错")
        check_command(self, command, [0], [-9], success, fail, memory_limit_fun=None)
        self.option("outtable", outputfile)

    def set_output(self):
        """
        设置输出文件路径
        :return:
        """
        pass

    def run(self):
        super(BlasrTool, self).run()
        self.run_blasr()
        self.set_output()
        self.end()
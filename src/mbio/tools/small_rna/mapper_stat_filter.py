# coding=utf-8
import os
import glob
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import unittest
import ConfigParser
from mbio.packages.small_rna.arf_to_wig import ArfToWig

# import pandas as pd
__author__ = 'liubinxu'


class MapperStatFilterAgent(Agent):
    """
    fasta files remove duplication
    """
    def __init__(self, parent):
        super(MapperStatFilterAgent, self).__init__(parent)
        options = [
            {'type': 'infile', 'name': 'config', 'format': 'small_rna.ini'},
            {'format': 'small_rna.fasta', 'type': 'infile', 'name': 'fasta'},
            {'format': 'small_rna.fasta', 'type': 'infile', 'name': 'ref'},
            {'default': '', 'type': 'string', 'name': 'gtf'},
            {'default': None, 'type': 'string', 'name': 'samples'},
            {'format': 'small_rna.common', 'type': 'infile', 'name': 'arf'},
            {'default': 30, "type": "int", 'name': "chr_num"},
            {'format': 'small_rna.common_dir', 'type': 'infile', 'name': 'map_stat'},
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option('config').is_set:
            raise OptionError('必须设置config参数:基因组bowtie索引前缀')
        else:
            pass

    def set_resource(self):
        self._cpu = "1"
        self._memory = "{}G".format('10')

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", ""]
            ])
        """
        # more detail
        result_dir.add_regexp_rules([
            [r"*.xls", "xls", "xxx"],
            [r"*.list", "", "xxx"],
            ])
        """
        super(MapperStatFilterAgent, self).end()


class MapperStatFilterTool(Tool):
    """
    和基因组比对统计结果
    """
    def __init__(self, config):
        super(MapperStatFilterTool, self).__init__(config)
        software_dir = self.config.SOFTWARE_DIR
        self.perl_path = 'program/perl-5.24.0/bin/perl'
        self.perl = software_dir + '/program/perl-5.24.0/bin'
        self.mapper = software_dir + '/bioinfo/miRNA/mirdeep2/mapper.pl'
        self.mapper_stat = self.config.PACKAGE_DIR + "/small_rna/sRNA_genome_distribution_merge.pl"
        # self.mapper_filter = software_dir + '/bioinfo/miRNA/mirdeep2/'
        self.bowtie = self.config.SOFTWARE_DIR + '/bioinfo/align/bowtie-1.1.2'
        self.mirdeep_dir = software_dir + '/bioinfo/miRNA/mirdeep2'
        self.set_environ(PATH=self.bowtie)
        self.set_environ(PATH=self.mirdeep_dir)
        self.set_environ(PATH=self.perl)
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/bioinfo/WGS/circos-0.69-6/bin')
        self.circo_script_path = self.config.PACKAGE_DIR + '/small_rna/draw_smallrna_mapping.circos.pl'
        self.perl2_path = 'miniconda2/bin/perl'
        self.sample_names = list()
        self.chr_read = dict

    def get_top_chr(self, sample_stat, n):
        with open(sample_stat, 'r') as map_f:
            chr_read = dict()
            for line in map_f.readlines()[3:]:
                cols = line.split("\t")
                if cols[0].startswith("chromosome_"):
                    cols[0] = cols[0].replace("chromosome_", "")
                chr_read[cols[0]] = int(cols[5])
            chr_sorted = sorted(chr_read.items(), key=lambda x:x[1], reverse=True)

        return [x[0] for x in chr_sorted[:n]]


    def run_map_stat(self):
        '''
        比对结果统计
        '''
        fasta = os.path.basename(self.option("fasta").prop['path'])

        cmd = '{} {} '.format(self.perl_path, self.mapper_stat)
        cmd += '-{} {} '.format("config", self.option("config").prop['path'])
        cmd += '-{} {} '.format("fa", fasta)
        cmd += '-{} {} '.format("arf", self.option("arf").prop['path'])
        cmd_name = 'map_stat'
        command = self.add_command(cmd_name, cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("{} Finished successfully".format(cmd_name))
        elif command.return_code is None:
            self.logger.warn("{} Failed and returned None, we will try it again.".format(cmd_name))
            command.rerun()
            self.wait()
            if command.return_code is 0:
                self.logger.info("{} Finished successfully".format(cmd_name))
            else:
                self.set_error("{} Failed. >>>{}".format(cmd_name, cmd))
        else:
            self.set_error("{} Failed. >>>{}".format(cmd_name, cmd))

    def get_chr_list(self):
        chrom = list()
        with open(self.option('ref').prop['path'] + '.fai', 'r') as f:
            for line in f.readlines():
                chrom.append(line.split()[0])
        return chrom

    def filter_chr(self, chr_list):
        with open(self.option('ref').prop['path'] + '.fai', 'r') as chr_f, open("chr.txt", 'w') as chr_w:
            for line in chr_f:
                if line.strip().split("\t")[0] in chr_list:
                    chr_w.write(line)
        return "chr.txt"

    def arf_to_circos(self, config_file):
        '''
        转换arf文件为circos的输入
        '''

        config = ConfigParser.ConfigParser()
        config.read(config_file)
        sample2name = dict(config.items('NAME'))

        for sample,name in sample2name.items():
            if self.option("samples") and  (not name in self.option("samples")):
                continue
            chr_list = self.get_top_chr(os.path.join(self.option("map_stat").prop['path'], name + "_map_stat.xls"), self.option("chr_num"))
            chr_file = self.filter_chr(chr_list)
            A = ArfToWig()
            A.step_width = 500000
            A.set_fa(self.option('ref').prop['path'])
            A.set_arf(self.option('arf').prop['path'])
            A.init_mapping()
            A.parser_arf()
            # A.chr_list=['I', 'II', 'III', 'IV', 'V', 'X', 'MtDNA']
            # A.chr_list = self.get_chr_list()
            A.chr_list = chr_list
            A.sample = sample.upper()
            A.write_window(sample)
            self.sample_names.append(name)
            cmd = '{} {} --windows 500000 --pos {}.pos.window --neg {}.neg.window --chrlist {} --gff {} --outdir {}'.format(
                self.perl2_path, self.circo_script_path, sample, sample, chr_file, self.option('gtf'), name )
            cmd_name = 'map_circo' + name.lower()
            command = self.add_command(cmd_name, cmd)
            command.run()
            self.wait()
            if command.return_code == 0:
                self.logger.info("{} Finished successfully".format(cmd_name))
            elif command.return_code is None:
                self.logger.warn("{} Failed and returned None, we will try it again.".format(cmd_name))
                command.rerun()
                self.wait()
                if command.return_code is 0:
                    self.logger.info("{} Finished successfully".format(cmd_name))
                else:
                    self.set_error("{} Failed. >>>{}".format(cmd_name, cmd))
            else:
                self.set_error("{} Failed. >>>{}".format(cmd_name, cmd))


    def set_output(self):
        if os.path.exists(os.path.join(self.output_dir, 'qc_file.config')):
            os.remove(os.path.join(self.output_dir, 'qc_file.config'))
        os.link(self.option('config').prop['path'], os.path.join(self.output_dir, 'qc_file.config'))
        for name in self.sample_names:
            fname = name + '/circos.svg'
            link = os.path.join(self.output_dir, name + '_circos.svg')
            if os.path.exists(link):
                os.remove(link)
            os.link(fname, link)

            fname = name + '/circos.png'
            link = os.path.join(self.output_dir, name + '_circos.png')
            if os.path.exists(link):
                os.remove(link)
            os.link(fname, link)

    def run(self):
        super(MapperStatFilterTool, self).run()
        self.arf_to_circos(self.option('config').prop['path'])
        self.set_output()
        self.end()

class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        import datetime
        test_dir='/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/test_small_RNA/data0'

        data = {
            "id": "MapperStatFilter" + datetime.datetime.now().strftime('%H-%M-%S'),
            "type": "tool",
            "name": "small_rna.mapper_stat_filter",
            "instant": False,
            "options": dict(
                config=test_dir + "/" + "qc_file.config",
                fasta=test_dir + "/FastaUniq_out/" + "uniq.fasta",
                chr_num=3,
                map_stat = test_dir + "/mapper_result/All_map_stat.xls",
                arf=test_dir + "/mapper_result/reads_vs_genome.arf",
                gtf=test_dir + "/" + "Caenorhabditis_elegans.WBcel235.36.gtf",
                ref=test_dir + "/" + "Caenorhabditis_elegans.WBcel235.dna_rm.toplevel.fa"
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()

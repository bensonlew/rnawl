# coding=utf-8
import os
import glob
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import unittest
import ConfigParser
from mbio.packages.small_rna.arf_to_wig import ArfToWig
from Bio import SeqIO
# import pandas as pd
__author__ = 'liubinxu'


class MapperAndStatAgent(Agent):
    """
    fasta files remove duplication
    """
    def __init__(self, parent):
        super(MapperAndStatAgent, self).__init__(parent)
        options = [
            {'type': 'infile', 'name': 'config', 'format': 'small_rna.ini'},
            {'format': 'small_rna.fasta', 'type': 'infile', 'name': 'fasta'},
            {'format': 'small_rna.fasta', 'type': 'infile', 'name': 'ref'},
            {'format': 'gene_structure.gtf', 'type': 'infile', 'name': 'gtf'},
            {'default': '', 'type': 'string', 'name': 'index'},
            {'default': 'reads_vs_genome.arf', 'type': 'string', 'name': 'arf'},
            {'default': 10, "type": "int", 'name': "chr_num"},
            {'default': 10, 'type': 'int', 'name': 'p'},
        ]
        self.add_option(options)
        self._memory_increase_step = 60

    def check_options(self):
        if self.option('index') == '':
            raise OptionError('必须设置index参数:基因组bowtie索引前缀')
        else:
            if os.path.exists(self.option('index') + ".rev.1.ebwt") or os.path.exists(self.option('index') + ".rev.1.ebwtl"):
                pass
            else:
                pass
                # raise OptionError('bowtie 索引不存在')
        pass

    def set_resource(self):
        self._cpu = 2
        file_size = float(os.path.getsize(self.option('ref').prop['path'])) / 1024 / 1024 / 1024
        self._memory = "{}G".format(max(int(file_size * 2 +3), 10))

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
        super(MapperAndStatAgent, self).end()


class MapperAndStatTool(Tool):
    """
    和基因组比对统计结果
    """
    def __init__(self, config):
        super(MapperAndStatTool, self).__init__(config)
        software_dir = self.config.SOFTWARE_DIR
        self.perl_path = 'program/perl-5.24.0/bin/perl'
        self.perl = software_dir + '/program/perl-5.24.0/bin'
        self.mapper = software_dir + '/bioinfo/miRNA/mirdeep2/mapper.pl'
        self.mapper_stat = self.config.PACKAGE_DIR + "/small_rna/sRNA_genome_distribution_merge.pl"
        # self.mapper_filter = software_dir + '/bioinfo/miRNA/mirdeep2/'
        # self.bowtie = self.config.SOFTWARE_DIR + '/bioinfo/align/bowtie-1.1.2'
        self.bowtie = self.config.SOFTWARE_DIR + '/bioinfo/align/bowtie-1.2.3-linux-x86_64'
        self.mirdeep_dir = software_dir + '/bioinfo/miRNA/mirdeep2'
        self.python_path = software_dir + "/miniconda2/bin/"
        self.set_environ(PATH=self.python_path)
        self.set_environ(PATH=self.bowtie)
        self.set_environ(PATH=self.mirdeep_dir)
        self.set_environ(PATH=self.perl)
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/bioinfo/WGS/circos-0.69-6/bin')
        self.circo_script_path = self.config.PACKAGE_DIR + '/small_rna/draw_smallrna_mapping.circos.pl'
        self.perl2_path = 'miniconda2/bin/perl'
        self.sample_names = list()
        self.bowtie_build =  'bioinfo/align/bowtie-1.1.2/bowtie-build'

    def get_top_chr(self, n):
        with open(self.work_dir + "/All_map_stat.xls", 'r') as map_f:
            chr_read = dict()
            for line in map_f.readlines()[3:]:
                cols = line.split("\t")
                if cols[0].startswith("chromosome_"):
                    cols[0] = cols[0].replace("chromosome_", "")
                chr_read[cols[0]] = int(cols[1])
            chr_sorted = sorted(chr_read.items(), key=lambda x:x[1], reverse=True)

        return [x[0] for x in chr_sorted[:n]]

    def run_map_fasta(self):
        self.logger.info("choose mapped fasta")
        map_list = set()
        with open(self.option("arf"), 'r') as f:
            for line in f:
                s = line.strip()
                s_id = s.split("\t")[0]
                map_list.add(s_id)
        with open(self.work_dir + "/uniq_mapped.fasta", 'w') as fo:
            for seq in SeqIO.parse(self.option("fasta").prop['path'], "fasta"):
                seq_id = seq.id
                if seq_id in map_list:
                    fo.write(">{}\n{}\n".format(seq_id, seq.seq))


    def run_map_stat(self):
        '''
        比对结果统计
        '''
        fasta = os.path.basename(self.option("fasta").prop['path'])

        cmd = '{} {} '.format(self.perl_path, self.mapper_stat)
        cmd += '-{} {} '.format("config", self.option("config").prop['path'])
        cmd += '-{} {} '.format("fa", fasta)
        cmd += '-{} {} '.format("arf", self.option("arf"))
        cmd += '-{} {} '.format("ref_fai", self.option('ref').prop['path'] + '.fai')
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

    def run_bowtie_index(self):
        '''
        bowtie index
        '''
        cmd = '{} {} {}'.format(self.bowtie_build, self.option('ref').prop['path'], self.option("index"))

        cmd_name = 'bowtie index'
        command = self.add_command(cmd_name, cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("{} Finished successfully".format(cmd_name))
        else:
            self.set_error("{} Failed. >>>{}".format(cmd_name, cmd))

    def run_mapper(self):
        '''
        bowtie比对
        '''
        fasta = os.path.basename(self.option("fasta").prop['path'])
        if os.path.exists(fasta):
            os.remove(fasta)
        os.link(self.option("fasta").prop['path'], fasta)
        cmd = '{} {} '.format(self.perl_path, self.mapper)
        cmd += '{} '.format(fasta)
        cmd += '-{} {} '.format("p", self.option("index"))
        cmd += '-{} {} '.format("t", self.option("arf"))
        cmd += '-{} {} '.format("o", self.option("p"))
        cmd += '-c -j -q -v -n'
        cmd_name = 'mapper'
        command = self.add_command(cmd_name, cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            if os.path.exists(self.option("arf")) and os.path.getsize(self.option("arf")) > 0:
                pass
            else:
                self.add_state("memory_limit", "memory is low!")
                self.set_error("{} Failed. >>>{}".format(cmd_name, cmd))
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
        chr_list = self.get_top_chr(self.option("chr_num"))
        chr_file = self.filter_chr(chr_list)
        for sample,name in sample2name.items():
            A = ArfToWig()
            A.step_width = 500000
            A.set_fa(self.option('ref').prop['path'])
            A.set_arf(self.option('arf'))
            A.init_mapping()
            A.parser_arf()
            # A.chr_list=['I', 'II', 'III', 'IV', 'V', 'X', 'MtDNA']
            A.chr_list = self.get_chr_list()
            A.sample = sample.upper()
            A.write_window(sample)
            self.sample_names.append(name)
            cmd = '{} {} --windows 500000 --pos {}.pos.window --neg {}.neg.window --chrlist {} --gff {} --outdir {}'.format(
                self.perl2_path, self.circo_script_path, sample, sample, chr_file, self.option('gtf').prop['path'], name )
            cmd_name = 'mapcirco' + name.lower()
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

        if os.path.exists(os.path.join(self.output_dir, 'uniq_mapped.fasta')):
            os.remove(os.path.join(self.output_dir, 'uniq_mapped.fasta'))
        os.link(self.work_dir + '/uniq_mapped.fasta', self.output_dir + '/uniq_mapped.fasta')

        all_files = [self.option('arf')] + glob.glob(r'*_map_stat.xls')
        for each in all_files:
            fname = os.path.basename(each)
            link = os.path.join(self.output_dir, fname)
            if os.path.exists(link):
                os.remove(link)
            os.link(each, link)
        if self.option("gtf").is_set:
            for name in self.sample_names:
                fname = name + '/circos.svg'
                link = os.path.join(self.output_dir, name + '_circos.svg')
                if os.path.exists(link):
                    os.remove(link)
                os.link(fname, link)

    def run(self):
        super(MapperAndStatTool, self).run()
        if os.path.exists(self.option("index") + '.1.ebwt'):
            pass
        elif os.path.exists(self.option("index") + '.1.ebwtl') and os.path.exists(self.option("index") + '.rev.1.ebwtl'):
            pass
        else:
            self.run_bowtie_index()
        self.run_mapper()
        self.run_map_fasta()
        self.run_map_stat()
        if self.option("gtf").is_set:
            self.arf_to_circos(self.option('config').prop['path'])
        else:
            pass
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
            "id": "MapperAndStat" + datetime.datetime.now().strftime('%H-%M-%S'),
            "type": "tool",
            "name": "small_rna.mapper_and_stat",
            "instant": False,
            "options": dict(
                config=test_dir + "/" + "qc_file.config",
                fasta=test_dir + "/FastaUniq_out/" + "uniq.fasta",
                index=test_dir + "/" + "Caenorhabditis_elegans.WBcel235.dna_rm_index",
                gtf=test_dir + "/" + "Caenorhabditis_elegans.WBcel235.36.gtf",
                ref=test_dir + "/" + "Caenorhabditis_elegans.WBcel235.dna_rm.toplevel.fa"
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()

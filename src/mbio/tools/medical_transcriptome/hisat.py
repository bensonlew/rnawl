# -*- coding: utf-8 -*-
# __author__ = 'zengjing,shicaiping'

import os, glob
import unittest

from biocluster.agent import Agent
from biocluster.config import Config
from biocluster.core.exceptions import OptionError
from biocluster.tool import Tool
from pymongo.errors import ServerSelectionTimeoutError
from pymongo.errors import NetworkTimeout
import time

class HisatAgent(Agent):
    """
    对客户传入的参考基因组建索引，reads比对参考基因组
    version = 'hisat2-2.0.0'
    last_modify: by shicaiping at 20180508
    """

    def __init__(self, parent):
        super(HisatAgent, self).__init__(parent)
        options = [
            {"name": "ref_genome", "type": "string"},  # 参考基因组参数
            {"name": "dna_index", "type": "string"},  # 参考基因组参数
            {"name": "genome_version", "type": "string", "default": "Custom"},  # 参考基因组版本
            {"name": "genome_annot_version", "type": "string", "default": "Custom"},  # 参考基因组注释版本
            {"name": "mapping_method", "type": "string"},  # 比对软件，tophat or hisat
            {"name": "seq_method", "type": "string"},  # 测序方式，PE or SE
            {"name": "single_end_reads", "type": "infile", "format": "sequence.fastq"},
            {"name": "left_reads", "type": "infile", "format": "sequence.fastq"},
            {"name": "right_reads", "type": "infile", "format": "sequence.fastq"},
            {"name": "bam_output", "type": "outfile", "format": "align.bwa.bam"},
            {"name": "assemble_method", "type": "string"},
            {"name": "sample", "type": "string"},
            {"name": "strand_specific", "type": "bool", "default": False},
            {'name': 'strand_direct', 'type': 'string', 'default': 'none'},
        ]
        self.add_option(options)
        self.step.add_steps('hisat')
        self.on('start', self.step_start)
        self.on('end', self.step_end)
        self._memory_increase_step = 40

    def step_start(self):
        self.step.hisat.start()
        self.step.update()

    def step_end(self):
        self.step.hisat.finish()
        self.step.update()

    def check_option(self):
        if not self.option("seq_method") in ["PE", "SE"]:
            raise OptionError("请选择是双端测序还是单端测序", code="33707306")
        else:
            if self.option("seq_method") == "PE":
                if not self.option("single_end_reads").is_set:
                    raise OptionError("请传入单端测序文件", code="33707307")
            else:
                if not self.option("left_reads").is_set:
                    raise OptionError("请传入左端测序文件", code="33707308")
                if not self.option("right_reads").is_set:
                    raise OptionError("请传入右端测序文件", code="33707309")
        if not self.option("assemble_method").lower() in ["cufflinks", "stringtie", "none"]:
            raise OptionError("请选择拼接软件", code="33707310")

    def set_resource(self):
        self._cpu = 1
        self._memory = '20G'


class HisatTool(Tool):
    def __init__(self, config):
        super(HisatTool, self).__init__(config)
        self.program = {
            'hisat2': 'bioinfo/align/hisat2/hisat2-2.1.0/hisat2',
            'samtools': 'miniconda2/bin/samtools'
        }
        self.hisat_path = 'bioinfo/align/hisat2/hisat2-2.1.0/'
        self.ht2_idx = str()
        self.process_rerun = 0

    def run(self):
        super(HisatTool, self).run()
        self.hisat_build()
        self.end()

    def hisat_build(self):
        try:
            self.ht2_idx = os.path.join(
                self.config.SOFTWARE_DIR, 'database/Genome_DB_finish/{}'.format(self.option("dna_index")))
        except (ServerSelectionTimeoutError, NetworkTimeout): # 捕获因为mongo服务器问题导致的异常后重运行此方法
            if self.process_rerun < 5:
                self.process_rerun += 1
                self.logger.info("检测到TimeoutError, 第{}次重运行方法".format(self.process_rerun))
                time.sleep(5)
                self.hisat_build()
            else:
                self.add_state('memory_limit', '检测到TimeoutError, 重运行tool')

        self.hisat_mapping()

    def hisat_mapping(self):
        sort_tmp_bam = glob.glob(os.path.join(self.work_dir, 'accepted_hits.bam.tmp.*.bam'))
        if len(sort_tmp_bam) >= 1:
            for file in sort_tmp_bam:
                os.remove(file)
        if self.option("seq_method") == "PE":
            if self.option("assemble_method") == "cufflinks":
                cmd = "{} -p 8 -q --dta-cufflinks -x {} -1 {} -2 {} -S accepted_hits.unsorted.sam".format(
                    self.program['hisat2'], self.ht2_idx, self.option("left_reads").prop["path"],
                    self.option("right_reads").prop["path"])
            elif self.option("assemble_method") == "stringtie":
                cmd = "{} -p 8 -q --dta -x {} -1 {} -2 {} -S accepted_hits.unsorted.sam".format(
                    self.program['hisat2'], self.ht2_idx, self.option("left_reads").prop["path"],
                    self.option("right_reads").prop["path"])
            else:
                cmd = "{} -p 8 -q -x {} -1 {} -2 {} -S accepted_hits.unsorted.sam".format(
                    self.program['hisat2'], self.ht2_idx, self.option("left_reads").prop["path"],
                    self.option("right_reads").prop["path"])
        else:
            if self.option("assemble_method") == "cufflinks":
                cmd = "{} -p 8 -q --dta-cufflinks -x {} -U {} -S accepted_hits.unsorted.sam".format(
                    self.program['hisat2'], self.ht2_idx, self.option("single_end_reads").prop["path"])
            elif self.option("assemble_method") == "stringtie":
                cmd = "{} -p 8  -q --dta -x {} -U {} -S accepted_hits.unsorted.sam".format(
                    self.program['hisat2'], self.ht2_idx, self.option("single_end_reads").prop["path"])
            else:
                cmd = "{} -p 8 -q -x {} -U {} -S accepted_hits.unsorted.sam".format(
                    self.program['hisat2'], self.ht2_idx, self.option("single_end_reads").prop["path"])
        if self.option("strand_specific"):
            if self.option('strand_direct') == 'firststrand':
                cmd += " --rna-strandness RF"
            else:
                cmd += " --rna-strandness FR"
        self.logger.info("start running hisat2 for alignment")
        command = self.add_command("hisat_mapping", cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.samtools_sort()
        elif command.return_code == 1:
            self.logger.info("return code: {}".format(command.return_code))
            self.add_state('memory_limit', 'memory is low!')
        elif command.return_code is None:
            command.rerun()
            if command.return_code == 0:
                self.logger.info("succeed in running hisat2")
                self.samtools_sort()
        else:
            self.set_error("fail to run hisat2")

    def samtools_sort(self):
        cmd = '{} sort -m 2G -o {} -@ 4 {}'.format(
            self.program['samtools'],
            os.path.join(self.work_dir, "accepted_hits.bam"),
            os.path.join(self.work_dir, "accepted_hits.unsorted.sam")
        )
        command = self.add_command("samtools_sort", cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.set_output()
        elif command.return_code == 1:
            self.logger.info("return code: {}".format(command.return_code))
            self.add_state('memory_limit', 'memory is low!')
        elif command.return_code is None:
            command.rerun()
            if command.return_code == 0:
                self.logger.info("succeed in running samtools")
                self.set_output()
        else:
            self.set_error("fail to run samtools")

    def set_output(self):
        bam_fp = os.path.join(self.output_dir, '{}.bam'.format(self.option("sample")))
        if os.path.exists(bam_fp):
            os.remove(bam_fp)
        os.link(os.path.join(self.work_dir, "accepted_hits.bam"), bam_fp)
        self.option('bam_output').set_path(bam_fp)


class TestFunction(unittest.TestCase):
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'hisat_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'ref_rna_v2.hisat',
            'instant': False,
            'options': dict(
                ref_genome='Arabidopsis_thaliana',
                genome_version='TAIR10',
                genome_annot_version='Ensembl_43',
                seq_method='PE',
                left_reads='/mnt/ilustre/users/sanger-dev/workspace/20200313/Refrna_tsg_37039/FastpRna/output/fastq'
                           '/TR2_1.clean.1.fastq',
                right_reads='/mnt/ilustre/users/sanger-dev/workspace/20200313/Refrna_tsg_37039/FastpRna/output/fastq'
                            '/TR2_1.clean.2.fastq',
                assemble_method='stringtie',
                sample='TR2_1'
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)

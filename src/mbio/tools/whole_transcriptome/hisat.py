# -*- coding: utf-8 -*-
# __author__ = 'zengjing'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import shutil
import subprocess
import json
from biocluster.config import Config
import unittest

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
        self._memory_increase_step = 35

    def step_start(self):
        self.step.hisat.start()
        self.step.update()

    def step_end(self):
        self.step.hisat.finish()
        self.step.update()

    def check_option(self):
        """
        检查参数
        """
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
        """
        设置所需资源
        """
        self._cpu = 10
        if self.option("seq_method") == "PE":
            self._memory = '{}G'.format((os.path.getsize(self.option("left_reads").prop["path"]) + os.path.getsize(self.option("right_reads").prop["path"])) / 1024 ** 3 * 4 + 40)

class HisatTool(Tool):
    def __init__(self, config):
        super(HisatTool, self).__init__(config)
        self.hisat_path = 'bioinfo/align/hisat2/hisat2-2.1.0/'
        # self.samtools_path = self.config.SOFTWARE_DIR + '/bioinfo/align/samtools-1.3.1/'
        self.samtools_path = self.config.SOFTWARE_DIR + '/bioinfo/ref_rna_v2/miniconda2/bin/'
        # self.sort_path = self.config.SOFTWARE_DIR + '/bioinfo/align/samtools-1.3.1/'
        self.sort_path = self.config.SOFTWARE_DIR + '/bioinfo/ref_rna_v2/miniconda2/bin/'

    def hisat_build(self):
        """
        参考基因组准备
        """
        if self.option("ref_genome") == "customer_mode":

            cmd = "{}hisat2-build -f {} ref_index".format(self.hisat_path, self.
                                                          option("ref_genome_custom").prop['path'])
            self.logger.info("开始运行hisat2-build，进行建索引")
            command = self.add_command("hisat_build", cmd)
            command.run()
            self.wait(command)
            if command.return_code == 0:
                return True
            else:
                command.rerun()
                self.wait(command)
                if command.return_code == 0:
                    return True
                else:
                    self.set_error("建立索引出错", code="33707305")
        else:
            db = Config().get_mongo_client(mtype="ref_rna_v2", dydb_forbid=True)[Config().get_mongo_dbname("ref_rna_v2", dydb_forbid=True)]
            col = db["sg_genome_db"]
            self.logger.debug('name = {}'.format(self.option('ref_genome')))
            self.logger.debug('assembly = {}'.format(self.option('genome_version')))
            self.logger.debug('annot_version = {}'.format(self.option('genome_annot_version')))
            genome_info = col.find_one({
                "name": self.option("ref_genome"),
                "assembly": self.option("genome_version"),
                "annot_version": self.option("genome_annot_version")
            })
            rel_index = genome_info["dna_index"]
            global index_ref
            index_ref = self.config.SOFTWARE_DIR + "/database/Genome_DB_finish/" + rel_index

    def hisat_mapping(self):
        if self.option("ref_genome") == "customer_mode":
            ref_path = os.path.join(self.work_dir, "ref_index")
        else:
            ref_path = index_ref
        if self.option("seq_method") == "PE":
            if self.option("assemble_method") == "cufflinks":
                cmd = "{}hisat2 -p 8 -q --dta-cufflinks -x {} -1 {} -2 {} -S \
                accepted_hits.unsorted.sam". \
                    format(self.hisat_path, ref_path, self.option("left_reads").prop["path"],
                           self.option("right_reads").prop["path"])
            elif self.option("assemble_method") == "stringtie":
                cmd = "{}hisat2 -p 8 -q --dta -x {} -1 {} -2 {} -S accepted_hits.unsorted.sam". \
                    format(self.hisat_path, ref_path, self.option("left_reads").prop["path"],
                           self.option("right_reads").prop["path"])
            else:
                cmd = "{}hisat2 -p 8 -q -x {} -1 {} -2 {} -S accepted_hits.unsorted.sam". \
                    format(self.hisat_path, ref_path, self.option("left_reads").prop["path"],
                           self.option("right_reads").prop["path"])
        else:
            if self.option("assemble_method") == "cufflinks":
                cmd = "{}hisat2 -p 8 -q --dta-cufflinks -x {} {} -S accepted_hits.unsorted.sam". \
                    format(self.hisat_path, ref_path, self.option("single_end_reads").prop["path"])
            elif self.option("assemble_method") == "stringtie":
                cmd = "{}hisat2 -p 8  -q --dta -x {} {} -S accepted_hits.unsorted.sam". \
                    format(self.hisat_path, ref_path, self.option("single_end_reads").prop["path"])
            else:
                cmd = "{}hisat2 -p 8 -q -x {} {} -S accepted_hits.unsorted.sam". \
                    format(self.hisat_path, ref_path, self.option("single_end_reads").prop["path"])
        if self.option("strand_specific"):
            if self.option('strand_direct') == 'firststrand':
                cmd += " --rna-strandness RF"
            else:
                cmd += " --rna-strandness FR"
        self.logger.info("开始运行hisat2，进行比对")
        command = self.add_command("hisat_mapping", cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.sam_bam()
        elif command.return_code == 1:
            self.logger.info("return code: %s" % command.return_code)
            self.add_state('memory_limit', 'memory is low!')
        elif command.return_code == None:
            command.rerun()
            if command.return_code == 0:
                self.logger.info("hisat运行完成")
                self.sam_bam()
        else:
            self.logger.error("hisat运行出错")
            self.set_error("hisat运行出错", code="33707306")
            self.set_error("运行hisat出错", code="33707307")

    def sam_bam(self):
        """
        运行samtools，将hisat比对的结果文件sam文件转成bam文件
        """
        sam_path = os.path.join(self.work_dir, "accepted_hits.unsorted.sam")
        cmd = "{}samtools view -bS {} > accepted_hits.unsorted.bam".format(self.samtools_path, sam_path)
        self.logger.info("开始运行samtools，将sam文件转为bam文件")
        self.logger.info(cmd)
        try:
            subprocess.check_output(cmd, shell=True)
        except subprocess.CalledProcessError:
            self.set_error("samtools运行出错!", code="33707308")
        bam_path = os.path.join(self.work_dir, "accepted_hits.unsorted.bam")
        sort_cmd = "{}samtools sort {} > accepted_hits.bam".format(self.sort_path, bam_path)
        self.logger.info("开始运行samtools sort，将转换的bam文件进行排序")
        subprocess.check_output(sort_cmd, shell=True)
        self.wait()
        output = os.path.join(self.work_dir, "accepted_hits.bam")
        pre = self.option("sample")
        if pre.find("_sickele") != -1:
            pre = pre[:-9]
        if os.path.exists(self.output_dir + "/" + pre + ".bam"):
            os.remove(self.output_dir + "/" + pre + ".bam")
        os.link(output, self.output_dir + "/" + pre + ".bam")

    def run(self):
        """
        运行
        """
        super(HisatTool, self).run()
        self.hisat_build()
        self.hisat_mapping()
        self.sam_bam()
        self.end()

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
                ref_genome='Homo_sapiens',
                genome_version='GRCh38.p10',
                genome_annot_version='Ensemble_release_89',
                seq_method='PE',
                left_reads='/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/lnc_rna/test_data/23Y_1115_S6_L002_1P.fastq',
                right_reads='/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/lnc_rna/test_data/23Y_1115_S6_L002_2P.fastq',
                assemble_method= 'stringtie',
                sample='Con1'
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)

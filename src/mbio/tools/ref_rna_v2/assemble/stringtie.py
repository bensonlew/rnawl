# -*- coding: utf-8 -*-
# __author__ = 'wangzhaoyue'

import os
import unittest

from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
from biocluster.tool import Tool


class StringtieAgent(Agent):
    """
    有参转录组stringtie拼接
    version v1.0.1
    author: wangzhaoyue
    last_modify: 2016.09.06
    """

    def __init__(self, parent):
        super(StringtieAgent, self).__init__(parent)
        options = [
            {"name": "sample_bam", "type": "infile", "format": "align.bwa.bam"},  # 所有样本比对之后的bam文件
            {"name": "ref_fa", "type": "infile", "format": "ref_rna_v2.fasta"},  # 参考基因文件
            {"name": "ref_gtf", "type": "infile", "format": "gene_structure.gtf"},  # 参考基因的注释文件
            {"name": "cpu", "type": "int", "default": 10},  # stringtie软件所分配的cpu数量
            {"name": "fr_stranded", "type": "string", "default": "fr-unstranded"},  # 是否链特异性
            {"name": "strand_direct", "type": "string", "default": "none"},  # 链特异性时选择正负链
            {"name": "sample_gtf", "type": "outfile", "format": "gene_structure.gtf"},  # 输出的gtf文件
            {"name": "min_coverage", "type": "int", "default": 3},  # minimum junction coverage
            {"name": "min_read", "type": "int", "default": 5},  # minimum reads per bp coverage
        ]
        self.add_option(options)
        self.step.add_steps("stringtie")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)
        self._memory_increase_step = 20

    def stepstart(self):
        self.step.stringtie.start()
        self.step.update()

    def stepfinish(self):
        self.step.stringtie.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option('sample_bam'):
            raise OptionError('必须输入样本文件为bam格式', code="33703901")
        if not self.option('ref_fa'):
            raise OptionError('必须输入参考序列ref.fa', code="33703902")
        if not self.option('ref_gtf'):
            raise OptionError('必须输入参考序列ref.gtf', code="33703903")
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 11
        self._memory = "15G"

    def end(self):
        # result_dir = self.add_upload_dir(self.output_dir)
        # result_dir.add_relpath_rules([
        #     [".", "", "结果输出目录"],
        # ])
        # result_dir.add_regexp_rules([
        #     ["_out.gtf", "gtf", "样本拼接之后的gtf文件"]
        # ])
        super(StringtieAgent, self).end()


class StringtieTool(Tool):
    def __init__(self, config):
        super(StringtieTool, self).__init__(config)
        self._version = "v1.0.1"
        # self.stringtie_path = 'bioinfo/rna/stringtie-1.3.3b/'
        # self.stringtie_path = 'bioinfo/rna/stringtie-1.3.4d/'
        # self.stringtie_path = 'miniconda2/bin/'
        self.stringtie_path = "miniconda2/bin/"
        self.gffread_path = "bioinfo/rna/cufflinks-2.2.1/"

    def run(self):
        """
        运行
        :return:
        """
        super(StringtieTool, self).run()
        self.run_stringtie()
        self.run_gtf_to_fa()
        self.set_output()
        self.end()

    def run_stringtie(self):
        """
        运行stringtie软件，进行拼接组装
        """
        sample_name = os.path.basename(self.option('sample_bam').prop['path']).split('.bam')[0]
        if self.option("fr_stranded") == "fr-unstranded":  # 判断连特异性
            cmd = self.stringtie_path + 'stringtie %s -p %s -j %s -c %s -G %s  -o ' % (
                self.option('sample_bam').prop['path'], self.option('cpu'), self.option('min_coverage'),
                self.option('min_read'), self.option('ref_gtf').prop['path']) + sample_name + "_out.gtf"
        else:
            if self.option("strand_direct") == "firststrand":
                cmd = self.stringtie_path + 'stringtie %s --rf -p %s -j %s -c %s -G %s  -o ' % (
                    self.option('sample_bam').prop['path'], self.option('cpu'), self.option('min_coverage'),
                    self.option('min_read'), self.option('ref_gtf').prop['path']) + sample_name + "_out.gtf"
            else:
                cmd = self.stringtie_path + 'stringtie %s --fr -p %s -j %s -c %s -G %s  -o ' % (
                    self.option('sample_bam').prop['path'], self.option('cpu'), self.option('min_coverage'),
                    self.option('min_read'), self.option('ref_gtf').prop['path']) + sample_name + "_out.gtf"
        self.logger.info('运行stringtie软件，进行组装拼接')
        command = self.add_command("stringtie_cmd", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("stringtie运行完成")
        elif command.return_code == 1:
            self.logger.info("return code: %s" % command.return_code)
            self.add_state('memory_limit', 'memory is low!')
        else:
            self.set_error("stringtie运行出错!", code="33703904")

    def run_gtf_to_fa(self):
        """
        运行gtf_to_fasta，转录本gtf文件转fa文件
        """
        sample_name = os.path.basename(self.option('sample_bam').prop['path']).split('.bam')[0]
        cmd = self.gffread_path + "gffread %s -g %s -w %s_out.fa" % (
            self.work_dir + "/" + sample_name + "_out.gtf", self.option('ref_fa').prop['path'],
            self.work_dir + "/" + sample_name)
        self.logger.info('运行gtf_to_fasta，形成fasta文件')
        command = self.add_command("gtf_to_fa_cmd", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("gtf_to_fasta运行完成")
        else:
            self.set_error("gtf_to_fasta运行出错!", code="33703905")

    def set_output(self):
        """
        将结果文件link到output文件夹下面
        :return:
        """
        for root, dirs, files in os.walk(self.output_dir):
            for names in files:
                os.remove(os.path.join(root, names))
        self.logger.info("设置结果目录")
        try:
            sample_name = os.path.basename(self.option('sample_bam').prop['path']).split('.bam')[0]
            os.link(self.work_dir + "/" + sample_name + "_out.fa", self.output_dir +
                    "/" + sample_name + "_out.fa")
            os.link(self.work_dir + "/" + sample_name + "_out.gtf", self.output_dir + "/" + sample_name + "_out.gtf")
            self.option('sample_gtf').set_path(self.work_dir + "/" + sample_name + "_out.gtf")
            self.logger.info("设置组装拼接分析结果目录成功")

        except Exception as e:
            self.logger.info("设置组装拼接分析结果目录失败{}".format(e))
            self.set_error("设置组装拼接分析结果目录失败%s", variables=(e), code="33703906")


class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'stringtie_{}_{}'.format(random.randint(1000, 10000), random.randint(1000, 10000)),
            'type': 'tool',
            'name': 'ref_rna_v2.assemble.stringtie',
            'instant': False,
            'options': {
                'sample_bam': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/data/test_bam/Con1.bam',
                'ref_gtf': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/Ensemble_release_89/gtf/Homo_sapiens.GRCh38.89.gtf',
                'ref_fa': '/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Homo_sapiens/Ensemble_release_89/dna/Homo_sapiens.GRCh38.dna_rm.toplevel.clean.fa',
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTests([TestFunction('test')])
    unittest.TextTestRunner(verbosity=2).run(suite)

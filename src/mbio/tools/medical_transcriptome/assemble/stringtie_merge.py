# -*- coding: utf-8 -*-
# __author__ = 'wangzhaoyue'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import shutil
import re
import regex
import unittest

class StringtieMergeAgent(Agent):
    """
    有参转录组stringtie合并
    version v1.0.1
    author: wangzhaoyue
    last_modify: 2016.09.13
    """
    def __init__(self, parent):
        super(StringtieMergeAgent, self).__init__(parent)
        options = [
            {"name": "assembly_GTF_list.txt", "type": "infile", "format": "assembly.merge_txt"},
            # 所有样本比对之后的bam文件路径列表
            {"name": "ref_fa", "type": "infile", "format": "ref_rna_v2.fasta"},  # 参考基因文件
            {"name": "ref_gtf", "type": "infile", "format": "gene_structure.gtf"},  # 参考基因的注释文件
            {"name": "cpu", "type": "int", "default": 10},  # stringtie软件所分配的cpu数
            {"name": "merged_gtf", "type": "outfile", "format": "gene_structure.gtf"},  # 输出的合并文件
            {"name": "min_iso", "type": "float", "default": 0.1}, # minimum isoform fraction
            {"name": "min_cov", "type": "int", "default": 5}, # minimum input transcript coverage
            {"name": "min_tpm", "type": "int", "default": 1}, # minimum input transcript TPM
        ]
        self.add_option(options)
        self.step.add_steps("stringtie_merge")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.stringtie_merge.start()
        self.step.update()

    def stepfinish(self):
        self.step.stringtie_merge.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option('assembly_GTF_list.txt'):
            raise OptionError('必须输入所有样本gtf路径文件为txt格式', code = "33704001")
        if not self.option('ref_fa'):
            raise OptionError('必须输入参考序列ref.fa', code = "33704002")
        if not self.option('ref_gtf'):
            raise OptionError('必须输入参考序列ref.gtf', code = "33704003")
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 2
        self._memory = "10G"

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["merged.gtf", "gtf", "样本合并之后的gtf文件"]
        ])
        super(StringtieMergeAgent, self).end()


class StringtieMergeTool(Tool):
    def __init__(self, config):
        super(StringtieMergeTool, self).__init__(config)
        self._version = "v1.0.1"
        # self.stringtie_merge_path = 'bioinfo/rna/stringtie-1.3.4d/'
        self.stringtie_merge_path = 'bioinfo/ref_rna_v2/miniconda2/bin/'
        self.gffread_path = "bioinfo/rna/cufflinks-2.2.1/"
        tmp = os.path.join(self.config.SOFTWARE_DIR, self.stringtie_merge_path)
        tmp_new = tmp + ":$PATH"
        self.logger.debug(tmp_new)
        self.set_environ(PATH=tmp_new)

    def run(self):
        """
        运行
        :return:
        """
        super(StringtieMergeTool, self).run()
        self.run_stringtie_merge()
        self.run_filter_gtf()
        self.run_gtf_to_fa()
        self.set_output()
        self.end()

    def run_stringtie_merge(self):
        """
        运行stringtie软件，进行拼接合并
        """
        cmd = self.stringtie_merge_path + 'stringtie --merge {} -T {} -c {} -f {} -p {} -G {} -s {} -o {}merged_raw.gtf ' .format(
            self.option('assembly_GTF_list.txt').prop['path'], self.option('min_tpm'), self.option('min_cov'), self.option('min_iso'),
            self.option('cpu'), self.option('ref_gtf').prop['path'], self.option('ref_fa').prop['path'], self.work_dir+"/")
        self.logger.info('运行stringtie软件，进行拼接合并')
        command = self.add_command("stringtie_merge_cmd", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("stringtie_merge运行完成")
        else:
            self.set_error("stringtie_merge运行出错!", code = "33704004")
    def run_filter_gtf(self):
        """
        过滤stringtie merge出的异常转录本，新基因存在已知转录本
        """
        known_tran2gene = {}
        known_gtf = self.option('ref_gtf').prop['path']
        merge_raw_gtf = self.work_dir + "/merged_raw.gtf"
        merge_gtf = self.work_dir + "/merged.gtf"

        with open(known_gtf, "r") as file:
            for line in file:
                line = line.strip()
                content_m = regex.match(
                    r'^([^#]\S*?)\t+((\S+)\t+){7}(.*;)*((transcript_id|gene_id)\s+?\"(\S+?)\");.*((transcript_id|gene_id)\s+?\"(\S+?)\");(.*;)*$',
                    line)
                if content_m:
                    if 'transcript_id' in content_m.captures(6):
                        tran_id = content_m.captures(7)[0]
                        gene_id = content_m.captures(10)[0]
                    else:
                        tran_id = content_m.captures(10)[0]
                        gene_id = content_m.captures(7)[0]
                    known_tran2gene[tran_id] = gene_id

        with open(merge_raw_gtf, "r") as merge_raw_file, open(merge_gtf, "w") as merge_file:
            for line in merge_raw_file:
                line = line.strip()
                content_m = regex.match(
                    r'^([^#]\S*?)\t+((\S+)\t+){7}(.*;)*((transcript_id|gene_id)\s+?\"(\S+?)\");.*((transcript_id|gene_id)\s+?\"(\S+?)\");(.*;)*$',
                    line)
                if content_m:
                    if 'transcript_id' in content_m.captures(6):
                        tran_id = content_m.captures(7)[0]
                        gene_id = content_m.captures(10)[0]
                    else:
                        tran_id = content_m.captures(10)[0]
                        gene_id = content_m.captures(7)[0]
                    if known_tran2gene.has_key(tran_id):
                        if known_tran2gene[tran_id] == gene_id:
                            merge_file.write(line + "\n")
                        else:
                            self.logger.info('merge gene_id异常: {}为已知转录本，不属于基因{}'.format(tran_id,gene_id))
                    else:
                        merge_file.write(line + "\n")
                else:
                    merge_file.write(line + "\n")

    def run_gtf_to_fa(self):
        """
        运行gtf_to_fasta，转录本gtf文件转fa文件
        """
        cmd = self.gffread_path + "gffread %s -g %s -w merged.fa" % (
        self.work_dir + "/"+"merged.gtf", self.option('ref_fa').prop['path'])
        self.logger.info('运行gtf_to_fasta，形成fasta文件')
        command = self.add_command("gtf_to_fa_cmd", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("gtf_to_fasta运行完成")
        else:
            self.set_error("gtf_to_fasta运行出错!", code = "33704005")

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        if os.path.exists(self.output_dir + "/merged.gtf"):
            os.remove(self.output_dir + "/merged.gtf")
        if os.path.exists(self.output_dir + "/merged.fa"):
            os.remove(self.output_dir + "/merged.fa")
        try:
            os.link(self.work_dir + "/merged.gtf", self.output_dir + "/merged.gtf")
            self.option('merged_gtf').set_path(self.work_dir + "/merged.gtf")
            os.link(self.work_dir + "/merged.fa", self.output_dir + "/merged.fa")
            self.logger.info("设置拼接合并分析结果目录成功")

        except Exception as e:
            self.logger.info("设置拼接合并分析结果目录失败{}".format(e))
            self.set_error("设置拼接合并分析结果目录失败%s", variables = (e), code = "33704006")

class TestFunction(unittest.TestCase):
    '''
    This is test for the tool. Just run this script to do test.
    '''
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'stringtie_merge_{}_{}'.format(random.randint(1000, 9999), random.randint(1000, 9999)),
            'type': 'tool',
            'name': 'ref_rna_v2.assemble.stringtie_merge',
            'instant': False,
            'options': {
                'assembly_GTF_list.txt': '/mnt/ilustre/users/sanger-dev/sg-users/qinjincheng/lnc_rna/assemble/assembly_GTF_list.txt',
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

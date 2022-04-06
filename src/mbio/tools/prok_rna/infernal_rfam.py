# -*- coding: utf-8 -*-
# __author__ = 'zhangyitong'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re
import ConfigParser
import unittest
import glob


class InfernalRfamAgent(Agent):
    """
    已知miRNA鉴定
    """
    def __init__(self, parent):
        super(InfernalRfamAgent, self).__init__(parent)
        options = [
            {"name": "query", "type": "infile", "format": "prok_rna.fasta"},  # 输入序列文件
            {'name': 'sample_name', 'type': 'string', 'default': 'genome'},
            {"name": "evalue", "type": "float", "default": 1e-5},
            {"name": "rfam_result", "type": "outfile", "format": 'prok_rna.common'},
        ]
        self.add_option(options)
        self.step.add_steps("infernal")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.infernal.start()
        self.step.update()

    def stepfinish(self):
        self.step.infernal.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not (self.option("query").is_set or self.option("query_fq")):
            raise OptionError("必须提供输入文件", code = "35002401")
        self.set_resource()
        return True

    def set_resource(self):
        """
        设置所需资源，需在类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 8
        self._memory = '10G'

    def end(self):
        super(InfernalRfamAgent, self).end()


class InfernalRfamTool(Tool):
    def __init__(self, config):
        super(InfernalRfamTool, self).__init__(config)
        self.set_environ(PATH=os.path.join(self.config.SOFTWARE_DIR, 'bioinfo/prok_rna/miniconda3/bin/'))
        self.input_file = self.option("query").prop["path"]
        self.program = {
            'infernal': 'bioinfo/prok_rna/miniconda3/bin/cmscan',
            'python': 'miniconda2/bin/python',
            'seqkit': 'bioinfo/seq/seqkit',
        }
        self.script = {
            'rfam_stat': os.path.join(self.config.PACKAGE_DIR, "prok_rna/rfam_stat.py"),
        }
        self.file = {
            # 'rfam_seed': os.path.join(self.config.SOFTWARE_DIR, "database/prok_rna/sRNA_Anno/Rfam.seed"),
            # 'rfam_cm': os.path.join(self.config.SOFTWARE_DIR, "database/prok_rna/sRNA_Anno/Rfam.cm-14.1/Rfam.cm"),
            'rfam_seed': os.path.join(self.config.SOFTWARE_DIR, "database/prok_rna/sRNA_Anno/Rfam-14.6/Rfam.seed"),
            'rfam_cm': os.path.join(self.config.SOFTWARE_DIR, "database/prok_rna/sRNA_Anno/Rfam-14.6/Rfam.cm"),
            'tblout': os.path.join(self.work_dir, self.option('sample_name') + '_Rfam.tblout'),
        }

    def run(self):
        """
        运行
        :return:
        """
        super(InfernalRfamTool, self).run()
        self.length_count()
        self.set_output()
        self.end()

    def length_count(self):
        self.logger.info('开始统计输入序列长度')
        cmd = "{} fx2tab ".format(os.path.join(self.config.SOFTWARE_DIR, self.program['seqkit']))
        cmd += "--length --name "
        cmd += "{} > {}_length.txt".format(self.input_file, os.path.join(self.work_dir, self.option('sample_name')))
        self.logger.info("使用seqkit统计序列长度")
        command = self.add_command("seq_length", cmd, shell=True)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("使用seqkit统计序列长度完成!")
        else:
            self.set_error("使用seqkit统计序列长度出错！")
        self.infernal_rfam()

    def infernal_rfam(self):
        self.logger.info('开始进行Rfam预测')
        out_path = os.path.join(self.work_dir, self.option('sample_name'))
        total = 0
        with open('{}_length.txt'.format(out_path), 'r') as l:
            for line in l:
                length = line.strip().split('\t')[-1]
                total += int(length)
        self.logger.info(total)
        z_size = float(total)*2/1000000
        cmd = "{} -Z {} ".format(self.program['infernal'], z_size)
        cmd += "-o {}_Rfam.txt --tblout {}_Rfam.tblout ".format(out_path, out_path)
        cmd += "--fmt 2 --noali --rfam --nohmmonly --cpu 8 "
        cmd += "-E {} {} {}".format(self.option('evalue'), self.file['rfam_cm'], self.input_file)
        self.logger.info("使用Infernal进行Rfam预测")
        command = self.add_command("infernal_rfam", cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("使用Infernal进行Rfam预测完成!")
        else:
            self.set_error("使用Infernal进行Rfam预测出错！")

    def set_output(self):
        all_files = os.listdir(self.work_dir)
        all_files = [self.work_dir + '/' + each for each in all_files]
        for each in all_files:
            if each.endswith('.tblout') or each.endswith('Rfam.txt'):
                fname = os.path.basename(each)
                link = os.path.join(self.output_dir, fname)
                if os.path.exists(link):
                    os.remove(link)
                os.link(each, link)
        result = glob.glob(os.path.join(self.output_dir, '*.tblout'))
        if result:
            self.option("rfam_result").set_path(result[0])



class TestFunction(unittest.TestCase):

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "InfernalRfam" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "prok_rna.infernal_rfam",
            "instant": False,
            "options": dict(
                query="/mnt/ilustre/users/sanger-dev/workspace/20210602/Single_Srna_3999/Srna2/Rockhopper__1/Rockhopper_Results/genome.predicted_RNA.fa",
                evalue='1e-5'
                #config="/mnt/ilustre/users/sanger-dev/sg-users/shicaiping/miRNA/quantifier_test/Uniq.cfg.ini"
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()

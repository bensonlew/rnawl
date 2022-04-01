# -*- coding: utf-8 -*-
# __author__ = 'wangzhaoyue'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import shutil


class CuffcompareAgent(Agent):
    """
    有参转录组cuffcompare合并
    version v1.0.1
    author: wangzhaoyue
    last_modify: 2016.09.12
    """
    def __init__(self, parent):
        super(CuffcompareAgent, self).__init__(parent)
        options = [
            {"name": "merged_gtf", "type": "infile", "format": "gene_structure.gtf"},  # 拼接合并之后的转录本文件
            {"name": "ref_fa", "type": "infile", "format": "sequence.fasta"},  # 参考基因文件
            {"name": "ref_gtf", "type": "infile", "format": "gene_structure.gtf"},  # 参考基因的注释文件
            {"name": "cuff_gtf", "type": "outfile", "format": "gene_structure.gtf"},  # compare后的combined.gtf文件
        ]
        self.add_option(options)
        self.step.add_steps("cuffcompare")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.cuffcompare.start()
        self.step.update()

    def stepfinish(self):
        self.step.cuffcompare.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option('merged_gtf'):
            raise OptionError('必须输入所有样本拼接合并后gtf文件')
        if not self.option('ref_fa'):
            raise OptionError('必须输入参考序列ref.fa')
        if not self.option('ref_gtf'):
            raise OptionError('必须输入参考序列ref.gtf')
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 1
        self._memory = "5G"

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            [".gtf", "gtf", ""]
        ])
        super(CuffcompareAgent, self).end()


class CuffcompareTool(Tool):
    def __init__(self, config):
        super(CuffcompareTool, self).__init__(config)
        self._version = "v1.0.1"
        self.cuffcompare_path = 'bioinfo/rna/cufflinks-2.2.1/'
        tmp = os.path.join(self.config.SOFTWARE_DIR, self.cuffcompare_path)
        tmp_new = tmp + ":$PATH"
        self.logger.debug(tmp_new)
        self.set_environ(PATH=tmp_new)

    def run(self):
        """
        运行
        :return:
        """
        super(CuffcompareTool, self).run()
        self.run_cuffcompare()
        self.set_output()
        self.end()

    def run_cuffcompare(self):
        """
        运行cufflinks软件，进行拼接合并
        """
        cmd = self.cuffcompare_path + 'cuffcompare -s {} -C -o {}cuffcmp -r {} {}'.format(
            self.option('ref_fa').prop['path'], self.work_dir+"/", self.option('ref_gtf').prop['path'],
            self.option('merged_gtf').prop['path'])
        self.logger.info('运行cufflinks软件，进行比较')
        command = self.add_command("cuffcompare_cmd", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("cuffcompare运行完成")
        else:
            self.set_error("cuffcompare运行出错!")

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        try:
            shutil.copy2(self.work_dir + "/cuffcmp.combined.gtf", self.output_dir + "/cuffcmp.combined.gtf")
            self.option('cuff_gtf').set_path(self.work_dir + "/cuffcmp.combined.gtf")
            self.logger.info("设置拼接比较结果目录成功")
        except Exception as e:
            self.logger.info("设置拼接比较分析结果目录失败{}".format(e))
            self.set_error("设置拼接比较分析结果目录失败{}".format(e))

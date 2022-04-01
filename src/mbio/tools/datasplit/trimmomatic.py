# -*- coding: utf-8 -*-
# __author__ = 'wangzhaoyue'

"""Trimmomatic 去低值修剪过滤 """
import os
from biocluster.tool import Tool
from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError


class TrimmomaticAgent(Agent):
    """
    Trimmomatic
    """
    def __init__(self, parent=None):
        super(TrimmomaticAgent, self).__init__(parent)
        options = [
            {'name': 'fq1', 'type': "infile", "format": "sequence.fastq"},
            {'name': 'fq2', 'type': "infile", "format": "sequence.fastq"},
            {'name': 'fq_type', 'type': "string", "default": "PE"},  # PE or SE
            {'name': 'illuminaclip', 'type': "string"},  # 2:30:10
            {'name': 'leading', 'type': "string", "default": "0"},  # 切除首端碱基质量小于0的碱基或者N
            {'name': 'tailing', 'type': "string", "default": "20"},  # 切除末端碱基质量小于20的碱基或者N
            {'name': 'sliding_window', 'type': "string", "default": "50:20"},  # 例50:20  Windows的size是50个碱基，其平均碱基质量小于20，则切除
            {'name': 'minlen', 'type': "string", "default": "50"},  # 最低reads长度
            {'name': 'out_fq1', 'type': "outfile", "format": "sequence.fastq"},
            {'name': 'out_fq2', 'type': "outfile", "format": "sequence.fastq"},  # 去低值之后的fq序列
            {'name': 'lib_name', 'type': "string"},  # 文库/样本名

        ]
        self.add_option(options)
        self.step.add_steps("trimmomatic")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.trimmomatic.start()
        self.step.update()

    def stepfinish(self):
        self.step.trimmomatic.finish()
        self.step.update()

    def check_options(self):
        """
        参数检测
        """
        if not self.option('fq1'):
            raise OptionError('必须输入1端序列', code="31802301")
        if not self.option('fq2'):
            raise OptionError('必须输入2端序列', code="31802302")
        if not self.option('lib_name'):
            raise OptionError('必须填写输出文件的名称，文库名或样本名', code="31802303")
        return True

    def set_resource(self):
        """
        设置所需要的资源
        """
        self._cpu = 2
        self._memory = '20G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        super(TrimmomaticAgent, self).end()


class TrimmomaticTool(Tool):
    """
    """
    def __init__(self, config):
        super(TrimmomaticTool, self).__init__(config)
        self.javaPath = "program/sun_jdk1.8.0/bin/java"
        self.trimmomatic_path = self.config.SOFTWARE_DIR + '/bioinfo/seq/trimmomatic-0.36/trimmomatic-0.36.jar '
        self.database = self.config.SOFTWARE_DIR + '/database/datasplit/adaptor.list'   # 微生物基因组使用
        # self.lib_name = ''  # 文库/样本名

    def run(self):
        super(TrimmomaticTool, self).run()
        self.run_trimmomatic()
        self.set_output()
        self.end()

    def run_trimmomatic(self):
        """
        运行trimmomatic-0.36.jar
        version 0.36
        """
        # self.lib_name = os.path.basename(self.option('fq1').prop['path']).split('.all.raw.valid.1.fq')[0]  # 获取文库/样本名
        cmd = self.javaPath + " -jar " + self.trimmomatic_path + ' %s -phred33 %s %s %s %s %s %s LEADING:%s TRAILING:%s SLIDINGWINDOW:%s MINLEN:%s' % (
            self.option("fq_type"), self.option("fq1").prop["path"], self.option("fq2").prop["path"],
            self.work_dir + "/" + self.option('lib_name') + ".trim.1.fq", self.work_dir + "/" + self.option('lib_name') + ".unpair.1.fq",
            self.work_dir + "/" + self.option('lib_name') + ".trim.2.fq", self.work_dir + "/" + self.option('lib_name') + ".unpair.2.fq",
            self.option("leading"), self.option("tailing"), self.option("sliding_window"), self.option("minlen"))
        if self.option('illuminaclip'):
            cmd += " ILLUMINACLIP:%s:%s" % (self.database, self.option('illuminaclip'))
        self.logger.info('运行trimmomatic，进行过滤')
        command = self.add_command("trimmomatic_cmd", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("trimmomatic运行完成")
        else:
            self.set_error("trimmomatic运行出错!", code="31802301")

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        for f in os.listdir(self.output_dir):
            os.remove(os.path.join(self.output_dir, f))
        try:
            os.link(self.work_dir + "/" + self.option('lib_name') + ".trim.1.fq",
                    self.output_dir + '/' + self.option('lib_name') + ".trim.1.fq")
            os.link(self.work_dir + "/" + self.option('lib_name') + ".trim.2.fq",
                    self.output_dir + '/' + self.option('lib_name') + ".trim.2.fq")
            self.option('out_fq1').set_path(self.output_dir + '/' + self.option('lib_name') + ".trim.1.fq")
            self.option('out_fq2').set_path(self.output_dir + '/' + self.option('lib_name') + ".trim.2.fq")
            self.logger.info("设置trimmomatic分析结果目录成功")

        except Exception as e:
            self.logger.info("设置trimmomatic分析结果目录失败{}".format(e))
            self.set_error("设置trimmomatic分析结果目录失败%s", variables=(e), code="31802302")

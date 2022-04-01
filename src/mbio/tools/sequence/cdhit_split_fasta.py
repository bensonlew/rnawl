# -*- coding: utf-8 -*-
# __author__ = 'zouxuan'
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
import shutil


class CdhitSplitFastaAgent(Agent):
    """
    cd-hit-div
    version v1.0
    author: zouxuan
    last modified:2017.9.12
    """

    def __init__(self, parent):
        super(CdhitSplitFastaAgent, self).__init__(parent)
        options = [
            {"name": "gene_tmp_fa", "type": "infile", "format": "sequence.fasta"},  # 输入序列
            {"name": "number", "type": "int", "default": 1},  # 切分为几份
            {"name": "ou_dir", "type": "string", "default": ""},  # 输出路径
            {"name": "order", "type": "int", "default": 1},  # 0切分时不排序，1切分后从长到短排序
            {"name": "pre", "type": "string", "default": "gene.geneset.tmp.fa.div"}  # 文件前缀
        ]
        self.add_option(options)
        self.step.add_steps('cdhitsplitfasta')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.cdhitsplitfasta.start()
        self.step.update()

    def step_end(self):
        self.step.cdhitsplitfasta.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        """
        if not self.option("gene_tmp_fa").is_set:
            raise OptionError("必须设置参数gene_tmp_fa", code="34000201")
        if self.option("number") <= 0:
            raise OptionError("切割份数必须大于等于1", code="34000202")
        if not self.option("ou_dir").strip():
            raise OptionError("必须设置参数ou_dir", code="34000203")
        if not self.option("order") in [0, 1]:
            raise OptionError("参数order只能为0或1", code="34000204")

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 16
        self._memory = str(os.path.getsize(self.option("gene_tmp_fa").prop['path']) / 100000000 + 2) + 'G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        super(CdhitSplitFastaAgent, self).end()


class CdhitSplitFastaTool(Tool):
    def __init__(self, config):
        super(CdhitSplitFastaTool, self).__init__(config)
        self._version = '1.0'
        self.div_path = 'bioinfo/uniGene/cd-hit-v4.6.1-2012-08-27/cd-hit-div'  # 该软件切分前先排序
        self.perl_path = '/program/perl/perls/perl-5.24.0/bin/perl '
        self.div_path1 = self.config.SOFTWARE_DIR + '/bioinfo/uniGene/cd-hit-v4.6.1-2012-08-27/cd-hit-div_new.pl'  # 改脚本不按长短切分

    def run(self):
        super(CdhitSplitFastaTool, self).run()
        self.div()
        self.set_output()

    #    def div_num(self):
    #        filesize = os.path.getsize(self.option("gene_tmp_fa").prop['path'])
    #        n = filesize/500000000 + 1
    #        return n

    def div(self):
        n = self.option("number")
        if os.path.exists(self.option("ou_dir")):
            pass
        else:
            os.mkdir(self.option("ou_dir"))
        if self.option("order") == 1:
            cmd = '%s -i %s -o %s -div %s' % (
                self.div_path, self.option("gene_tmp_fa").prop['path'],
                os.path.join(self.option("ou_dir"), self.option("pre")), n)
        elif self.option("order") == 0:
            self.option("gene_tmp_fa").get_info()
            base = self.option("gene_tmp_fa").prop['bases']
            cmd = '%s %s %s %s %s %s' % (self.perl_path, self.div_path1, self.option("gene_tmp_fa").prop['path'],
                                         os.path.join(self.option("ou_dir"), self.option("pre")), n, base)
        self.logger.info(cmd)
        command1 = self.add_command('cmd_1', cmd)
        command1.run()
        self.wait(command1)
        if command1.return_code == 0:
            self.logger.info("div succeed")
        else:
            self.set_error("div failed", code="34000201")
            raise Exception("div failed")

    def set_output(self):
        self.linkdir(self.option("ou_dir"), "")
        self.end()

    def linkdir(self, dirpath, dirname):
        """
        link文件夹下的所有文件到本module的output目录
        """
        allfiles = os.listdir(dirpath)
        newdir = os.path.join(self.output_dir, dirname)
        if not os.path.exists(newdir):
            os.mkdir(newdir)
        oldfiles = [os.path.join(dirpath, i) for i in allfiles if os.path.isfile(os.path.join(dirpath, i))]
        newfiles = [os.path.join(newdir, i) for i in allfiles if os.path.isfile(os.path.join(dirpath, i))]
        for newfile in newfiles:
            if os.path.exists(newfile):
                os.remove(newfile)
        for i in range(len(newfiles)):
            os.link(oldfiles[i], newfiles[i])

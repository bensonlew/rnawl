# -*- coding: utf-8 -*-
# __author__ = 'zouxuan'
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
import subprocess
from biocluster.core.exceptions import OptionError
import shutil


class PanClusterAgent(Agent):
    """
    cd-hit-est
    version v1.0
    author: zouxuan
    last modified:2018.7.17
    """

    def __init__(self, parent):
        super(PanClusterAgent, self).__init__(parent)
        options = [
            {"name": "query", "type": "infile", "format": "sequence.fasta"},  # 输入fasta文件
            {"name": "qunum", "type": "int", "default": 0},  # fasta编号
            {"name": "identity", "type": "float", "default": 0.95},  ##给出cdhit的参数identity
            {"name": "coverage", "type": "float", "default": 0.9},  # 给出cdhit的参数coverage
            {"name": "memory_limit", "type": "int", "default": 10000},  # 内存大小，0为无限制
            {"name": "compare_method", "type": "int", "default": 0},  # 1为全局比对，0为局部比对
            {"name": "direction", "type": "int", "default": 1},  # 1为双向比对，0为单向比对
            {"name": "num_threads", "type": "int", "default": 2},  # cpu数
            {"name": "select", "type": "int", "default": 1},  # 1为聚类到最相似的类中，0为聚类到第一个符合阈值的类
            {"name": "compare", "type": "string", "default": ""},  # 比对结果输出路径
            {"name": "pre", "type": "string", "default": "gene.geneset.tmp.fa.div-"},  # 文件前缀
            {"name": "output", "type": "outfile", "format": "sequence.fasta"}, # 输出fasta文件
            {"name": "ana_type", "type": "string", "default": "nucl"},  # 输入分析类型，是对核酸聚类还是随蛋白聚类
            {"name": "method", "type": "string"}, #输入方法类型{'orthofinder', 'orthomcl', 'get_homologus', 'pgap', 'roary'}
        ]
        self.add_option(options)
        self.step.add_steps('cluster')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.cluster.start()
        self.step.update()

    def step_end(self):
        self.step.cluster.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        """
        if not self.option("query").is_set:
            raise OptionError("必须设置参数query")
        if not self.option("compare").strip():
            raise OptionError("必须设置输出路径compare")
        if not 0.75 <= self.option("identity") <= 1:
            raise OptionError("identity必须在0.75，1之间")
        if not 0 <= self.option("coverage") <= 1:
            raise OptionError("coverage必须在0,1之间")
        if not self.option("method"):
            raise OptionError("必须提供method方法")

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = self.option("num_threads")
        self._memory = str(self.option("memory_limit") / 1000) + 'G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        super(PanClusterAgent, self).end()


class PanClusterTool(Tool):
    def __init__(self, config):
        super(PanClusterTool, self).__init__(config)
        self._version = '1.0'
        self.cdhit_est_path = 'bioinfo/uniGene/cd-hit-v4.6.1-2012-08-27/cd-hit-est'
        self.cdhit_prot_path = 'bioinfo/uniGene/cd-hit-v4.6.1-2012-08-27/cd-hit'

    def run(self):
        """
        开始运行
        :return:
        """
        self.logger.info('开始运行')
        super(PanClusterTool, self).run()
        self.run_cdhit()
        self.set_output()

    def word_len(self):
        """
        设置word_length长度
        :return:
        """
        word_length = 8
        if self.option("identity") >= 0.9:
            word_length = 8
        elif 0.88 <= self.option("identity") < 0.9:
            word_length = 7
        elif 0.85 <= self.option("identity") < 0.88:
            word_length = 6
        elif 0.8 <= self.option("identity") < 0.85:
            word_length = 5
        elif 0.75 <= self.option("identity") < 0.8:
            word_length = 4
        return word_length

    def run_cdhit(self):
        length = self.word_len()
        out_dir = self.option("compare") + '/' + self.option("pre") + str(self.option("qunum")) + "-"
        if os.path.exists(out_dir):
            pass
        else:
            os.mkdir(out_dir)
        if self.option("ana_type") == "nucl":
            cmd = '%s -i %s -o %s -c %s -aS %s -n %s -G %s -M %s -d %s -r %s -g %s -T %s' % (
                self.cdhit_est_path, self.option("query").prop['path'], out_dir + "/o", self.option("identity"),
                self.option("coverage"), 11, 1, self.option("memory_limit"), 0,
                self.option("direction"), self.option("select"), self.option("num_threads"))
        else:
            cmd = '%s -i %s -o %s -c %s -aS %s -n %s -G %s -M %s -d %s -g %s -T %s' % (
                self.cdhit_prot_path, self.option("query").prop['path'], out_dir + "/o", self.option("identity"),
                self.option("coverage"), 5, self.option("compare_method"), self.option("memory_limit"), 0, self.option("select"), self.option("num_threads"))
        self.logger.info(cmd)
        command1 = self.add_command('cmd_1', cmd)
        command1.run()
        self.wait(command1)
        if command1.return_code == 0:
            # f=open(out_dir + "/o",'r')
            # line=f.readline()
            # if line.startswith('>'):
            #     self.logger.info("compare single succeed")
            # 修改读文件的方法为，fasta检查 by ghd @ 20180717
            try:
                self.option('output', out_dir + '/o')
                self.logger.info("compare single succeed at fist time")
            except:
            # else:
                command1.rerun()
                self.wait(command1)
                if command1.return_code == 0:
                    self.option('output', out_dir + "/o")
                    self.logger.info("compare single succeed at second time")
                else:
                    self.set_error("compare single failed", code="31600201")
                    raise Exception("compare single failed")
        else:
            self.set_error("compare single failed", code="31600202")
            raise Exception("compare single failed")

    def run_orthofinder(self):
        """
        运行orthofinder进行聚类
        :return:
        """


    def run_orthomcl(self):
        """
        运行orthomcl进行聚类
        :return:
        """


    def run_homologus(self):
        """
        运行Get_homologus软件进行聚类
        :return:
        """


    def run_pgap(self):
        """
        运行PGAP软件进行聚类
        :return:
        """

    def run_roary(self):
        """
        运行Roary软件进行聚类
        :return:
        """

    def set_output(self):
        # self.linkdir(self.option("compare") + '/gene.geneset.tmp.fa.div-' + str(self.option("qunum")) + "-",
        #              "gene.geneset.tmp.fa.div-" + str(self.option("qunum")) + "-")
        newdir = os.path.join(self.output_dir, self.option("pre") + str(self.option("qunum")) + "-")
        if not os.path.exists(newdir):
            os.mkdir(newdir)
        oldfiles = os.path.join(self.option("compare") + '/' + self.option("pre") + str(self.option("qunum")) + "-", 'o')
        newfiles = os.path.join(newdir, 'o')
        if os.path.exists(newfiles):
            os.remove(newfiles)
        os.link(oldfiles,newfiles)
        self.end()
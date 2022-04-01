# -*- coding: utf-8 -*-
# __author__ = 'yuguo'

"""RDP taxon 物种分类工具"""

from biocluster.tool import Tool
from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError
import os
import subprocess


class QiimeAssignAgent(Agent):
    """
    Qiime taxon_assign.py
    version v1.0
    """
    def __init__(self, parent=None):
        """
        """
        super(QiimeAssignAgent, self).__init__(parent)
        options = [
            {'name': 'fasta', 'type': 'infile', 'format': 'sequence.fasta'},  # 输入fasta文件
            {'name': 'revcomp', 'type': 'bool'},  # 序列是否翻转
            {'name': 'confidence', 'type': 'float', 'default': 0.7},  # 置信度值
            # {"name": "customer_mode", "type": "bool", "default": False},  # customer 自定义数据库
            {'name': 'database', 'type': 'string'},  # 数据库选择
            {'name': 'ref_fasta', 'type': 'infile', 'format': 'sequence.fasta'},  # 参考fasta序列
            {'name': 'ref_taxon', 'type': 'infile', 'format': 'taxon.seq_taxon'},  # 参考taxon文件
            {'name': 'taxon_file', 'type': 'outfile', 'format': 'taxon.seq_taxon'}  # 输出序列的分类信息文件
        ]
        self.add_option(options)
        self.step.add_steps('qiime_assign')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.qiime_assign.start()
        self.step.update()

    def step_end(self):
        self.step.qiime_assign.finish()
        self.step.update()

    def check_options(self):
        """
        检查参数设置
        """
        if not self.option("fasta").is_set:
            raise OptionError("必须设置参数fasta")
        if self.option("revcomp") not in [True, False]:
            raise OptionError("必须设置序列是否翻转")
        if self.option('database') == "custom_mode":
            if not self.option("ref_fasta").is_set or not self.option("ref_taxon").is_set:
                raise OptionError("数据库自定义模式必须设置参考fasta序列和参考taxon文件")
        else:
            if self.option("database") not in ['silva119/16s_bacteria', 'silva119/16s_archaea', 'silva119/16s', 'silva119/18s_eukaryota', 'unite7.0/its_fungi', 'fgr/amoA', 'fgr/nosZ', 'fgr/nirK', 'fgr/nirS', 'fgr/nifH', 'fgr/pmoA', 'fgr/mmoX']:
                raise OptionError("数据库{}不被支持".format(self.option("database")))

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
            ["seqs_tax_assignments.txt", "xls", "OTU的分类学信息"]
        ])
        super(QiimeAssignAgent, self).end()

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 10
        self._memory = '50000M'


class QiimeAssignTool(Tool):
    """
    Qiime Taxon Classify tool
    """
    def __init__(self, config):
        super(QiimeAssignTool, self).__init__(config)
        self.qiime_path = "program/Anaconda2/bin/"

    def run_prepare(self):
        if self.option('revcomp'):
            self.logger.info("revcomp 输入的fasta文件")
            try:
                subprocess.check_output(self.config.SOFTWARE_DIR+"/bioinfo/seq/fastx_toolkit_0.0.14/revcomp "+self.option('fasta').prop['path']+" > seqs.fasta", shell=True)
                self.logger.info("OK")
                return True
            except subprocess.CalledProcessError:
                self.logger.info("revcomp 出错")
                return False
        else:
            self.logger.info("链接输入文件到工作目录")
            if os.path.exists(self.work_dir+'/seqs.fasta'):
                os.remove(self.work_dir+'/seqs.fasta')
            os.link(self.option('fasta').prop['path'], self.work_dir+"/seqs.fasta")
            self.logger.info("OK")
            return True

    def run_assign(self):
        ref_fas = self.config.SOFTWARE_DIR+"/database/taxon_db/"+self.option('database')+'.fasta'
        ref_tax = self.config.SOFTWARE_DIR+"/database/taxon_db/"+self.option('database')+'.tax'
        if self.option('database') == "custom_mode":
            ref_fas = self.option('ref_fasta').prop['path']
            ref_tax = self.option('ref_taxon').prop['path']
        # export RDP_JAR_PATH=$HOME/app/rdp_classifier_2.2/rdp_classifier-2.2.jar"
        self.set_environ(RDP_JAR_PATH=self.config.SOFTWARE_DIR+"/bioinfo/taxon/rdp_classifier_2.2/rdp_classifier-2.2.jar")
        cmd = self.qiime_path+"assign_taxonomy.py  -m rdp -i seqs.fasta -c "+str(self.option('confidence'))+"  -r "+ref_fas+" -t "+ref_tax+" -o .  --rdp_max_memory 50000"
        self.logger.info(u"生成命令: "+cmd)
        assign = self.add_command("assign", cmd)
        self.logger.info("开始运行assign")
        assign.run()
        self.wait(assign)
        if assign.return_code == 0:
            self.logger.info("assign运行完成")
            subprocess.check_output("cat " + self.work_dir + "/seqs_tax_assignments.txt|sed  's/Unclassified/d__Unclassified/' > " + self.work_dir + "/seqs_tax_assignments.fix.txt", shell=True)
            os.system('rm -rf '+ self.output_dir)
            os.system('mkdir '+ self.output_dir)
            os.link(self.work_dir + '/seqs_tax_assignments.fix.txt', self.output_dir+'/seqs_tax_assignments.txt')
            self.option('taxon_file').set_path(self.output_dir+'/seqs_tax_assignments.txt')
        else:
            self.set_error("assign运行出错!")

    def run(self):
        super(QiimeAssignTool, self).run()
        if self.run_prepare():
            self.run_assign()
            self.end()
        else:
            self.set_error("run_prepare运行出错!")

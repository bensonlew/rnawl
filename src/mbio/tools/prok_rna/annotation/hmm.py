# -*- coding: utf-8 -*-
# __author__ = 'litangjian,zoujiaxun'

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
import shutil
import pandas as pd
from biocluster.core.exceptions import OptionError
from mbio.packages.prok_rna.pfam_domtblout import pfam_out
from mbio.packages.rna.annot_config import AnnotConfig

class HmmAgent(Agent):
    """
    transdecoder:Hmm,cds预测软件
    hmmscan：Pfam数据库比对工具
    version 1.0
    author: qindanhua
    last_modify: 2016.07.11

    now
    hmmscan  Version 3.2.1
    last_modify: 20200110
    """

    def __init__(self, parent):
        super(HmmAgent, self).__init__(parent)
        options = [
            {"name": "pep", "type": "infile", "format": "sequence.fasta"},
            {"name": "bed", "type": "outfile", "format": "gene_structure.bed"},  # 输出结果
            {"name": "cds", "type": "outfile", "format": "sequence.fasta"},  # 输出结果
            {"name": "e_value", "type": "float", "default": 1e-3},
            {"name": "pfam.domtblout", "type": "outfile", "format": "denovo_rna_v2.common"},
            {"name": "pfam_version", "type": "string", "default": "32"},
            # 输出结果
        ]
        self.add_option(options)
        self.step.add_steps('hmmscan')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

# self.step.xx.{finish, start}和self.step.add_steps("xx),里面的xx需要保持一致
    def step_start(self):
        self.step.hmmscan.start()
        self.step.update()

    def step_end(self):
        self.step.hmmscan.finish()
        self.step.update()

    def check_options(self):
        """
        检查参数
        """
        if not self.option("pep").is_set:
            raise OptionError("请传入fasta序列文件", code = "35000901")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 3
        self._memory = '10G'
    # 这个end就是下面run里面的self.end()调用的end函数
    def end(self):
        # result_dir = self.add_upload_dir(self.output_dir)
        # result_dir.add_relpath_rules([
        #     [".", "", "结果输出目录"],
        #     # ["./estimators.xls", "xls", "alpha多样性指数表"]
        # ])
        # result_dir.add_regexp_rules([
        #     [r"transdecoder.pep$", "fasta", "蛋白质序列文件"],
        #     [r"transdecoder.cds$", "fasta", "cds序列文件"],
        #     [r"transdecoder.bed$", "bed", "orf位置信息bed格式文件"]
        # ])
        #
        # if self.option("species_type").lower() == "animal":
        #     result_dir.add_regexp_rules([
        #         ["./animal_transcription_analysis_details", "", "动物转录因子分析详情"]
        #     ])
        #
        # if self.option("species_type").lower() == "plant":
        #     result_dir.add_regexp_rules([
        #         ["./plant_transcription_analysis_details", "", "植物转录因子分析详情"]
        #     ])

        super(HmmAgent, self).end()


class HmmTool(Tool):
    """
    version 1.0
    """

    def __init__(self, config):
        super(HmmTool, self).__init__(config)
        self.transdecoder_path = "bioinfo/gene-structure/TransDecoder-3.0.0/"
        # self.hmmscan_path = "bioinfo/align/hmmer-3.1b2-linux-intel-x86_64/binaries/"
        self.hmmscan_path = "miniconda2/bin/"
        # self.pfam_db = self.config.SOFTWARE_DIR + "/database/Pfam/Pfam-A.hmm"# hmm参考库
        self.pfam_db = AnnotConfig().get_file_path(
            file ="Pfam-A",
            db = "pfam",
            version = self.option("pfam_version"))
        # self.pfam_db = self.config.SOFTWARE_DIR + "/database/Annotation/other2019/pfam32/Pfam-A.hmm"  # hmm参考库
        # self.fasta_name = self.option("fasta").prop["path"].split("/")[-1]# 得到fasta文件名
        self.gcc = self.config.SOFTWARE_DIR + '/gcc/5.1.0/bin'
        self.gcc_lib = self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64'
        # TransDecoder-3.0.0的运行依赖perl的环境，所以程序必须设置一个环境变量，直接在我们的服务器程序里测试单条是不行的，我们没有perl_env
        self.perl_lib = self.config.SOFTWARE_DIR + '/miniconda2/bin'
        self.set_environ(PATH=self.perl_lib)
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)
        self.tfdb_path = self.config.SOFTWARE_DIR + "/database/TFDB/"

    def hmmscan(self, pep):
        cmd = "{}hmmscan -E {} --domE {} --cpu 2 --noali --acc --notextw --domtblout {} --tblout {} {} {}".format(
            self.hmmscan_path, self.option("e_value"), self.option("e_value"), "pfam.domtblout", "pfam.tblout", self.pfam_db, pep)
        print(cmd)
        self.logger.info("开始运行hmmscan")
        self.logger.info(cmd)
        command = self.add_command("hmmscan", cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("运行hmmscan结束")
        else:
            self.set_error("运行hmmscan出错")
        pfam_out("pfam.domtblout")

    def run(self):
        super(HmmTool, self).run()
        """将Blast和Pfam搜索结果整合到编码区域选择
        TransDecoder借助上面生成的输出结果来确定将这些被blast命中的和结构域命中的多肽保留在报告的编码区集合中。像这样运行TransDecoder.Predict：
        TransDecoder.Predict -t target_transcripts.fasta --retain_pfam_hits pfam.domtblout --retain_blastp_hits blastp.outfmt6
        最终的编码区预测结果将包含与编码区域一致的序列字符以及blast得到的直系同源结果或pfam结构域的内容。这里选择是否比对pfam库,
        这一步是为了做注释使用，hmmscan其实做转录因子其实是不需要做domtblout的结果输出的，我们流程默认的Predict是要整合pfam数据库的结果"""

        with open(self.option("pep").prop['path'], "r") as fr:
            content = fr.read()
            if not content:
                self.logger.info("没有找到orf")
                self.end()
        # 这里如果直接写self.hmmscan(self.option(
        # "pep"))会报出命令里面有重定向符，因为pep文件里面有>，当然你肯定需要路径啊
        self.hmmscan(self.option("pep").prop['path'])
        self.end()

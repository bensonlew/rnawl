# -*- coding: utf-8 -*-
# __author__ = 'qindanhua'

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
import shutil
from biocluster.core.exceptions import OptionError
from mbio.packages.denovo_rna.gene_structure.pfam_domtblout import pfam_out
from mbio.packages.rna.annot_config import AnnotConfig

class OrfAgent(Agent):
    """
    transdecoder:orf预测软件
    hmmscan：Pfam数据库比对工具
    version 1.0
    author: qindanhua
    last_modify: 2016.07.11
    """

    def __init__(self, parent):
        super(OrfAgent, self).__init__(parent)
        options = [
            {"name": "fasta", "type": "infile", "format": "sequence.fasta"},  # 输入文件
            {"name": "search_pfam", "type": "bool", "default": True},  # 是否比对Pfam数据库
            {"name": "p_length", "type": "int", "default": 50},  # 最小蛋白长度
            {"name": "Markov_length", "type": "int", "default": 3000},  # 马尔科夫训练长度
            {"name": "bed", "type": "outfile", "format": "gene_structure.bed"},  # 输出结果
            {"name": "cds", "type": "outfile", "format": "sequence.fasta"},  # 输出结果
            {"name": "pep", "type": "outfile", "format": "sequence.fasta"},  # 输出结果
            {"name": "pfam_version", "type": "string", "default": "32"},
        ]
        self.add_option(options)
        self.step.add_steps('orf')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.orf.start()
        self.step.update()

    def step_end(self):
        self.step.orf.finish()
        self.step.update()

    def check_options(self):
        """
        检查参数
        """
        if not self.option("fasta").is_set:
            raise OptionError("请传入fasta序列文件")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 10
        self._memory = '5G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
            # ["./estimators.xls", "xls", "alpha多样性指数表"]
        ])
        result_dir.add_regexp_rules([
            [r"transdecoder.pep$", "fasta", "蛋白质序列文件"],
            [r"transdecoder.cds$", "fasta", "cds序列文件"],
            [r"transdecoder.bed$", "bed", "orf位置信息bed格式文件"]
        ])
        if self.option("search_pfam") is True:
            result_dir.add_relpath_rules([
                ["./pfam_domain", "", "Pfam比对蛋白域结果信息"]
            ])
        # print self.get_upload_files()
        super(OrfAgent, self).end()


class OrfTool(Tool):
    """
    version 1.0
    """

    def __init__(self, config):
        super(OrfTool, self).__init__(config)
        self.transdecoder_path = "bioinfo/gene-structure/TransDecoder-3.0.0/"
        self.hmmscan_path = "bioinfo/align/hmmer-3.1b2-linux-intel-x86_64/binaries/"
        self.pfam_db = AnnotConfig().get_file_path(
            file ="Pfam-A",
            db = "pfam",
            version = self.option("pfam_version"))
        # self.pfam_db = self.config.SOFTWARE_DIR + "/database/Pfam/Pfam-A.hmm"
        self.fasta_name = self.option("fasta").prop["path"].split("/")[-1]
        self.gcc = self.config.SOFTWARE_DIR + '/gcc/5.1.0/bin'
        self.gcc_lib = self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64'
        self.perl_lib = self.config.SOFTWARE_DIR + '/miniconda2/bin'
        self.set_environ(PATH=self.perl_lib)
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)

    def td_longorfs(self):
        self.logger.info(self.option("p_length"))
        cmd = "{}TransDecoder.LongOrfs -t {} -m {}".format(self.transdecoder_path, self.option("fasta").prop["path"],
                                                           self.option("p_length"))
        print(cmd)
        self.logger.info("开始提取长orf")
        command = self.add_command("transdecoder_longorfs", cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("提取结束！")
        else:
            self.set_error("提取orf过程出错")

    def td_predict(self, hmm_out=None):
        if hmm_out is None:
            hmm_out = ""
        else:
            hmm_out = "--retain_pfam_hits {}".format(hmm_out)
        cmd = "{}TransDecoder.Predict -t {} -T {} {}".format(
            self.transdecoder_path, self.option("fasta").prop["path"], self.option("Markov_length"), hmm_out)
        print(cmd)
        self.logger.info("开始预测编码区域")
        command = self.add_command("transdecoder_predict", cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("预测编码区域完成！")
        else:
            self.set_error("预测过程过程出错")

    def hmmscan(self, pep):
        cmd = "{}hmmscan --cpu 20 --noali --cut_nc --acc --notextw --domtblout {} {} {}".format(
            self.hmmscan_path, "pfam.domtblout", self.pfam_db, pep)
        print(cmd)
        self.logger.info("开始运行hmmscan")
        command = self.add_command("hmmscan", cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("运行hmmscan结束")
        else:
            self.set_error("运行hmmscan出错")
        pfam_out("pfam.domtblout")

    def set_output(self):
        self.logger.info("set out put")
        for f in os.listdir(self.output_dir):
            os.remove(os.path.join(self.output_dir, f))
        pep = '{}.transdecoder.pep'.format(self.fasta_name)
        bed = '{}.transdecoder.bed'.format(self.fasta_name)
        cds = '{}.transdecoder.cds'.format(self.fasta_name)
        os.link(self.work_dir+"/"+pep, self.output_dir+"/"+pep)
        self.option('pep').set_path(self.output_dir+"/"+pep)
        os.link(self.work_dir+"/"+bed, self.output_dir+"/"+bed)
        self.option('bed').set_path(self.output_dir+"/"+bed)
        os.link(self.work_dir+"/"+cds, self.output_dir+"/"+cds)
        self.option('cds').set_path(self.output_dir+"/"+cds)
        if self.option("search_pfam") is True:
            os.link(self.work_dir+"/"+"pfam_domain", self.output_dir+"/"+"pfam_domain")
        fasta_for_len = os.path.join(self.work_dir, "ORF_fasta")
        if os.path.exists(fasta_for_len):
            shutil.rmtree(fasta_for_len)
        os.makedirs(fasta_for_len)
        os.link(self.work_dir+"/"+cds, fasta_for_len+"/"+cds+".fasta")
        self.logger.info("done")

    def run(self):
        """
        运行
        """
        super(OrfTool, self).run()
        if self.option("search_pfam") is True:
            self.td_longorfs()
            with open(self.work_dir + "/{}.transdecoder_dir/longest_orfs.pep".format(self.fasta_name), "r") as fr:
                content = fr.read()
                if not content:
                    self.logger.info("没有找到orf")
                    self.end()
            self.hmmscan("{}.transdecoder_dir/longest_orfs.pep".format(self.fasta_name))
            self.td_predict("pfam.domtblout")
        else:
            self.td_longorfs()
            # pfam_out("pfam.domtblout")
            # self.set_output()
            self.td_predict()
        self.set_output()
        self.end()

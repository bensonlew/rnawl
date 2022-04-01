# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

import os,re
from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool
from Bio import SeqIO

class UnicyclerAgent(Agent):
    """
    进行Unicycler拼接,主要使用二代数据，三代数据组装，有利于完成图组装（包含环型和线型基因组）
    version: v1.0
    author: gaohao
    last_modify: 2021.05.18
    """
    def __init__(self, parent):
        super(UnicyclerAgent, self).__init__(parent)
        options = [
            {"name": "read1", "type": "infile", "format": "sequence.fastq"},  # 输入fastq文件
            {"name": "read2", "type": "infile", "format": "sequence.fastq"},  # 输入fastq文件
            {"name": "unpaired", "type": "infile", "format": "sequence.fastq"},  #双端fq的单端reads
            {"name": "long", "type": "infile", "format": "sequence.fastq"},  # 三代数据只支持fq格式
            {"name": "min_len", "type": "int", "default": 200},  # 拼接最小长度100bp
            {"name": "sample_name", "type": "string", "required": True},
            {"name": "data_type", "type": "string"}
        ]
        self.add_option(options)
        self.queue = 'SANGER'
        self._memory_increase_step = 50  # 每次重运行增加内存50G

    def check_options(self):
        """
        检查参数
        :return:
        """
        pass

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 8
        self._memory = "100G"

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(UnicyclerAgent, self).end()


class UnicyclerTool(Tool):
    def __init__(self, config):
        super(UnicyclerTool, self).__init__(config)
        self.lib = self.config.SOFTWARE_DIR + '/program/Python35/lib:' + self.config.SOFTWARE_DIR + '/gcc/7.2.0/lib64'
        self.path = self.config.SOFTWARE_DIR + '/bioinfo/align/bowtie2-2.4.4:' + \
                    self.config.SOFTWARE_DIR + '/program/Python35/bin:' + self.config.SOFTWARE_DIR + '/gcc/7.2.0/bin:' + \
                    self.config.SOFTWARE_DIR + '/bioinfo/Genomic/Sofware/racon/build/bin:' + \
                    self.config.SOFTWARE_DIR + '/bioinfo/align/bcftools-1.6:' + \
                    self.config.SOFTWARE_DIR + '/bioinfo/metaGenomic/SPAdes-3.15.3-Linux/bin:' + \
                    self.config.SOFTWARE_DIR + '/bioinfo/align/ncbi-blast-2.12.0+/bin:' + \
                    self.config.SOFTWARE_DIR + '/bioinfo/Genomic/Sofware/Pilon1.23:' + \
                    self.config.SOFTWARE_DIR + '/bioinfo/align/samtools-1.12'
        self.set_environ(PATH=self.path, LD_LIBRARY_PATH=self.lib)
        self.set_environ(JAVA_HOME=self.config.SOFTWARE_DIR + '/program/sun_jdk1.8.0')
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/program/sun_jdk1.8.0/bin')
        self.set_environ(CLASSPATH='.:' + self.config.SOFTWARE_DIR + '/program/sun_jdk1.8.0/lib/dt.jar:' + self.config.SOFTWARE_DIR + '/program/sun_jdk1.8.0/lib/tools.jar')
        self.ass_path = '/program/Python35/bin/unicycler'

    def unicycler_run(self):
        """
        进行Unicycler拼接
        :return:
        """
        cmd = "{} ".format(self.ass_path)
        if self.option("read1").is_set:
            cmd += "-1 {} -2 {} ".format(self.option("read1").prop['path'], self.option('read2').prop["path"])
        if self.option("unpaired").is_set:
            cmd += "-s {} ".format(self.option("unpaired").prop['path'])
        if self.option("long").is_set:
            cmd += "-l {} ".format(self.option("long").prop['path'])
        cmd += "-o {}".format(self.work_dir + "/result")
        self.logger.info("运行Unicycler拼接")
        command = self.add_command("unicycler", cmd, ignore_error=True)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行Unicycler完成")
        else:
            self.logger.info("return code: %s" % command.return_code)
            self.set_error("Unicycler运行出错!")

    def get_status(self):
        with open(self.output_dir+"/all.status.txt", "w") as g:
            list1 = []
            for fa_iterator in SeqIO.parse(self.work_dir + "/result/assembly.fasta", "fasta"):
                id = fa_iterator.id
                if re.search("circular=true", fa_iterator.description):
                    g.write("{}\t{}\n".format("Scaffold" + str(id), "circle"))
                else:
                    g.write("{}\t{}\n".format("Scaffold" + str(id), "linear"))
                fa_iterator.id = "Scaffold" + str(id)
                list1.append(fa_iterator)
            SeqIO.write(list1, self.output_dir + "/" + self.option("sample_name") + ".scaffold.fna", "fasta")


    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        if os.path.exists(self.output_dir+"/"+ self.option("sample_name") + ".raw.assembly.log"):
            os.remove(self.output_dir+"/"+ self.option("sample_name") + ".raw.assembly.log")
        os.link(self.work_dir+"/result/unicycler.log", self.output_dir+"/"+ self.option("sample_name") + ".raw.assembly.log")
        self.logger.info("设置Unicycler分析结果目录成功")

    def run(self):
        """
        运行
        :return:
        """
        super(UnicyclerTool, self).run()
        self.unicycler_run()
        self.get_status()
        self.set_output()
        self.end()
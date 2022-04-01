# -*- coding: utf-8 -*-
# __author__ = 'sj'
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
import json


class TophatAgent(Agent):
    """
    tophat  
    version 2.0
    author: sj
    last_modify: 2016.9.8
    """
    def __init__(self, parent):
        super(TophatAgent, self).__init__(parent)
        options = [
            {"name": "ref_genome", "type": "string"},
            {"name": "ref_genome_custom", "type": "infile", "format": "sequence.fasta"},
            {"name": "mapping_method", "type": "string"},
            {"name": "seq_method", "type": "string"},
            {"name": "single_end_reads", "type": "infile", "format": "sequence.fastq"},
            {"name": "left_reads", "type": "infile", "format": "sequence.fastq"},
            {"name": "right_reads", "type": "infile", "format": "sequence.fastq"},
            {"name": "bam_output", "type": "outfile", "format": "align.bwa.bam"},
            {"name": "assemble_method", "type": "string", "default": "none"},
            {"name": "sample", "type": "string"},
            {"name": "mate_std", "type": "int", "default": 50},  # 末端配对插入片段长度标准差
            {"name": "mid_dis", "type": "int", "default": 50},  # 两个成对引物间的距离中间值
            {"name": "result_reserved", "type": "int", "default": 1},  # 最多保留的比对结果数目
            {"name": "strand_specific", "type": "bool", "default": "False"}
            ]
        self.add_option(options)
        self.step.add_steps('Tophat')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.Tophat.start()
        self.step.update()

    def step_end(self):
        self.step.Tophat.finish()
        self.step.update()

    def check_options(self):
        self.logger.info("开始运行check步骤")
        if self.option("ref_genome") == "customer_mode" and not self.option("ref_genome_custom").is_set:
            raise OptionError("请上传自定义参考基因组")
        if self.option("seq_method") == "PE":
            if self.option("single_end_reads").is_set:
                raise OptionError("上传的是单端测序的序列，请上传双端序列")
            elif not (self.option("left_reads").is_set and self.option("right_reads").is_set):
                raise OptionError("缺少某端序列")
            else:
                pass
        else:
            if not self.option("single_end_reads").is_set:
                raise OptionError("请上传单端序列")
            elif self.option("left_reads").is_set or self.option("right_reads").is_set:
                raise OptionError("只需要有单端的序列")
            else:
                pass
        if not self.option("assemble_method") in ["cufflinks", "stringtie", "none"]:
            raise OptionError("请选择拼接软件")
        return True

    def set_resource(self):
        self._cpu = 10
        self._memory = '10G'

    def end(self):
        super(TophatAgent, self).end()


class TophatTool(Tool):
    def __init__(self, config):
        super(TophatTool, self).__init__(config)
        self.cmd_path = "bioinfo/align/tophat-2.1.1/tophat-2.1.1.Linux_x86_64/"
    
    def run_build_index(self):
        ref_name = os.path.basename(self.option("ref_genome_custom").prop['path'])
        new_ref_genome_custom_path = os.path.join(self.work_dir, ref_name)
        if os.path.exists(new_ref_genome_custom_path):
            os.remove(new_ref_genome_custom_path)
        os.link(self.option("ref_genome_custom").prop['path'], new_ref_genome_custom_path)
        ref_genome_custom_path = os.path.split(new_ref_genome_custom_path)[0]
        ref_path = os.path.join(ref_genome_custom_path, "ref_index")  # 保证生成的索引都在和基因组同一级的目录下
        cmd = "{}bowtie2-build {} {}".format(self.cmd_path, new_ref_genome_custom_path, ref_path)
        index_command_obj = self.add_command("build_index", cmd).run()
        self.wait(index_command_obj)
        if index_command_obj.return_code == 0:
            self.logger.info("索引建立完成")
            self.run_tophat(ref_path)
        else:
            self.set_error("索引建立出错")
            raise Exception("索引建立出错")

    def run_tophat(self, index_ref):
        self.logger.info(self.option("sample"))
        name = self.option("sample")
        if name.find("_sickle") != -1:
            name = name[:-9]
        cmd = "{}tophat2  --mate-std-dev {} -r {} -g {}".format(self.cmd_path, self.option("mate_std"),
                                                                self.option("mid_dis"),
                                                                self.option("result_reserved")
                                                                )
        if self.option("seq_method") == "PE":
            cmd += " {} {} {}".format(index_ref, self.option("left_reads").prop['path'],
                                      self.option("right_reads").prop['path'])
        else:
            cmd += " {} {}".format(index_ref, self.option("single_end_reads").prop['path'])
        if self.option("strand_specific"):
            cmd += " --library-type fr-firststrand"
        tophat_command = self.add_command("tophat", cmd)
        self.logger.info("开始运行tophat")
        tophat_command.run()
        self.wait()
        if tophat_command.return_code == 0:
            output = os.path.join(self.work_dir, "tophat_out/accepted_hits.bam")
            outfile_path = os.path.split(output)[0]
            outfile = os.path.join(outfile_path, name + ".bam")
            os.rename(output, outfile)
            if os.path.exists(self.output_dir + "/" + name + ".bam"):
                os.remove(self.output_dir + "/" + name + ".bam")
            os.link(outfile, self.output_dir + "/" + name + ".bam")
        return True
    
    def run(self):
        """
        运行
        :return:
        """
        super(TophatTool, self).run()
        if self.option("ref_genome") == "customer_mode":
            self.logger.info("开始运行自定义模式")
            self.run_build_index()
        else:
            with open(self.config.SOFTWARE_DIR + "/database/Genome_DB_finish/annot_species.json", "r") as f:
                dict = json.loads(f.read())
                rel_index = dict[self.option("ref_genome")]["dna_index"]
                abs_index = self.config.SOFTWARE_DIR + "/database/Genome_DB_finish/" +  rel_index
                # index_ref = os.path.join(os.path.split(ref)[0], "ref_index")
                self.run_tophat(abs_index)
        self.end()

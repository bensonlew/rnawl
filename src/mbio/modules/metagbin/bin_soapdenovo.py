# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'
import os
import re
import time
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
from collections import defaultdict
from mbio.packages.metagbin.common_function import link_dir


class BinSoapdenovoModule(Module):
    """
    微生物宏基因组binning
    """
    def __init__(self, work_id):
        super(BinSoapdenovoModule, self).__init__(work_id)
        options = [
            {"name": "fastq1", "type": "infile", "format": "sequence.fastq"},  # 输入文件,sample.sickle.l.fastq
            {"name": "fastq2", "type": "infile", "format": "sequence.fastq"},  # 输入文件,sample.sickle.r.fastq
            {"name": "fastqs", "type": "infile", "format": "sequence.fastq"},  # 输入文件,sample.sickle.s.fastq
            {'name': 'sample_name', "type": "string"},  # 基因组名称对应项目的genome名称
            {"name": "max_rd_len", "type": "string"},  # read最大读长
            {"name": "insert_size", "type": "string"},  # 平均插入片段长度
            {"name": "reverse_seq", "type": "string", "default": "0"},  # 配置文件的其他参数
            {"name": "asm_flags", "type": "string", "default": "3"},  # 配置文件的其他参数
            {"name": "rank", "type": "string", "default": "1"},  # 配置文件的其他参数
            {"name": "scaffold", "type": "outfile", "format": "sequence.fasta"},  # 输出文件
        ]
        self.add_option(options)
        self.step.add_steps('assemble','scaf_select','gapcloser', 'scaf_agp_contig')
        self.scaf_select = self.add_tool('assemble.scaf_select')
        self.gapcloser = self.add_tool('assemble.gapcloser_scaf_bin')
        self.scaf_agp_contig = self.add_tool('assemble.scaf_agp_contig')
        self.sample_path = {}
        self.modules = []
        self.assembly_result_path = ''

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option('fastq1').is_set:
            raise OptionError('必须输入fastq1序列', code="")
        if not self.option('fastq2'):
            raise OptionError('必须输入fastq2序列', code="")
        if not self.option('fastqs'):
            raise OptionError('必须输入fastqs序列', code="")

    def finish_update(self, event):
        step = getattr(self.step, event['data'])
        step.finish()
        self.step.update()

    def run_config(self):
        """
        生成config文件
        :return:
        """
        config_file = self.work_dir + "/" + self.option('sample_name')+ ".config"
        with open(config_file, "w+") as fw:
            first = "max_rd_len=" + self.option('max_rd_len') + "\n"
            second = "[LIB]" + "\n"
            third = "avg_ins=" + self.option('insert_size') + "\n"
            forth = "reverse_seq=" + self.option('reverse_seq') + "\n"
            fifth = "asm_flags=" + self.option('asm_flags') + "\n"
            sixth = "rank=" + self.option('rank') + "\n"
            q1 = "q1=" + self.option('fastq1').prop['path'] + "\n"
            q2 = "q2=" + self.option('fastq2').prop['path']
            if not self.option('fastqs').is_set:
                fw.write(first + second + third + forth + fifth + sixth + q1 + q2)
            else:
                qs = "\nq=" + self.option('fastqs').prop['path']
                fw.write(first + second + third + forth + fifth + sixth + q1 + q2 + qs)
        #self.config_file.run()
        #self.step.config_file.finish()
        self.step.assemble.start()
        self.step.update()

    def run_assemble(self):
        """
        用soapdenovo进行组装
        :return:
        """
        num = self.get_size()
        self.assembly_result_path = os.path.join(self.work_dir, "Assembly")
        if os.path.exists(self.assembly_result_path):
            pass
        else:
            os.mkdir(self.assembly_result_path)
        n = 0
        for k in [21,23,25,27,29,31,33,35,37,39,41]:
            for d in [3,5,10]:
                self.assemble = self.add_tool('assemble.bac_denova_ass')
                config_file = self.work_dir + "/" + self.option('sample_name')+ ".config"
                self.step.add_steps('assemble{}'.format(n))
                opts = {
                    "config": config_file,
                    "kmerFreqCutoff": int(d),
                    "sample_name": self.option('sample_name'),
                    "kmer": int(k),
                    "mem": num
                }
                self.assemble.set_options(opts)
                step = getattr(self.step, 'assemble{}'.format(n))
                step.start()
                self.step.update()
                self.assemble.on('end', self.finish_update, 'assemble{}'.format(n))
                self.modules.append(self.assemble)
                n += 1
        self.logger.info(self.modules)
        self.on_rely(self.modules, self.run_select)
        self.step.update()
        for module in self.modules:
            module.run()

    def run_select(self):
        """
        对soapdenovo组装的结果根据N50筛选
        :return:
        """
        assemble_dir = self.assembly_result_path
        for i in self.modules:
            for f in os.listdir(i.output_dir):
                file_path = os.path.join(i.output_dir, f)
                new_path = os.path.join(self.assembly_result_path, os.path.basename(file_path))
                if os.path.exists(new_path):
                    os.remove(new_path)
                os.link(file_path, new_path)
        opts = {
            "seq_dir": assemble_dir,
        }
        self.scaf_select.set_options(opts)
        self.scaf_select.run()
        self.step.scaf_select.finish()
        self.step.gapcloser.start()
        self.step.update()

    def run_gaploser(self):
        """
        先过滤掉小于300bp的序列，再对筛选结果根据参考序列进行补洞
        :return:
        """
        config_file = self.work_dir + "/" + self.option('sample_name')+ ".config"
        opts = {
            "seq_scaf": self.scaf_select.option('scf_seq'),
            "config": config_file,
            "sample_name": self.option('sample_name'),
        }
        self.gapcloser.set_options(opts)
        self.gapcloser.on('end', self.set_output, 'gapcloser')
        self.gapcloser.run()
        self.step.gapcloser.finish()
        self.step.update()

    def run_scaf_agp_contig(self):
        """
        过滤掉长度小于200bp的序列，挑选最终的scaffold和contigs
        :return:
        """
        opts = {
            "seq_scaf": self.gapcloser.option('seq'),
            "sample_name": self.option('sample_name'),
        }
        self.scaf_agp_contig.set_options(opts)
        self.scaf_agp_contig.on('end', self.set_output, 'scaf_agp_contig')
        self.scaf_agp_contig.run()
        self.step.scaf_agp_contig.finish()
        self.step.update()

    def run(self):
        """
        运行
        :return:
        """
        super(BinSoapdenovoModule, self).run()
        self.run_config()
        self.run_assemble()
        self.scaf_select.on('end', self.run_gaploser)
        self.gapcloser.on('end', self.run_scaf_agp_contig)
        self.scaf_agp_contig.on('end', self.end)

    def set_output(self, event):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        if event["data"] == "gapcloser":
            self.logger.info('gapcloser运行成功')
            #link_dir(self.gapcloser.output_dir,self.output_dir + '/scf')
            #self.option('scaffold').set_path(self.output_dir + '/scf/' + self.option('sample_name') + '.scaffold.fna')
        if event["data"] == "scaf_agp_contig":
            if os.path.exists(self.output_dir +'/'+self.option('sample_name')+'_scaffold.fa'):
                os.remove(self.output_dir +'/'+self.option('sample_name')+'_scaffold.fa')
            os.link(self.scaf_agp_contig.option('scaffold').prop['path'] , self.output_dir + '/'+self.option('sample_name')+'_scaffold.fa')
            #link_dir(self.scaf_agp_contig.output_dir , self.output_dir + '/assembly')
            self.option('scaffold',self.output_dir + '/'+self.option('sample_name')+'_scaffold.fa')
        self.logger.info('设置结果目录成功')

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(BinSoapdenovoModule, self).end()

    def get_size(self):
        r1_bases =os.path.getsize(self.option('fastq1').prop['path'])
        rs_bases = os.path.getsize(self.option('fastqs').prop['path'])
        num = float(r1_bases*2+rs_bases)/1024/1024/1024
        return num
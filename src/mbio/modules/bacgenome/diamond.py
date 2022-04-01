# -*- coding: utf-8 -*-
# __author__ = 'shijin'
# last_modify:2017.03.17

from biocluster.module import Module
import os
import shutil
from biocluster.core.exceptions import OptionError
from mbio.packages.bac_comp_genome.common_function import link_dir,link_file


class DiamondModule(Module):
    def __init__(self, work_id):
        super(DiamondModule, self).__init__(work_id)
        self._fasta_type = {'Protein': 'prot', 'DNA': 'nucl'}
        options = [
            {"name": "query", "type": "infile", "format": "sequence.fasta"},  # 输入文件
            {"name": "lines", "type": "int", "default": 500},  # 将fasta序列拆分此行数的多个文件
            {"name": "query_type", "type": "string"},  # 输入的查询序列的格式，为nucl或者prot
            {"name": "database", "type": "string", "default": "nr"},
            # 比对数据库 nt nr string swissprot kegg customer_mode
            {"name": "outfmt", "type": "int", "default": 6},  # 输出格式，只为6
            {"name": "blast", "type": "string"},
            {"name": "reference", "type": "infile", "format": "sequence.fasta"},  # 参考序列  选择customer时启用
            {"name": "reference_type", "type": "string"},  # 参考序列(库)的类型  为nucl或者prot
            {"name": "evalue", "type": "float", "default": 1e-5},  # evalue值
            {"name": "num_threads", "type": "int", "default": 9},  # cpu数
            {"name": "outtable", "type": "outfile", "format": "align.blast.blast_table"},  # 输出格式为6时输出
            # 当输出格式为非5，6时，只产生文件不作为outfile
            {"name": "sensitive", "type": "int", "default": 2}
        ]
        self.add_option(options)
        self.splitfasta = self.add_tool("sequence.split_fasta")
        self.step.add_steps('blast', 'split_fasta')
        self.blast_tools = []
        self.catblast_tools = []

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def check_options(self):
        if not self.option("query").is_set:
            raise OptionError("必须设置参数query", code="21100201")
        if self.option('query_type') not in ['nucl', 'prot']:
            raise OptionError('query_type查询序列的类型为nucl(核酸)或者prot(蛋白):%s', variables=(self.option('query_type')), code="21100202")
        if self.option('outfmt') not in [6]:
            raise OptionError('outfmt遵循blast+输出规则，目前只支持6输出格式：%s', variables=(self.option('outfmt')), code="21100203")
        if not 1 > self.option('evalue') >= 0:
            raise OptionError('E-value值设定必须为[0-1)之间：%s' , variables=(self.option('evalue')), code="21100204")
        if not isinstance(self.option('lines'), int):
            raise OptionError("行数必须为整数", code="21100205")
        return True

    def run_splitfasta(self):
        self.splitfasta.set_options({
            "fasta": self.option("query"),
            "lines": self.option("lines"),
        })
        self.splitfasta.on('start', self.set_step, {'start': self.step.split_fasta})
        self.splitfasta.on('end', self.set_step, {'end': self.step.split_fasta})
        self.splitfasta.on('end', self.set_step, {'start': self.step.blast})
        self.splitfasta.on('end', self.run_blast)
        self.splitfasta.run()

    def run_blast(self):
        opts = {
            "query_type": self.option("query_type"),
            "database": self.option("database"),
            "outfmt": 6,
            "blast": self.option("blast"),
            "evalue": self.option("evalue"),
            "num_threads": self.option("num_threads"),
            "sensitive": self.option("sensitive")
        }
        for f in os.listdir(self.splitfasta.output_dir):
            opts['query'] = os.path.join(self.splitfasta.output_dir, f)
            blast_tool = self.add_tool('bacgenome.diamond')
            blast_tool.set_options(opts)
            self.blast_tools.append(blast_tool)
        if len(self.blast_tools) == 1:
            self.blast_tools[0].on("end", self.set_output)
        else:
            self.on_rely(self.blast_tools, self.set_output)
        for tool in self.blast_tools:
            tool.run()

    def set_output(self):
        if len(self.blast_tools) == 1:
            for i in os.listdir(self.catblast_tools[0].output_dir):
                link_file(self.catblast_tools[0].output_dir+"/"+i, self.output_dir + "/blast.table")
            self.option('outtable').set_path(self.output_dir + "/blast.table")
        else:
            list1 = []
            for module in self.blast_tools:
                for i in os.listdir(module.output_dir):
                    list1.append(module.output_dir+"/"+i)
            self.logger.info(" ".join(list1))
            os.system("cat {} > {}".format(" ".join(list1), self.output_dir + "/blast.table"))
            self.option('outtable').set_path(self.output_dir + "/blast.table")
        self.end()

    def run(self):
        super(DiamondModule, self).run()
        self.run_splitfasta()

    def end(self):
        repaths = [
            [".", "", "blast输出目录"],
            ["blast.table.", "xls", "blast xls输出结果文件"]
        ]
        sdir = self.add_upload_dir(self.output_dir)
        sdir.add_relpath_rules(repaths)
        super(DiamondModule, self).end()
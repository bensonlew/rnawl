# -*- coding: utf-8 -*-
# __author__ = 'shijin'
# last_modify:2017.09.11
# last_modify by shaohua.yuan

from biocluster.module import Module
import os
import shutil
from biocluster.core.exceptions import OptionError


class MetaDiamondModule(Module):
    def __init__(self, work_id):
        super(MetaDiamondModule, self).__init__(work_id)
        self._fasta_type = {'Protein': 'prot', 'DNA': 'nucl'}
        options = [
            {"name": "query", "type": "infile", "format": "sequence.fasta"},  # 输入文件
            {"name": "lines", "type": "int", "default": 50000},  # 将fasta序列拆分此行数的多个文件
            {"name": "query_type", "type": "string"},  # 输入的查询序列的格式，为nucl或者prot
            {"name": "database", "type": "string", "default": "nr"},
            # 比对数据库 nt nr string swissprot kegg customer_mode
            {"name": "outfmt", "type": "int", "default": 5},  # 输出格式，只为5
            {"name": "blast", "type": "string"},
            {"name": "reference", "type": "infile", "format": "sequence.fasta"},  # 参考序列  选择customer时启用
            {"name": "reference_type", "type": "string"},  # 参考序列(库)的类型  为nucl或者prot
            {"name": "evalue", "type": "float", "default": 1e-5},  # evalue值
            {"name": "num_threads", "type": "int", "default": 10},  # cpu数
            {"name": "outxml", "type": "outfile", "format": "align.blast.blast_xml"},  # 输出格式为5时输出
            {"name": "outxml_dir", "type": "outfile", "format": "align.blast.blast_xml_dir"},    #输出拆分文件比对的xml文件夹
            {"name": "outtable", "type": "outfile", "format": "align.blast.blast_table"},  # 输出格式为6时输出
            {"name": "target_num", "type": "int", "default": 1}         ### maximum number of target sequences to report alignments for
            # 当输出格式为非5，6时，只产生文件不作为outfile
        ]
        self.add_option(options)
        self.catblast = self.add_tool("align.cat_blastout")
        self.splitfasta = self.add_tool("sequence.split_fasta")
        self.step.add_steps('blast', 'split_fasta', 'cat_blastout')
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
            raise OptionError("必须设置参数query", code="21100501")
        if self.option('query_type') not in ['nucl', 'prot']:
            raise OptionError('query_type查询序列的类型为nucl(核酸)或者prot(蛋白):%s', variables=(self.option('query_type')), code="21100502")
        if self.option('outfmt') not in [5, 6]:
            raise OptionError('outfmt遵循blast+输出规则，目前只支持5，6输出格式：%s', variables=(self.option('outfmt')), code="21100503")
        if not 1 > self.option('evalue') >= 0:
            raise OptionError('E-value值设定必须为[0-1)之间：%s', variables=(self.option('evalue')), code="21100504")
        if not isinstance(self.option('lines'), int):
            raise OptionError("行数必须为整数", code="21100505")
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
            "outfmt": 5,
            "blast": self.option("blast"),
            "evalue": self.option("evalue"),
            "num_threads": self.option("num_threads"),
            "target_num": self.option("target_num"),
        }
        #if self.option("target_num") in self.get_option_object().keys():
            #opts["target_num"] = self.option("target_num")
        for f in os.listdir(self.splitfasta.output_dir):
            opts['query'] = os.path.join(self.splitfasta.output_dir, f)
            blast_tool = self.add_tool('align.meta_diamond')
            blast_tool.set_options(opts)
            # blast_tool.run()
            self.blast_tools.append(blast_tool)
        if len(self.blast_tools) == 1:
            self.blast_tools[0].on("end", self.run_catblastout)
        else:
            self.on_rely(self.blast_tools, self.run_catblastout)
        for tool in self.blast_tools:  # edited by shijin on 20170627
            tool.run()

    def run_catblastout(self):
        self.set_step(event={'data': {'end': self.step.blast}})
        if os.path.exists(self.work_dir + '/blast_tmp'):
            shutil.rmtree(self.work_dir + '/blast_tmp')
        os.mkdir(self.work_dir + '/blast_tmp')
        for i in self.blast_tools:
            if len(os.listdir(i.output_dir)) != 1:
                self.logger.info(str(os.listdir(i.output_dir)))
                for f in os.listdir(i.output_dir):
                    if f.endswith("xml"):
                        file = f
            else:
                file = os.listdir(i.output_dir)[0]
            _path = os.path.join(i.output_dir, file)
            if os.path.exists(self.work_dir + '/blast_tmp/' + os.path.basename(_path)):
                os.remove(self.work_dir + '/blast_tmp/' + os.path.basename(_path))
            os.link(_path, self.work_dir + '/blast_tmp/' + os.path.basename(_path))
        self.catblast.set_options({"blastout": self.work_dir + '/blast_tmp'})
        self.catblast_tools.append(self.catblast)
        self.catblast.on('start', self.set_step, {'start': self.step.cat_blastout})
        if len(self.catblast_tools) == 1:
            self.catblast.on('end', self.set_step, {'end': self.step.cat_blastout})
            self.catblast.on('end', self.set_output)
        else:
            self.on_rely(self.catblast_tools, self.set_output)
        for i in self.catblast_tools:
            i.run()

    def set_output(self):
        self.set_step(event={'data': {'end': self.step.cat_blastout}})
        # for root, dirs, files in os.walk(self.output_dir):
        #     for names in files:
        #         os.remove(os.path.join(root, names))
        for file in os.listdir(self.output_dir):
            file_path = os.path.join(self.output_dir, file)
            os.remove(file_path)
        if self.option('outfmt') == 6:
            from mbio.files.align.blast.blast_xml import BlastXmlFile
            xml = BlastXmlFile()
            xml.set_path(self.catblast.option('blast_result').prop['path'])
            xml.convert2table(self.output_dir + "/blast.table")
            self.option('outtable').set_path(self.output_dir + "/blast.table")
        os.link(self.catblast.option('blast_result').prop['path'], self.output_dir + '/blast.xml')
        self.option('outxml', self.output_dir + '/blast.xml')
        self.option('outxml_dir',self.work_dir + '/blast_tmp')
        self.end()

    def run(self):
        super(MetaDiamondModule, self).run()
        self.run_splitfasta()

    def end(self):
        repaths = [
            [".", "", "blast输出目录"],
            ["blast.xml", "xml", "blast xml输出结果文件"],
            ["blast_table.xls", "xls", "blast xls输出结果文件"]
        ]
        sdir = self.add_upload_dir(self.output_dir)
        sdir.add_relpath_rules(repaths)
        super(MetaDiamondModule, self).end()

# -*- coding: utf-8 -*-
# __author__ = 'shijin'
# last_modify:2017.03.17

from biocluster.module import Module
import os
import shutil
import re
from biocluster.core.exceptions import OptionError
from biocluster.config import Config
from mbio.packages.rna.annot_config import AnnotConfig

class DiamondModule(Module):
    def __init__(self, work_id):
        super(DiamondModule, self).__init__(work_id)
        self._fasta_type = {'Protein': 'prot', 'DNA': 'nucl'}
        options = [
            {"name": "query", "type": "infile", "format": "labelfree.common"},  # 输入文件
            {"name": "lines", "type": "int", "default": 500},  # 将fasta序列拆分此行数的多个文件
            {"name": "query_type", "type": "string", "default": "prot"},  # 输入的查询序列的格式，为nucl或者prot
            {"name": "database", "type": "string", "default": "nr"},
            {"name": "kegg_species", "type": "string", "default": None},
            # 比对数据库 nt nr string swissprot kegg customer_mode
            {"name": "outfmt", "type": "int", "default": 5},  # 输出格式，只为5
            {"name": "blast", "type": "string"},
            {"name": "reference", "type": "infile", "format": "sequence.fasta"},  # 参考序列  选择customer时启用
            {"name": "reference_type", "type": "string", "default": "prot"},  # 参考序列(库)的类型  为nucl或者prot
            {"name": "evalue", "type": "float", "default": 1e-5},  # evalue值
            {"name": "identity", "type": "float", "default": 0},  # identity过滤条件
            {"name": "similarity", "type": "float", "default": 0},  # similarity过滤条件
            {"name": "num_threads", "type": "int", "default": 10},  # cpu数
            {"name": "outxml", "type": "outfile", "format": "align.blast.blast_xml"},  # 输出格式为5时输出
            {"name": "outtable", "type": "outfile", "format": "align.blast.blast_table"},  # 输出格式为6时输出
            {'name': 'kegg_version', 'type': 'string', 'default': "202007"},
            {"name": "nr_version", "type": "string", "default": "2019"},
            {"name": "swissprot_version", "type": "string", "default": "2019"},
            {"name": "string_version", "type": "string", "default": "2019"},
            {"name": "eggnog_version", "type": "string", "default": "2019"},
            {"name": "version", "type": "string", "default": "2019"},
            {"name": "diamond_version", "type": "string", "default": "v0.9.24.125"},
            # 当输出格式为非5，6时，只产生文件不作为outfile
        ]
        self.add_option(options)
        self.catblast = self.add_tool("labelfree.annotation.cat_blastout")
        self.splitfasta = self.add_tool("labelfree.annotation.split_fasta")
        self.xmlfilter = self.add_tool("labelfree.annotation.filter_annot")
        self.step.add_steps('blast', 'split_fasta', 'cat_blastout', 'blast_filter')
        self.blast_tools = []
        self.catblast_tools = []
        if self.option('kegg_version'):
            self.kegg_files_dict = AnnotConfig().get_file_dict(db="kegg", version=self.option("kegg_version"))

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def check_options(self):
        if not self.option("query").is_set:
            raise OptionError("必须设置参数query")
        if self.option('query_type') not in ['nucl', 'prot']:
            raise OptionError('query_type查询序列的类型为nucl(核酸)或者prot(蛋白):{}'.format(self.option('query_type')))
        if self.option('outfmt') not in [5, 6]:
            raise OptionError('outfmt遵循blast+输出规则，目前只支持5，6输出格式：{}'.format(self.option('outfmt')))
        if not 1 > self.option('evalue') >= 0:
            raise OptionError('E-value值设定必须为[0-1)之间：{}'.format(self.option('evalue')))
        if not isinstance(self.option('lines'), int):
            raise OptionError("行数必须为整数")
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
            "version": self.option("version"),
            "num_threads": self.option("num_threads"),
            "kegg_version": self.option("kegg_version"),
            "diamond_version": self.option("diamond_version")
        }
        for opt in ["kegg_version", "nr_version", "swissprot_version", "string_version", "eggnog_version"]:
            opts.update({opt: self.option(opt)})


        # os.path.join(Config().SOFTWARE_DIR, "database/KEGG/genome2.xls")
        if self.option("kegg_species"):
            try:
                ref_path = AnnotConfig().get_kegg_fa(version=self.option("kegg_version"), species=self.option("kegg_species"))
            except:
                self.set_error("无法找到物种 {} 的数据".format(self.option("kegg_species")))
            opts['reference'] = ref_path
            opts['database'] = 'customer_mode'

        for f in os.listdir(self.splitfasta.output_dir):
            opts['query'] = os.path.join(self.splitfasta.output_dir, f)
            blast_tool = self.add_tool('labelfree.annotation.diamond')
            blast_tool.set_options(opts)
            # blast_tool.run()
            self.blast_tools.append(blast_tool)
        if len(self.blast_tools) == 1:
            self.blast_tools[0].on("end", self.run_catblastout)
        else:
            self.on_rely(self.blast_tools, self.run_catblastout)
        for tool in self.blast_tools:  # edited by shijin on 20170627
            tool.run()

    def run_blast_filter(self):
        options = {
            'xml': self.catblast.option('blast_result').prop['path'],
            'types': "xml",
            'evalue': self.option('evalue'),
            'identity': self.option('identity'),
            'similarity': self.option('similarity')
        }
        self.xmlfilter.set_options(options)
        self.xmlfilter.on('start', self.set_step, {'start': self.step.blast_filter})
        self.xmlfilter.on('end', self.set_step, {'end': self.step.blast_filter})
        self.xmlfilter.on('end', self.set_output)
        self.xmlfilter.run()

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
            self.catblast.on('end', self.run_blast_filter)
        else:
            self.on_rely(self.catblast_tools, self.run_blast_filter)
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
            xml.set_path(self.xmlfilter.option('outxml').prop['path'])
            xml.convert2table(self.output_dir + "/blast.table")
            self.option('outtable').set_path(self.output_dir + "/blast.table")
        os.link(self.xmlfilter.option('outxml').prop['path'], self.output_dir + '/blast.xml')
        self.option('outxml', self.output_dir + '/blast.xml')
        self.end()

    def run(self):
        super(DiamondModule, self).run()
        self.run_splitfasta()

    def end(self):
        repaths = [
            [".", "", "blast输出目录"],
            ["blast.xml", "xml", "blast xml输出结果文件"],
            ["blast_table.xls", "xls", "blast xls输出结果文件"]
        ]
        sdir = self.add_upload_dir(self.output_dir)
        sdir.add_relpath_rules(repaths)
        super(DiamondModule, self).end()

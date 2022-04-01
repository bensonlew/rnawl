# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'
# last_modify:2016.11.16

from biocluster.module import Module
import os
import shutil
from biocluster.core.exceptions import OptionError
import unittest


class BlastModule(Module):
    def __init__(self, work_id):
        super(BlastModule, self).__init__(work_id)
        self._fasta_type = {'Protein': 'prot', 'DNA': 'nucl'}
        self._blast_type = {'nucl': {'nucl': ['blastn', 'tblastn'],
                                     'prot': ['blastx']},
                            'prot': {'nucl': [],
                                     'prot': ['blastp']}}
        self._database_type = {'nt': 'nucl', 'nr': 'prot', 'kegg': 'prot', 'swissprot': 'prot', 'string': 'prot'}
        options = [
            {"name": "query", "type": "infile", "format": "small_rna.fasta"},  # 输入文件
            {"name": "lines", "type": "int", "default": 50000},  # 将fasta序列拆分此行数的多个文件
            {"name": "query_type", "type": "string", "default": "nucl"},  # 输入的查询序列的格式，为nucl或者prot
            {"name": "outfmt", "type": "int", "default": 5},  # 输出格式，只为5，6
            {"name": "blast", "type": "string", "default": "blastn"},  # 设定blast程序有blastn，blastx，blastp，tblastn，此处需要严格警告使用者必须选择正确的比对程序
            {"name": "evalue", "type": "float", "default": 1e-5},  # evalue值
            {"name": "num_threads", "type": "int", "default": 10},  # cpu数
            {"name": "num_alignment", "type": "int", "default": 5},  # 序列比对最大输出条数，默认500
            {"name": "outxml", "type": "outfile", "format": "align.blast.blast_xml"},  # 输出格式为5时输出
            {"name": "outtable", "type": "outfile", "format": "align.blast.blast_table"},  # 输出格式为6时输出
            {"name": "reference", "type": "infile", "format": "small_rna.fasta"}, # 序列比对文件
            {"name": "database", "type": "string", "default":""}, # 比对的数据库类型
            {"name": "reference_type", "type": "string", "default":"nucl"}, # 比对的数据库文件格式
        ]
        self.add_option(options)
        self.catblast = self.add_tool("small_rna.srna.cat_blastout")
        self.splitfasta = self.add_tool("small_rna.srna.split_fasta")
        self.step.add_steps('blast_index', 'blast', 'split_fasta', 'cat_blastout')
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
            raise OptionError("必须设置参数query")
        if self.option('query_type') not in ['nucl', 'prot']:
            raise OptionError('query_type查询序列的类型为nucl(核酸)或者prot(蛋白):{}'.format(self.option('query_type')))
        if not isinstance(self.option('lines'), int):
            raise OptionError("行数必须为整数")
        if self.option('lines') <= 0:
            raise OptionError("行数小于等于0，请重设！")
        if self.option("database").lower() not in ["rfam", "repeat", "intron", "exon"]:
            raise OptionError("必须指定比对数据库")
        if self.option("database").lower() != "rfam":
            if not self.option("reference").is_set:
                raise OptionError("进行repeat比对必须提供repeat文件")
            new_reference = os.path.join(self.work_dir, os.path.basename(self.option("reference").prop['path']))
            if not os.path.exists(new_reference):
                os.link(self.option("reference").prop['path'], new_reference)
            self.option("reference").set_path(new_reference)
        return True

    def run_blast_index(self):
        self.blast_index = self.add_tool("small_rna.srna.blast_index")
        self.blast_index.set_options({
            "reference": self.option("reference"),
            "reference_type": "nucl",
        })
        self.blast_index.on('start', self.set_step, {'start': self.step.blast_index})
        self.blast_index.on('end', self.set_step, {'end': self.step.blast_index})
        self.blast_index.on('end', self.run_splitfasta)
        self.blast_index.run()

    def run_splitfasta(self):
        self.splitfasta.set_options({
            "fasta": self.option("query"),
            "lines": self.option("lines"),
        })
        self.splitfasta.on('start', self.set_step, {'start': self.step.split_fasta})
        self.splitfasta.on('end', self.set_step, {'end': self.step.split_fasta})
        self.splitfasta.on('end', self.run_blast)
        self.splitfasta.run()

    def run_blast(self):
        opts = {
            "query": self.option("query"),
            "query_type": self.option("query_type"),
            "outfmt": self.option("outfmt"),
            "blast": self.option("blast"),
            "evalue": self.option("evalue"),
            "num_threads": self.option("num_threads"),
            "num_alignment": self.option("num_alignment"),
            "database": self.option("database"),
            "reference": self.option("reference")
        }
        for f in os.listdir(self.splitfasta.output_dir):
            opts['query'] = os.path.join(self.splitfasta.output_dir, f)
            blast_tool = self.add_tool('small_rna.srna.blast')
            blast_tool.set_options(opts)
            self.blast_tools.append(blast_tool)
        if len(self.blast_tools) == 1:
            self.blast_tools[0].on("end", self.run_catblastout)
        else:
            self.on_rely(self.blast_tools, self.run_catblastout)
        for tool in self.blast_tools:
            tool.run()

    def run_catblastout(self):
        self.set_step(event={'data': {'end': self.step.blast}})
        if os.path.exists(self.work_dir + '/blast_tmp'):
            shutil.rmtree(self.work_dir + '/blast_tmp')
        os.mkdir(self.work_dir + '/blast_tmp')
        for i in self.blast_tools:
            if self.option('outfmt') == 6:
                _path = i.option('outtable').prop['path']
                if os.path.exists(self.work_dir + '/xml_tmp'):
                    shutil.rmtree(self.work_dir + '/xml_tmp')
                os.mkdir(self.work_dir + '/xml_tmp')
                xml_path = i.option('outxml').prop['path']
                os.link(xml_path, self.work_dir + '/xml_tmp/' + os.path.basename(xml_path))
            else:
                _path = i.option('outxml').prop['path']
            os.link(_path, self.work_dir + '/blast_tmp/' + os.path.basename(_path))
        self.catblast.set_options({"blastout": self.work_dir + '/blast_tmp'})
        # self.catblast.run()
        self.catblast_tools.append(self.catblast)
        if self.option('outfmt') == 6:
            self.tmp_tool = self.add_tool('align.ncbi.cat_blastout')
            self.tmp_tool.set_options({"blastout": self.work_dir + '/xml_tmp'})
            # self.tmp_tool.run()
            self.catblast_tools.append(self.tmp_tool)
        self.catblast.on('start', self.set_step, {'start': self.step.cat_blastout})
        if len(self.catblast_tools) == 1:
            self.catblast.on('end', self.set_step, {'end': self.step.cat_blastout})
            self.catblast.on('end', self.set_output)
        else:
            # self.on_rely(self.catblast_tools, self.set_step, {'end': self.step.cat_blastout})
            self.on_rely(self.catblast_tools, self.set_output)
        for i in self.catblast_tools:
            i.run()

    def set_output(self):
        self.set_step(event={'data': {'end': self.step.cat_blastout}})
        for root, dirs, files in os.walk(self.output_dir):
            for names in files:
                os.remove(os.path.join(root, names))
        if self.option('outfmt') == 6:
            os.link(self.catblast.option('blast_result').prop['path'], self.output_dir + '/blast_table.xls')
            self.option('outtable', self.output_dir + '/blast_table.xls')
            self.option('outxml', self.tmp_tool.option('blast_result').prop['path'])
        else:
            os.link(self.catblast.option('blast_result').prop['path'], self.output_dir + '/blast.xml')
            self.option('outxml', self.output_dir + '/blast.xml')
        self.end()

    def run(self):
        super(BlastModule, self).run()
        if self.option("database").lower() == "rfam":
            self.run_splitfasta()
        else:
            self.run_blast_index()

    def end(self):
        repaths = [
            [".", "", "blast输出目录"],
            ["blast.xml", "xml", "blast xml输出结果文件"],
            ["blast_table.xls", "xls", "blast xls输出结果文件"]
        ]
        sdir = self.add_upload_dir(self.output_dir)
        sdir.add_relpath_rules(repaths)
        super(BlastModule, self).end()

class TestFunction(unittest.TestCase):

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "Blast" + str(random.randint(1, 10000)),
            "type": "module",
            "name": "small_rna.srna.blast",
            "instant": False,
            "options": dict(
                query="/mnt/ilustre/users/sanger-test/workspace/20181011/Single_srna10-33-54/Srna/KnownMirna/filtered.fa",
                database="rfam",
                #reference="/mnt/ilustre/users/sanger-dev/workspace/20180503/Single_IntronExon_8468/IntronExon/output/ref.exon.fa",
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()
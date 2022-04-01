# -*- coding: utf-8 -*-
# __author__ = 'qiuping'
# last_modify:2016.11.16

from biocluster.module import Module
import os
import shutil
from biocluster.core.exceptions import OptionError


class BlastModule(Module):
    def __init__(self, work_id):
        super(BlastModule, self).__init__(work_id)
        self._fasta_type = {'Protein': 'prot', 'DNA': 'nucl'}
        self._blast_type = {'nucl': {'nucl': ['blastn', 'tblastn'],
                                     'prot': ['blastx']},
                            'prot': {'nucl': [],
                                     'prot': ['blastp']}}
        self._database_type = {'nt': 'nucl', 'nr': 'prot', 'kegg': 'prot', 'swissprot': 'prot', 'string': 'prot',
                               'nt_v20200604': 'nucl','swissprot_v20200617': 'prot','nt_v20210917': 'nucl',"nr_v20210917":"prot"}
        options = [
            {"name": "query", "type": "infile", "format": "sequence.fasta"},  # 输入文件
            {"name": "lines", "type": "int", "default": 500},  # 将fasta序列拆分此行数的多个文件
            {"name": "query_type", "type": "string"},  # 输入的查询序列的格式，为nucl或者prot
            {"name": "database", "type": "string", "default": "nr"},
            # 比对数据库 nt nr string swissprot kegg customer_mode
            {"name": "outfmt", "type": "int", "default": 5},  # 输出格式，只为5，6
            {"name": "blast", "type": "string"},  # 设定blast程序有blastn，blastx，blastp，tblastn，此处需要严格警告使用者必须选择正确的比对程序
            {"name": "reference", "type": "infile", "format": "sequence.fasta"},  # 参考序列  选择customer时启用
            {"name": "reference_type", "type": "string"},  # 参考序列(库)的类型  为nucl或者prot
            {"name": "evalue", "type": "float", "default": 1e-5},  # evalue值
            {"name": "num_threads", "type": "int", "default": 10},  # cpu数
            {"name": "memory", "type": "int", "default": 20},  # 内存G
            {"name": "num_alignment", "type": "int", "default": 5},  # 序列比对最大输出条数，默认500
            {"name": "outxml", "type": "outfile", "format": "align.blast.blast_xml"},  # 输出格式为5时输出
            {"name": "outtable", "type": "outfile", "format": "align.blast.blast_table"},  # 输出格式为6时输出
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
            raise OptionError("必须设置参数query", code="21100101")
        if self.option('query_type') not in ['nucl', 'prot']:
            raise OptionError('query_type查询序列的类型为nucl(核酸)或者prot(蛋白):%s', variables=(self.option('query_type')), code="21100102")
        else:
            if self._fasta_type[self.option('query').prop['seq_type']] != self.option('query_type'):
                raise OptionError(
                    '文件检查发现查询序列为:%s, 而需要的文件类型为:%s', variables=(
                        self._fasta_type[self.option('query').prop['seq_type'], self.option('query_type')]), code="21100103")
        if self.option("database") == 'customer_mode':
            if not self.option("reference").is_set:
                raise OptionError("使用自定义数据库模式时必须设置reference", code="21100104")
            if self.option('reference_type') not in ['nucl', 'prot']:
                raise OptionError('reference_type参考序列的类型为nucl(核酸)或者prot(蛋白):%s', variables=(self.option('query_type')), code="21100105")
            else:
                if self._fasta_type[self.option('reference').prop['seq_type']] != self.option('reference_type'):
                    raise OptionError(
                        '文件检查发现参考序列为:%s, 而需要的文件类型为:%s', variables=(
                            self._fasta_type[self.option('reference').prop['seq_type'], self.option('reference_type')]), code="21100106")
        elif self.option("database").lower() not in ["nt", "nr", "string", 'kegg', 'swissprot', "nt_v20200604","nt_v20210917",
                                                     "swissprot_v20200617"] and self.option("database") not in ['fgr/amoA_archaea_202012',
                                                    'fgr/amoA_bacteria_202012', 'nt', "nt_v20200604","nr_v20210917", 'fgr/amoA_AOB_like_202012',
                                                     'fgr/amoA_comammox_202012', 'fgr/nosZ_202012',
                                                     'fgr/nosZ_atypical_1_202012', 'fgr/nosZ_atypical_2_202012',
                                                     'fgr/nirK_202012', 'fgr/nirS_202012', 'fgr/mcrA_202012', 'fgr/nifH_202012',
                                                     'fgr/pmoA_202012', 'fgr/mmoX_202012']:
            raise OptionError("数据库%s不被支持", variables=(self.option("database")), code="21100107")
        else:
            if self.option("database").lower() in ["nt", "nr", "string", 'kegg', 'swissprot', "nt_v20200604","nt_v20210917","nr_v20210917",
                                                     "swissprot_v20200617"]:
                self.option('reference_type', self._database_type[self.option("database").lower()])
        if self.option('outfmt') not in [5, 6]:
            raise OptionError('outfmt遵循blast+输出规则，目前只支持5，6输出格式：%s', variables=(self.option('outfmt')), code="21100108")
        if not 1 > self.option('evalue') >= 0:
            raise OptionError('E-value值设定必须为[0-1)之间：%s', variables=(self.option('evalue')), code="21100109")
        if not 0 < self.option('num_alignment') < 1001:
            raise OptionError('序列比对保留数必须设置在1-1000之间:%s', variables=(self.option('num_alignment')), code="21100110")
        if self.option("database").lower() in ["nt", "nr", "string", 'kegg', 'swissprot', "nt_v20200604","nt_v20210917","nr_v20210917",
                                               "swissprot_v20200617"]:
            if self.option('blast') not in self._blast_type[self.option('query_type')][self.option('reference_type')]:
                raise OptionError(
                    '程序不试用于提供的查询序列和库的类型，请仔细检查，核酸比对核酸库只能使用blastn或者tblastn，\
                    核酸比对蛋白库只能使用blastp， 蛋白比对蛋白库只能使用blastp, 或者没有提供blast参数', code="21100111")
        if not isinstance(self.option('lines'), int):
            raise OptionError("行数必须为整数", code="21100112")
        # if self.option('lines') % 2:
        #     raise OptionError("行数必须为整除2")
        if self.option('lines') <= 0:
            raise OptionError("行数小于等于0，请重设！", code="21100113")
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
            "outfmt": self.option("outfmt"),
            "blast": self.option("blast"),
            "reference": self.option("reference"),
            "reference_type": self.option("reference_type"),
            "evalue": self.option("evalue"),
            "num_threads": self.option("num_threads"),
            "memory": self.option("memory"),
            "num_alignment": self.option("num_alignment"),
        }
        for f in os.listdir(self.splitfasta.output_dir):
            opts['query'] = os.path.join(self.splitfasta.output_dir, f)
            blast_tool = self.add_tool('align.blast')
            blast_tool.set_options(opts)
            # blast_tool.run()
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
                if os.path.exists(self.work_dir + '/xml_tmp/' + os.path.basename(xml_path)):
                    os.remove(self.work_dir + '/xml_tmp/' + os.path.basename(xml_path))
                os.link(xml_path, self.work_dir + '/xml_tmp/' + os.path.basename(xml_path))
            else:
                _path = i.option('outxml').prop['path']
            if os.path.exists(self.work_dir + '/blast_tmp/' + os.path.basename(_path)):
                os.remove(self.work_dir + '/blast_tmp/' + os.path.basename(_path))
            os.link(_path, self.work_dir + '/blast_tmp/' + os.path.basename(_path))
        self.catblast.set_options({"blastout": self.work_dir + '/blast_tmp'})
        # self.catblast.run()
        self.catblast_tools.append(self.catblast)
        if self.option('outfmt') == 6:
            self.tmp_tool = self.add_tool('align.cat_blastout')
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
            if os.path.exists(self.output_dir + '/blast_table.xls'):
                os.remove(self.output_dir + '/blast_table.xls')
            os.link(self.catblast.option('blast_result').prop['path'], self.output_dir + '/blast_table.xls')
            self.option('outtable', self.output_dir + '/blast_table.xls')
            self.option('outxml', self.tmp_tool.option('blast_result').prop['path'])
        else:
            if os.path.exists(self.output_dir + '/blast.xml'):
                os.remove(self.output_dir + '/blast.xml')
            os.link(self.catblast.option('blast_result').prop['path'], self.output_dir + '/blast.xml')
            self.option('outxml', self.output_dir + '/blast.xml')
        self.end()

    def run(self):
        super(BlastModule, self).run()
        self.run_splitfasta()

    def end(self):
        repaths = [
            [".", "", "blast输出目录"],
            ["blast.xml", "xml", "blast xml输出结果文件"],
            ["blast_table.xls", "xls", "blast xls输出结果文件"]
        ]
        sdir = self.add_upload_dir(self.output_dir)
        sdir.add_relpath_rules(repaths)
        super(BlastModule, self).end()

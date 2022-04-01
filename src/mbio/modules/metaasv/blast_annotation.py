# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

from biocluster.module import Module
import os,re
import shutil
import gevent
from biocluster.core.exceptions import OptionError
from mbio.packages.metaasv.common_function import link_dir,link_file


class BlastAnnotationModule(Module):
    """
    metaasv
    功能：1. 能对identity和coverage进行筛选
    2. 数据库比较齐全（metaasv库全部共有32种数据库）
    nt用老流程，其他数据库用新流程
    比对 + 注释
    """
    def __init__(self, work_id):
        super(BlastAnnotationModule, self).__init__(work_id)
        self.DATABASE = ['unite7.2/its_fungi','unite8.0/its_fungi',
                        'fgr/amoA', 'fgr/nosZ', 'fgr/nirK', 'fgr/nirS','fgr/nifH', 'fgr/pmoA', 'fgr/mmoX', 'fgr/mcrA', 'fgr/amoA_archaea', 'fgr/amoA_bacteria',
                        'maarjam081/AM', 'Protist_PR2_v4.5',
                        'silva132/16s_archaea', 'silva132/16s_bacteria','silva132/18s_eukaryota', 'silva132/16s',
                        'silva138/16s_archaea', 'silva138/16s_bacteria','silva138/18s_eukaryota', 'silva138/16s',
                        'greengenes135/16s', 'greengenes135/16s_archaea', 'greengenes135/16s_bacteria',
                         'nt_v20210917/16s_archaea', 'nt_v20210917/16s_bacteria', 'nt_v20210917/16s',
                         'nt_v20210917/18s_eukaryota', 'nt_v20210917/its_fungi', 'nt_v20210917',
                        'rdp11.5/16s', 'rdp11.5/16s_bacteria', 'rdp11.5/16s_archaea', 'nt', 'nt/16s','nt/18s','nt/its','nt_v20200327/16s_archaea', 'nt_v20200327/16s_bacteria','nt_v20200327/16s','Human_HOMD_v15.2', 'nt_v20200327/18s_eukaryota', 'nt_v20200327/its_fungi', "nt_v20200604",
                        'fgr/amoA_archaea_202012', 'fgr/amoA_bacteria_202012', 'fgr/amoA_AOB_like_202012', 'fgr/amoA_comammox_202012', 'fgr/nosZ_202012', 'fgr/nosZ_atypical_1_202012',
                        'fgr/nosZ_atypical_2_202012', 'fgr/nirK_202012', 'fgr/nirS_202012', 'fgr/mcrA_202012', 'fgr/nifH_202012', 'fgr/pmoA_202012', 'fgr/mmoX_202012']
        self._fasta_type = {'Protein': 'prot', 'DNA': 'nucl'}
        self._blast_type = {'nucl': {'nucl': ['blastn', 'tblastn'],
                                     'prot': ['blastx']},
                            'prot': {'nucl': [],
                                     'prot': ['blastp']}}
        self._database_type = {'nt': 'nucl', 'nr': 'prot', 'kegg': 'prot', 'swissprot': 'prot', 'string': 'prot'}
        options = [
            {"name": "query", "type": "infile", "format": "sequence.fasta"},  # 输入文件
            {"name": "lines", "type": "int", "default": 5000},  # 将fasta序列拆分此行数的多个文件
            {"name": "query_type", "type": "string", "default": "nucl"},  # 输入的查询序列的格式，为nucl或者prot
            {"name": "database", "type": "string", "default": "nr"},# 比对数据库 nt nr string swissprot kegg custom_mode
            {"name": "outfmt", "type": "int", "default": 5},  # 输出格式，只为5，6
            {"name": "blast", "type": "string","default": "blastn"},  # 设定blast程序有blastn，blastx，blastp，tblastn，此处需要严格警告使用者必须选择正确的比对程序
            {"name": "ref_fasta", "type": "infile", "format": "sequence.fasta"},  # 参考序列  选择customer时启用
            {'name': 'ref_taxon', 'type': 'infile', 'format': 'taxon.seq_taxon'},  # 参考taxon文件
            {"name": "reference_type", "type": "string","default": "nucl" },  # 参考序列(库)的类型  为nucl或者prot
            {"name": "evalue", "type": "float", "default": 1e-5},  # evalue值
            {"name": "num_threads", "type": "int", "default": 10},  # cpu数
            {"name": "memory", "type": "int", "default": 20},  # 内存G
            {"name": "num_alignment", "type": "int", "default": 1},  # 序列比对最大输出条数，默认1
            {'name': 'identity', 'type': 'int'},  # blast的比对identity
            {'name': 'coverage', 'type': 'int'},  # blast的比对coverage
            {"name": "meta_pipeline", "type": "string", "default": ""},  # 用于多样性选择氨基酸进行物种分类
        ]
        self.add_option(options)
        self.catblast = self.add_tool("sequence.cat_table")##m6格式
        self.catblast2 = self.add_tool("sequence.cat_table")##m6格式
        self.splitfasta = self.add_tool("sequence.split_fasta")
        self.step.add_steps('split_fasta', 'blast', 'cat_blastout')
        self.anno_dir = os.path.join(self.work_dir, "anno_dir")
        self.blast_dir = os.path.join(self.work_dir, "blast_dir")

        self.blast_tools = []
        self.catblast_tools = []

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def check_options(self):
        """
        参数二次检查
        :return:
        """
        if not self.option("query").is_set:
            raise OptionError("必须设置参数query")
        if self.option('query_type') not in ['nucl', 'prot']:
            raise OptionError('query_type查询序列的类型为nucl(核酸)或者prot(蛋白):%s', variables=(self.option('query_type')))
        else:
            if self._fasta_type[self.option('query').prop['seq_type']] != self.option('query_type'):
                raise OptionError(
                    '文件检查发现查询序列为:%s, 而需要的文件类型为:%s', variables=(
                        self._fasta_type[self.option('query').prop['seq_type'], self.option('query_type')]))
        if self.option("database") == 'custom_mode':
            if not self.option("ref_fasta").is_set:
                raise OptionError("使用自定义数据库模式时必须设置ref_fasta")
            if self.option('reference_type') not in ['nucl', 'prot']:
                raise OptionError('reference_type参考序列的类型为nucl(核酸)或者prot(蛋白):%s', variables=(self.option('query_type')))
            else:
                if self._fasta_type[self.option('ref_fasta').prop['seq_type']] != self.option('reference_type'):
                    raise OptionError(
                        '文件检查发现参考序列为:%s, 而需要的文件类型为:%s', variables=(
                            self._fasta_type[self.option('ref_fasta').prop['seq_type'], self.option('reference_type')]))
        elif self.option("database").lower() not in self.DATABASE and (self.option("database") not in self.DATABASE):
            raise OptionError("数据库%s不被支持", variables=(self.option("database")))

        if self.option('outfmt') not in [5, 6]:
            raise OptionError('outfmt遵循blast+输出规则，目前只支持5，6输出格式：%s', variables=(self.option('outfmt')))
        if not 1 > self.option('evalue') >= 0:
            raise OptionError('E-value值设定必须为[0-1)之间：%s', variables=(self.option('evalue')))
        if not 0 < self.option('num_alignment') < 1001:
            raise OptionError('序列比对保留数必须设置在1-1000之间:%s', variables=(self.option('num_alignment')))
        if self.option("database").lower() in ["nt", "nr", "string", 'kegg', 'swissprot', "nt_v20200604","nt_v20210917",
                                               "swissprot_v20200617"]:
            if self.option('blast') not in self._blast_type[self.option('query_type')][self.option('reference_type')]:
                raise OptionError(
                    '程序不试用于提供的查询序列和库的类型，请仔细检查，核酸比对核酸库只能使用blastn或者tblastn，\
                    核酸比对蛋白库只能使用blastp， 蛋白比对蛋白库只能使用blastp, 或者没有提供blast参数')
        if not isinstance(self.option('lines'), int):
            raise OptionError("行数必须为整数")
        if self.option('lines') <= 0:
            raise OptionError("行数小于等于0，请重设！")
        return True

    def run_splitfasta(self):
        """
        切分文件
        :return:
        """
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
        """
        blast比对（diamond）+ 注释
        :return:
        """
        opts = {
            "query_type": self.option("query_type"),
            "database": self.option("database"),
            "blast": self.option("blast"),
            "reference_type": self.option("reference_type"),
            "evalue": self.option("evalue"),
            "num_threads": self.option("num_threads"),
            "top_num": self.option("num_alignment"),
            "identity": self.option("identity"),
            "coverage": self.option("coverage"),
        }
        if self.option("database") in ['custom_mode']:
            opts["ref"] = self.option("ref_fasta")
            opts['ref_taxon'] = self.option("ref_taxon")
        for f in os.listdir(self.splitfasta.output_dir):
            opts['query'] = os.path.join(self.splitfasta.output_dir, f)
            blast_tool = self.add_tool('metaasv.blast')
            blast_tool.set_options(opts)
            self.blast_tools.append(blast_tool)
        if len(self.blast_tools) == 1:
            self.blast_tools[0].on("end", self.run_catblastout)
        else:
            self.on_rely(self.blast_tools, self.run_catblastout)
        for tool in self.blast_tools:
            tool.run()
            gevent.sleep(0)

    def run_catblastout(self):
        """
        blast的结果合并注释结果表
        m6格式的注释结果表，已经注释完毕
        :return:
        """
        if os.path.exists(self.anno_dir):
            shutil.rmtree(self.anno_dir)
        os.mkdir(self.anno_dir)
        i = 1
        for tool in self.blast_tools:
            for file in os.listdir(tool.output_dir):
                if re.search(r"ASV_tax_assignments", file):
                    file_path = os.path.join(tool.output_dir, file)
                    new_file_path = os.path.join(self.anno_dir, file + str(i))
                    link_file(file_path, new_file_path)
                    i += 1
        if os.path.exists(self.blast_dir):
            shutil.rmtree(self.blast_dir)
        os.mkdir(self.blast_dir)
        i = 1
        for tool in self.blast_tools:
            for file in os.listdir(tool.output_dir):
                if re.search(r"blast.m6.xls", file):
                    file_path = os.path.join(tool.output_dir, file)
                    new_file_path = os.path.join(self.blast_dir, file + str(i))
                    link_file(file_path, new_file_path)
                    i += 1

        options = {
            "fa_dir": self.anno_dir,
            "prefix": 'ASV_tax_assignments'
        }
        self.catblast.set_options(options)
        self.catblast_tools.append(self.catblast)
        options2 = {
            "fa_dir": self.blast_dir,
            "prefix": 'blast_table'
        }
        self.catblast2.set_options(options2)
        self.catblast_tools.append(self.catblast2)
        self.on_rely(self.catblast_tools, self.set_output)
        for tool in self.catblast_tools:
            tool.run()

    def run_nt_blast(self):
        """
        nt 数据库比对 增加了功能基因蛋白数据的比对 20210520
        :return:
        """
        self.nt_blast = self.add_module("align.blast")
        self.nt_blast.set_options({
            "query": self.option("query"),
            "lines": self.option("lines"),
            "query_type": self.option("query_type"),
            "database": self.option("database"),
            "blast": self.option("blast"),
            "reference_type": self.option("reference_type"),
            "evalue": self.option("evalue"),
            "num_threads": self.option("num_threads"),
            "memory": self.option("memory"),
            "num_alignment": self.option("num_alignment")
        })
        self.nt_blast.on("end", self.run_taxon)
        self.nt_blast.run()

    def run_taxon(self):
        """
        nt注释比对结果
        :return:
        """
        self.ncbi_taxon = self.add_tool("annotation.mg_ncbi_taxon")
        self.ncbi_taxon.set_options({
            "blastout": self.nt_blast.option("outxml"),
            #"out_type": 1,
            "blastdb": self.option("database")
        })
        self.ncbi_taxon.on("end", self.set_output)
        self.ncbi_taxon.run()

    def run_remove_header(self,infile, outfile):
        """
        对合并后的blast结果做一步处理
        :return:
        """
        with open(infile, 'r') as f, open(outfile, 'w') as w:
            w.write("ASV_ID\tTaxonomy\tIdentity\n")
            for line in f:
                w.write(line)

    def set_output(self):
        """
        设置结果文件目录
        :return:
        """
        for root, dirs, files in os.walk(self.output_dir):
            for names in files:
                os.remove(os.path.join(root, names))
        if self.option('database') in ['nt', "nt_v20200604","nt_v20210917"]:
            link_file(self.ncbi_taxon.output_dir+"/query_taxons.xls", self.output_dir + '/ASV_tax_assignments.txt')
            link_file(self.ncbi_taxon.output_dir+"/blast_table.xls", self.output_dir + '/blast_table.xls')
        else:
            self.run_remove_header(self.catblast2.output_dir+"/blast_table.xls", self.work_dir+"/blast_table.xls")
            link_file(os.path.join(self.catblast.output_dir, "ASV_tax_assignments.xls"), self.output_dir + '/ASV_tax_assignments.txt')
        self.end()

    def run(self):
        """
        运行
        :return:
        """
        super(BlastAnnotationModule, self).run()
        if self.option("database") in ['nt', "nt_v20200604","nt_v20210917"]:
            self.run_nt_blast()
        elif self.option("meta_pipeline") == "functional_gene_prot":
            self.run_blast()
        else:
            self.run_splitfasta()

    def end(self):
        """
        结束运行
        :return:
        """
        super(BlastAnnotationModule, self).end()

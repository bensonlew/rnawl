# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'
# last_modify:2017.8.23

from biocluster.module import Module
import os
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.config import Config
from mbio.packages.annotation.mg_annotation.mg_taxon import mg_taxon


class MetaTaxModule(Module):
    def __init__(self, work_id):
        super(MetaTaxModule, self).__init__(work_id)
        options = [
            {'name': 'fasta', 'type': 'infile', 'format': 'sequence.fasta'},  # 输入fasta文件
            {'name': 'revcomp', 'type': 'bool', 'default': False},  # 序列是否翻转只针对用qiime比对的序列
            {'name': 'confidence', 'type': 'float', 'default': 0.7},  # 置信度值只针对用qiime比对的序列
            # {"name": "customer_mode", "type": "bool", "default": False},  # customer 自定义数据库
            {'name': 'database', 'type': 'string'},  # 数据库选择
            {'name': 'ref_fasta', 'type': 'infile', 'format': 'sequence.fasta'},  # 参考fasta序列
            {'name': 'ref_taxon', 'type': 'infile', 'format': 'taxon.seq_taxon'},  # 参考taxon文件
            {'name': 'taxon_file', 'type': 'outfile', 'format': 'taxon.seq_taxon'},  # 输出序列的分类信息文件
            {"name": "lines", "type": "int", "default": 500},  # 将fasta序列拆分此行数的多个文件
            {"name": "query_type", "type": "string"},  # 输入的查询序列的格式，为nucl或者prot
            {"name": "blast", "type": "string"},  # 设定blast程序有blastn，blastx，blastp，tblastn，此处需要严格警告使用者必须选择正确的比对程序
            {"name": "reference_type", "type": "string"},  # blast参考序列(库)的类型  为nucl或者prot
            {"name": "evalue", "type": "float", "default": 1e-5},  # blast比对evalue值
            {"name": "num_threads", "type": "int", "default": 10},  # blast比对cpu数
            {"name": "memory", "type": "int", "default": 30},  # 内存G
            {"name": "num_alignment", "type": "int", "default": 1},  # 序列比对最大输出条数，默认500
            {"name": "meta_pipeline", "type": "string", "default": ""}, # 用于多样性选择氨基酸进行物种分类
        ]
        self._fasta_type = {'Protein': 'prot', 'DNA': 'nucl'}
        self._blast_type = {'nucl': {'nucl': ['blastn', 'tblastn'],
                                     'prot': ['blastx']},
                            'prot': {'nucl': [],
                                     'prot': ['blastp']}}
        self._database_type = {'nt': 'nucl', 'nr': 'prot', 'kegg': 'prot', 'swissprot': 'prot', 'string': 'prot'}
        self.add_option(options)  ##### 检查option是否list格式，其中每个opt是否字典格式
        self.blast = self.add_module("align.blast")
        self.diamond2 = self.add_module("metaasv.diamond2")
        self.ncbi_taxon = self.add_tool("taxon.ncbi_taxon")
        self.qiime = self.add_tool("taxon.qiime_assign")

    def check_options(self):
        if not self.option("fasta").is_set:
            raise OptionError("必须设置参数fasta", code="21201001")
        if self.option("revcomp") not in [True, False]:
            raise OptionError("必须设置序列是否翻转", code="21201002")
        if self.option('database') == "custom_mode":
            if not self.option("ref_fasta").is_set or not self.option("ref_taxon").is_set:
                raise OptionError("数据库自定义模式必须设置参考fasta序列和参考taxon文件", code="21201003")
        elif self.option('database') == "customer_mode":
            if not self.option("ref_fasta").is_set or not self.option("ref_taxon").is_set:
                raise OptionError("数据库自定义模式必须设置参考fasta序列和参考taxon文件", code="21201003")
        else:
            if self.option("database") not in ['silva123/16s_bacteria', 'silva123/16s_archaea',
                                               'silva123/16s', 'silva123/18s_eukaryota', 'silva123',
                                               'silva119/16s_bacteria', 'silva119/16s_archaea',
                                               'silva119/16s', 'silva119/18s_eukaryota',
                                               'unite8.0/its_fungi', 'unite7.2/its_fungi', 'unite7.0/its_fungi',
                                               'fgr/amoA', 'fgr/nosZ', 'fgr/nirK', 'fgr/nirS',
                                               'fgr/nifH', 'fgr/pmoA', 'fgr/mmoX', 'fgr/mcrA',
                                               'fgr/amoA_archaea', 'fgr/amoA_bacteria',
                                               'maarjam081/AM', 'Human_HOMD', 'Human_HOMD_v15.2', 'Human_HPB', 'Protist_PR2_v4.5',
                                               'silva128/16s_archaea', 'silva128/16s_bacteria',
                                               'silva128/18s_eukaryota', 'silva128/16s',
                                               'silva132/16s_archaea', 'silva132/16s_bacteria',
                                               'silva132/18s_eukaryota', 'silva132/16s',
                                               'silva138/16s_archaea', 'silva138/16s_bacteria',
                                               'silva138/18s_eukaryota', 'silva138/16s',
                                               'greengenes135/16s', 'greengenes135/16s_archaea',
                                               'greengenes135/16s_bacteria',
                                               'rdp11.5/16s', 'rdp11.5/16s_bacteria', 'rdp11.5/16s_archaea',
                                               'nt_v20200327/16s_archaea', 'nt_v20200327/16s_bacteria','nt_v20200327/16s',"nr_v20210917",
                                               'nt_v20200327/18s_eukaryota', 'nt_v20200327/its_fungi','fgr/amoA_archaea_202012','fgr/amoA_bacteria_202012',
                                               'nt', "nt_v20200604",'fgr/amoA_AOB_like_202012', 'fgr/amoA_comammox_202012', 'fgr/nosZ_202012',
                                               'fgr/nosZ_atypical_1_202012', 'fgr/nosZ_atypical_2_202012', 'fgr/nirK_202012',
                                               'fgr/nirS_202012', 'fgr/mcrA_202012', 'fgr/nifH_202012', 'fgr/pmoA_202012','fgr/mmoX_202012',
                                               'nt_v20210917/16s_archaea','nt_v20210917/16s_bacteria', 'nt_v20210917/16s',
                                               'nt_v20210917/18s_eukaryota', 'nt_v20210917/its_fungi', "nt_v20210917"]:
                # 多样性目前的database
                raise OptionError("数据库%s不被支持", variables=(self.option("database")), code="21201004")
        if self.option("database") in ['nt', "nt_v20200604","nt_v20210917"]:
            if self.option('query_type') not in ['nucl', 'prot']:
                raise OptionError('query_type查询序列的类型为nucl(核酸)或者prot(蛋白):%s', variables=(self.option('query_type')), code="21201005")
            else:
                if self._fasta_type[self.option('fasta').prop['seq_type']] != self.option('query_type'):
                    raise OptionError(
                        '文件检查发现查询序列为:%s, 而需要的文件类型为:%s', variables=(
                            self._fasta_type[self.option('fasta').prop['seq_type'], self.option('query_type')]), code="21201006")
            self.option('reference_type', 'nucl')
            if not 1 > self.option('evalue') >= 0:
                raise OptionError('E-value值设定必须为[0-1)之间：%s', variables=(self.option('evalue')), code="21201007")
            if not 0 < self.option('num_alignment') < 1001:
                raise OptionError('序列比对保留数必须设置在1-1000之间:%s', variables=(self.option('num_alignment')), code="21201008")
            if self.option('blast') not in self._blast_type[self.option('query_type')][self.option('reference_type')]:
                raise OptionError(
                    '程序不试用于提供的查询序列和库的类型，请仔细检查，核酸比对核酸库只能使用blastn或者tblastn，\
                     核酸比对蛋白库只能使用blastp， 蛋白比对蛋白库只能使用blastp, 或者没有提供blast参数', code="21201009")
            if not isinstance(self.option('lines'), int):
                raise OptionError("行数必须为整数", code="21201010")
            # if self.option('lines') % 2:
            #     raise OptionError("行数必须为整除2")
            if self.option('lines') <= 0:
                raise OptionError("行数小于等于0，请重设！", code="21201011")
        elif self.option("database") in ["nr_v20210917"]:
            pass
        return True

    def run_qiime(self):
        self.qiime.set_options({
            "fasta": self.option("fasta"),
            "revcomp": self.option("revcomp"),
            "confidence": self.option("confidence"),
            "database": self.option("database"),
            "ref_fasta": self.option("ref_fasta"),
            "ref_taxon": self.option("ref_taxon")
        })
        self.qiime.on('end', self.set_output, 'qiime')
        self.qiime.run()

    def run_blast(self):
        self.blast.set_options({
            "query": self.option("fasta"),
            "lines": self.option("lines"),
            "query_type": self.option("query_type"),
            "database": self.option("database"),
            "blast": self.option("blast"),
            "reference": self.option("ref_fasta"),
            "reference_type": self.option("reference_type"),
            "evalue": self.option("evalue"),
            "num_threads": self.option("num_threads"),
            "memory": self.option("memory"),
            "num_alignment": self.option("num_alignment")
        })
        self.blast.on("end", self.run_taxon)
        self.blast.run()

    def run_diamond2(self):
        self.diamond2.set_options({
            "query": self.option("fasta"),
            "query_type": self.option("query_type"),
            "database": self.option("database"),
            "blast": self.option("blast"),
            "evalue": self.option("evalue")
        })
        self.diamond2.on("end", self.run_taxon)
        self.diamond2.run()

    def run_taxon(self):
        if self.option("database").split("/")[0] == "fgr":
            db_path = Config().SOFTWARE_DIR + "/database/Framebot/tax/" + self.option("database").split("/")[1] + ".tax"
            from mbio.packages.align.blast.xml2table import xml2table
            self.logger.info(self.blast.option("outxml"))
            table = xml2table(self.blast.option("outxml").prop["path"], self.work_dir + '/temp_blastable.xls')
            with open(table,"r") as f, open(db_path,"r") as v, open(self.output_dir + "/query_taxons.xls","w") as t:
                data1 = f.readlines()
                data2 = v.readlines()
                tax_dict = {}
                for x in data2:
                    tax_dict[x.strip().split("\t")[0]] = x.strip().split("\t")[1]
                for i in data1[1:]:
                    if i.strip().split()[10] in tax_dict:
                        t.write(i.strip().split()[5] + "\t" + tax_dict[i.strip().split()[10]] + "\n")

            if os.path.exists(self.output_dir + "/blast_table.xls"):
                os.remove(self.output_dir + "/blast_table.xls")
            os.link(self.work_dir + '/temp_blastable.xls', self.output_dir + "/blast_table.xls")
            self.option("taxon_file", self.output_dir + "/query_taxons.xls")
            self.end()
        elif self.option("meta_pipeline") == "functional_gene_prot" and self.option("database") == "customer_mode":
            db_path = self.option("ref_taxon").prop["path"]
            from mbio.packages.align.blast.xml2table import xml2table
            self.logger.info(self.blast.option("outxml"))
            table = xml2table(self.blast.option("outxml").prop["path"], self.work_dir + '/temp_blastable.xls')
            with open(table, "r") as f, open(db_path, "r") as v, open(self.output_dir + "/query_taxons.xls", "w") as t:
                data1 = f.readlines()
                data2 = v.readlines()
                tax_dict = {}
                for x in data2:
                    tax_dict[x.strip().split("\t")[0]] = x.strip().split("\t")[1]
                for i in data1[1:]:
                    if i.strip().split()[10] in tax_dict:
                        t.write(i.strip().split()[5] + "\t" + tax_dict[i.strip().split()[10]] + "\n")

            if os.path.exists(self.output_dir + "/blast_table.xls"):
                os.remove(self.output_dir + "/blast_table.xls")
            os.link(self.work_dir + '/temp_blastable.xls', self.output_dir + "/blast_table.xls")
            self.option("taxon_file", self.output_dir + "/query_taxons.xls")
            self.end()
        elif self.option("database") in ["nr_v20210917"]:
            from mbio.packages.align.blast.xml2table import xml2table
            table = xml2table(self.diamond2.output_dir+"/blast.xml", self.work_dir + '/temp_blastable.xls')
            query = self.filter_query(table)
            from mbio.packages.taxon.nr2taxon_v20210917 import accession_taxon
            id_taxons = accession_taxon(set(query.values()))
            for i in query:
                query[i] = (query[i], id_taxons[query[i]])
            with open(self.output_dir + '/query_taxons_detail.xls', 'w') as w:
                for item in query.iteritems():
                    if item[1][1]:
                        w.write(item[0] + '\t' + item[1][0] + '\t' + item[1][1] + '\n')

            tax = mg_taxon()
            tax.detail_to_level(self.output_dir + '/query_taxons_detail.xls', self.output_dir)
            if os.path.exists(self.output_dir + "/blast_table.xls"):
                os.remove(self.output_dir + "/blast_table.xls")
            os.link(self.work_dir + '/temp_blastable.xls', self.output_dir + "/blast_table.xls")
            self.option("taxon_file", self.output_dir + "/query_taxons.xls")
            self.end()
        else:
            self.ncbi_taxon.set_options({
                "blastout": self.blast.option("outxml"),
                "out_type": 1,
                "blastdb": self.option("database")
            })
            self.ncbi_taxon.on("end", self.set_output, 'ncbi')
            self.ncbi_taxon.run()

    def filter_query(self, fp):
        """生成query对gi的"""
        query = dict()
        openfp = open(fp)
        openfp.readline()
        for i in openfp:
            line_sp = i.strip().split('\t')
            if line_sp[5] in query:
                continue
            else:
                query[line_sp[5]] = line_sp[10]
        return query

    def set_output(self, event):
        obj = event["bind_object"]
        result = os.path.join(self.output_dir, 'seqs_tax_assignments.txt')
        if event['data'] == "qiime":
            self.option("taxon_file", obj.option("taxon_file"))
        else:
            self.option("taxon_file", obj.option("taxon_out"))
        if os.path.isfile(result):
            os.remove(result)  # 防止工作流重运行时报错 by ghd @ 20181012
        os.link(self.option("taxon_file").prop['path'],result)
        if self.option("database") in ["nt", 'nt_v20200604']:
            blast_file = self.ncbi_taxon.work_dir + '/temp_blastable.xls'
            blast_result = os.path.join(self.output_dir, 'blast_table.xls')
            if os.path.exists(blast_result):
                os.remove(blast_result)
            os.link(blast_file,blast_result)
        self.end()

    def run(self):
        super(MetaTaxModule, self).run()
        self.logger.info(self.option("database"))
        if self.option("database") in ["nr_v20210917"]:
            self.run_diamond2()
        elif self.option("database") in ["nt", "nt_v20200604","nt_v20210917"] or self.option("meta_pipeline") == "functional_gene_prot":
            self.run_blast()
        else:
            self.run_qiime()

    def end(self):
        super(MetaTaxModule, self).end()

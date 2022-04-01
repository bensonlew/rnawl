# -*- coding: utf-8 -*-
# __author__ = 'liubinxu1'
# last_modify:2018.05.08

from biocluster.module import Module
import os
import shutil
import re
from biocluster.core.exceptions import OptionError
from biocluster.config import Config
import unittest
from mbio.packages.rna.annot_config import AnnotConfig

class AnnotMapdbModule(Module):
    '''
    将序列与数据库比对 包括blast\pfam数据库和GO、COG、KEGG对应关系映射
    '''
    def __init__(self, work_id):
        super(AnnotMapdbModule, self).__init__(work_id)
        self._fasta_type = {'Protein': 'prot', 'DNA': 'nucl'}
        options = [
            {"name": "query", "type": "infile", "format": "ref_rna_v2.fasta"},  # 输入文件
            {"name": "querypep", "type": "infile", "format": "ref_rna_v2.fasta"},  # 输入文件
            {"name": "lines", "type": "int", "default": 10000},  # 将fasta序列拆分此行数的多个文件
            {"name": "query_type", "type": "string", "default": "nucl"},  # 输入的查询序列的格式，为nucl或者prot
            # {""},
            {"name": "database", "type": "string", "default": "nr;kegg;swissprot;eggnog"},

            {"name": "nr_db", "type": "string", "default": "nr"},
            {"name": "known_go", "type": "string", "default": None},
            {"name": "tax", "type": "bool", "default": False}, #要不要做物种分类
            {"name": "string_db", "type": "string", "default": "string"},
            {"name": "swissprot_db", "type": "string", "default": "swissprot"},
            {"name": "eggnog_db", "type": "string", "default": "eggnog"},
            {"name": "kegg_db", "type": "string", "default": "kegg"},
            {"name": "kegg_version", "type": "string", "default": ""},

            {"name": "kegg_species", "type": "string", "default": ""}, # 指定kegg数据库具体物种
            # 比对数据库 nt nr string swissprot kegg customer_mode
            {"name": "outfmt", "type": "int", "default": 5},  # 输出格式，只为5
            {"name": "blast", "type": "string"},
            {"name": "method", "type": "string", "default": "diamond"}, # 比对方法 diamond|blast
            {"name": "reference", "type": "infile", "format": "sequence.fasta"},  # 参考序列  选择customer时启用
            {"name": "reference_type", "type": "string", "default": "prot"},  # 参考序列(库)的类型  为nucl或者prot
            {"name": "evalue", "type": "float", "default": 1e-5},  # evalue值,blast参数
            # {"name": "identity", "type": "float", "default": 0},   # identity过滤条件
            # {"name": "similarity", "type": "float", "default": 0},  # similarity过滤条件
            {"name": "num_threads", "type": "int", "default": 10},  # blast cpu数
            {"name": "outxml", "type": "outfile", "format": "align.blast.blast_xml"},  # 输出格式为5时输出
            {"name": "outtable", "type": "outfile", "format": "align.blast.blast_table"},  # 输出格式为6时输出
            {"name": "blast2go_annot", "type": "infile", "format": "ref_rna_v2.blast2go_annot"},
            {"name": "merge_type", "type": "string", "default": "partial"},
            {"name": "go_version", "type": "string", "default": ""},
            {"name": "pfam_version", "type": "string", "default": "32"},
            {'name': 'kegg_version', 'type': 'string', 'default': "2019"},
            {"name": "nr_version", "type": "string", "default": "2019"},
            {"name": "swissprot_version", "type": "string", "default": "2019"},
            {"name": "string_version", "type": "string", "default": "2019"},
            {"name": "pir_version", "type": "string", "default": "2019"},
            {"name": "eggnog_version", "type": "string", "default": "2019"},
            {"name": "version", "type": "string", "default": "2019"},
            {"name": "diamond_version", "type": "string", "default": "v0.9.24.125"},

            # 当输出格式为非5，6时，只产生文件不作为outfile
        ]
        self.add_option(options)
        self.catblast = self.add_tool("ref_rna_v2.annotation.cat_blastout")
        self.splitfasta = self.add_tool("ref_rna_v2.annotation.split_fasta")
        self.xmlfilter = self.add_tool("ref_rna_v2.annotation.filter_annot")
        self.step.add_steps('blast2go', 'split_fasta',
                            'blast_nr', 'blast_string', 'blast_kegg', 'blast_swissprot', 'blast_eggnog',
                            'merge_nr', 'merge_string', 'merge_kegg', 'merge_swissprot', 'merge_eggnog')
        self.merge_step = {
            "nr": self.step.merge_nr,
            "string": self.step.merge_string,
            "kegg": self.step.merge_kegg,
            "swissprot": self.step.merge_swissprot,
            "eggnog": self.step.merge_eggnog,
        }
        # blast/diamond 分布运行结果
        self.blast_nr_tools = []
        self.blast_swissprot_tools = []
        self.blast_kegg_tools = []
        self.blast_string_tools = []
        self.blast_eggnog_tools = []

        self.hmm_pfam_tools = []

        self.catblast_tools = []
        self.string2cog_tools = []
        self.kegg2ko_tools = []
        self.nr2go_tools = []
        self.nr2ncbitax_tools = []
        self.merge_tools = {}

        self.blast_opts = {
            "query_type": self.option("query_type"),
            "outfmt": 5,
            "blast": self.option("blast"),
            "evalue": self.option("evalue"),
            "num_threads": self.option("num_threads"),
        }

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def check_options(self):
        if not self.option("query").is_set:
            raise OptionError("必须设置参数query", code = "23700301")
        if self.option('query_type') not in ['nucl', 'prot']:
            raise OptionError('query_type查询序列的类型为nucl(核酸)或者prot(蛋白):%s', variables = (self.option('query_type')), code = "23700302")
        if self.option('outfmt') not in [5, 6]:
            raise OptionError('outfmt遵循blast+输出规则，目前只支持5，6输出格式：%s', variables = (self.option('outfmt')), code = "23700303")
        if not 1 > self.option('evalue') >= 0:
            raise OptionError('E-value值设定必须为[0-1)之间：%s', variables = (self.option('evalue')), code = "23700304")
        if not isinstance(self.option('lines'), int):
            raise OptionError("行数必须为整数", code = "23700305")
        return True

    def run_splitfasta(self):
        '''
        fasta文件分割
        '''
        self.splitfasta.set_options({
            "fasta": self.option("query"),
            "lines": self.option("lines"),
        })
        self.splitfasta.on('start', self.set_step, {'start': self.step.split_fasta})
        self.splitfasta.on('end', self.set_step, {'end': self.step.split_fasta})
        self.splitfasta.on('end', self.set_step, {'start': self.step.blast_nr})
        self.splitfasta.on('end', self.set_step, {'start': self.step.blast_kegg})
        self.splitfasta.on('end', self.set_step, {'start': self.step.blast_string})
        self.splitfasta.on('end', self.set_step, {'start': self.step.blast_eggnog})
        self.splitfasta.on('end', self.set_step, {'start': self.step.blast_swissprot})
        self.splitfasta.on('end', self.run_blast)
        self.splitfasta.run()

    def run_blast_nr(self):
        '''
        设置NR GO NCBItaxon 运行逻辑
        '''
        opts = self.blast_opts.copy()
        opts.update({"database": self.option("nr_db")})
        opts.update({"version": self.option("version")})
        opts.update({"nr_version": self.option("nr_version")})
        opts.update({"diamond_version": self.option("diamond_version")})
        if self.option("diamond_version") >= "v2.0.13":
            opts.update({"sensitive": 1})
        i = 0
        for f in os.listdir(self.splitfasta.output_dir):
            opts['query'] = os.path.join(self.splitfasta.output_dir, f)
            if self.option('method') == "blast":
                blast_tool = self.add_tool('ref_rna_v2.annotation.blast')
            else:
                blast_tool = self.add_tool('ref_rna_v2.annotation.diamond')
            blast_tool.on("end", self.run_nr2go, i)
            if self.option("tax"):
                blast_tool.on("end", self.run_nr2ncbitax, i)
            blast_tool.set_options(opts)
            # blast_tool.run()
            self.blast_nr_tools.append(blast_tool)
            nr2go = self.add_tool('ref_rna_v2.annotation.nr2go')
            nr2ncbitaxon = self.add_tool('ref_rna_v2.annotation.ncbi_taxon')
            self.nr2go_tools.append(nr2go)
            self.nr2ncbitax_tools.append(nr2go)
            i = i + 1
        if len(self.blast_nr_tools) == 1:
            self.blast_nr_tools[0].on("end", self.run_catblastout, "nr")
            self.nr2go_tools[0].on("end", self.run_catgo)
            if self.option("tax"):
                self.nr2ncbitax_tools[0]("end", self.run_catncbitax)
        else:
            self.on_rely(self.blast_nr_tools, self.run_catblastout, "nr")
            self.on_rely(self.nr2go_tools, self.run_catgo)
            if self.option("tax"):
                self.on_rely(self.nr2ncbitax_tools, self.run_catncbitax)
        for tool in self.blast_nr_tools:
            tool.run()
    def run_catncbitax(self):
        '''
        暂时未实现
        '''

    def run_nr2go(self, event):
        '''
        获取GO注释
        '''
        obj = event['bind_object']
        options = {
            'blastout': obj.option("outxml"),
        }
        options.update({"pir_version": self.option("pir_version")})

        # 删除，在module中合并GO注释
        # if self.option("known_go"):
        #     options.update({
        #         'known_go': self.option("go_annot")
        #     })
        i = event['data']
        self.nr2go_tools[i].set_options(options)
        self.nr2go_tools[i].run()

    def run_nr2ncbitax(self, event):
        '''
        获取NCBI物种分类信息
        '''
        obj = event['bind_object']
        options = {
            'blastout': obj.output_dir + "/",
        }
        i = event['data']
        self.nr2ncbitax_tools[i].set_options(options)
        self.nr2ncbitax_tools[i].run()

    def run_blast_swissprot(self):
        '''
        swissprot blast比对
        '''
        opts = self.blast_opts.copy()
        opts.update({"database": self.option("swissprot_db")})
        opts.update({"version": self.option("version")})
        opts.update({"swissprot_version": self.option("swissprot_version")})
        opts.update({"diamond_version": self.option("diamond_version")})
        if self.option("diamond_version") >= "v2.0.13":
            opts.update({"sensitive": 1})
        i = 0
        for f in os.listdir(self.splitfasta.output_dir):
            opts['query'] = os.path.join(self.splitfasta.output_dir, f)
            if self.option('method') == "blast":
                blast_tool = self.add_tool('ref_rna_v2.annotation.blast')
            else:
                blast_tool = self.add_tool('ref_rna_v2.annotation.diamond')
            blast_tool.set_options(opts)
            self.blast_swissprot_tools.append(blast_tool)

        if len(self.blast_swissprot_tools) == 1:
            self.blast_swissprot_tools[0].on("end", self.run_catblastout, "swissprot")
        else:
            self.on_rely(self.blast_swissprot_tools, self.run_catblastout, "swissprot")
        for tool in self.blast_swissprot_tools:
            tool.run()

    def run_blast_string(self):
        '''
        string数据库比对 与 COG流程
        '''
        opts = self.blast_opts.copy()
        opts.update({"database": self.option("string_db")})
        opts.update({"version": self.option("version")})
        opts.update({"string_version": self.option("string_version")})
        opts.update({"diamond_version": self.option("diamond_version")})
        if self.option("diamond_version") >= "v2.0.13":
            opts.update({"sensitive": 1})
        i = 0
        for f in os.listdir(self.splitfasta.output_dir):
            opts['query'] = os.path.join(self.splitfasta.output_dir, f)
            if self.option('method') == "blast":
                blast_tool = self.add_tool('ref_rna_v2.annotation.blast')
            else:
                blast_tool = self.add_tool('ref_rna_v2.annotation.diamond')
            blast_tool.on("end", self.run_string2cog, i)
            blast_tool.set_options(opts)
            self.blast_string_tools.append(blast_tool)
            string2cog_tool = self.add_tool('ref_rna_v2.annotation.string2cog')
            self.string2cog_tools.append(string2cog_tool)
            i = i + 1
        if len(self.blast_string_tools) == 1:
            self.blast_string_tools[0].on("end", self.run_catblastout, "string")
            self.string2cog_tools[0].on("end", self.run_catcog)
        else:
            self.on_rely(self.blast_string_tools, self.run_catblastout, "string")
            self.on_rely(self.string2cog_tools, self.run_catcog)
        for tool in self.blast_string_tools:
            tool.run()

    def run_blast_eggnog(self):
        '''
        egg nog 数据库比对
        '''
        opts = self.blast_opts.copy()
        opts.update({"database": self.option("eggnog_db")})
        opts.update({"version": self.option("version")})
        opts.update({"eggnog_version": self.option("eggnog_version")})
        opts.update({"diamond_version": self.option("diamond_version")})
        if self.option("diamond_version") >= "v2.0.13":
            opts.update({"sensitive": 1})
        for f in os.listdir(self.splitfasta.output_dir):
            opts['query'] = os.path.join(self.splitfasta.output_dir, f)
            if self.option('method') == "blast":
                blast_tool = self.add_tool('ref_rna_v2.annotation.blast')
            else:
                blast_tool = self.add_tool('ref_rna_v2.annotation.diamond')
            blast_tool.set_options(opts)
            # blast_tool.run()
            self.blast_eggnog_tools.append(blast_tool)
        if len(self.blast_eggnog_tools) == 1:
            self.blast_eggnog_tools[0].on("end", self.run_catblastout, "eggnog")
        else:
            self.on_rely(self.blast_eggnog_tools, self.run_catblastout, "eggnog")
        for tool in self.blast_eggnog_tools:
            tool.run()

    def run_string2cog(self, event):
        '''
        string xml 转化为cog的脚本，暂时不在该处实现
        '''
        obj = event['bind_object']
        i = event['data']
        options = {
            'blastout': obj.output_dir + "/",
            'string_table': self.option('blast_string_table')
        }
        '''
        self.string_cog.set_options(options)
        self.string_cog.on('start', self.set_step, {'start': self.step.cog_annot})
        self.string_cog.on('end', self.set_step, {'end': self.step.cog_annot})
        self.string_cog.on('end', self.set_output, 'string_cog')
        self.string_cog.run()
        '''

    def run_blast_kegg(self):
        '''
        kegg 比对方法
        '''
        opts = self.blast_opts.copy()
        opts.update({"database": self.option("kegg_db")})
        opts.update({"version": self.option("version")})
        opts.update({"kegg_version": self.option("kegg_version")})
        opts.update({"diamond_version": self.option("diamond_version")})
        if self.option("diamond_version") >= "v2.0.13":
            opts.update({"sensitive": 1})
        genome_path = os.path.join(Config().SOFTWARE_DIR, "database/KEGG/genome2.xls")
        if self.option("kegg_species"):
            opts["database"] = 'customer_mode'
            genome_id = ''
            genome_abr = ''
            genome = ''
            with open (genome_path, 'r' ) as f:
                lines = f.readlines()
                for line in lines:
                    genome_id = re.sub("gn:", "", line.split("\t")[0])
                    genome_abr = line.split("\t")[1].split(',')[0]
                    genome = line.split("\t")[1].split(';')[-1].strip()
                    if self.option("kegg_species") == genome_abr or self.option("kegg_species") == genome:
                        break
            if genome_id:
                ref_path = os.path.join(Config().SOFTWARE_DIR, "database/KEGG/kegg_2017-05-01/kegg/genes/organisms", genome_abr, genome_id + '.pep' )
                opts['reference'] = ref_path

        for f in os.listdir(self.splitfasta.output_dir):
            opts['query'] = os.path.join(self.splitfasta.output_dir, f)
            if self.option('method') == "blast":
                blast_tool = self.add_tool('ref_rna_v2.annotation.blast')
            else:
                blast_tool = self.add_tool('ref_rna_v2.annotation.diamond')

            blast_tool.set_options(opts)
            self.blast_kegg_tools.append(blast_tool)
        if len(self.blast_kegg_tools) == 1:
            self.blast_kegg_tools[0].on("end", self.run_catblastout, 'kegg')
        else:
            self.on_rely(self.blast_kegg_tools, self.run_catblastout, 'kegg')
        for tool in self.blast_kegg_tools:
            tool.run()

    def run_kegg2ko(self):
        '''
        kegg xml转为ko的脚本暂时未实现
        '''

    def run_blast(self):
        '''
        blast比对步骤
        '''
        self.logger.info("blast 比对开始")
        database = self.option("database").split(";")
        for db in database:
            if db == 'nr':
                self.run_blast_nr()
            elif db == 'string':
                self.run_blast_string()
            elif db == 'swissprot':
                self.run_blast_swissprot()
            elif db == 'kegg':
                self.run_blast_kegg()
            elif db == 'pfam':
                pass
                # self.run_hmm_pfam()
            elif db == 'eggnog':
                self.run_blast_eggnog()
            else:
                self.logger.info("不支持该类型{}的数据库".format(db))

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
        #self.xmlfilter.on('end', self.set_output)
        self.xmlfilter.run()

    def run_catblastout(self, event):
        '''
        合并blast结果
        '''
        obj = event['bind_object']
        tools = []
        if event['data'] == 'nr':
            # self.set_step(event={'data': {'end': self.step.blast}})
            tools = self.blast_nr_tools
        elif event['data'] == 'string':
            tools = self.blast_string_tools
        elif event['data'] == 'swissprot':
            tools = self.blast_swissprot_tools
        elif event['data'] == 'eggnog':
            tools = self.blast_eggnog_tools
        elif event['data'] == 'kegg':
            tools = self.blast_kegg_tools

        if os.path.exists(self.work_dir + '/blast_tmp_' + event['data']):
            shutil.rmtree(self.work_dir + '/blast_tmp_' + event['data'])
        os.mkdir(self.work_dir + '/blast_tmp_' + event['data'])

        for i in tools:
            if len(os.listdir(i.output_dir)) != 1:
                self.logger.info(str(os.listdir(i.output_dir)))
                for f in os.listdir(i.output_dir):
                    if f.endswith("xml"):
                        file = f
            else:
                file = os.listdir(i.output_dir)[0]
            _path = os.path.join(i.output_dir, file)
            if os.path.exists(self.work_dir + '/blast_tmp_'+ event['data'] + "/" + os.path.basename(_path)):
                os.remove(self.work_dir + '/blast_tmp_'+ event['data'] + "/" + os.path.basename(_path))
            os.link(_path, self.work_dir + '/blast_tmp_'+ event['data'] + "/" + os.path.basename(_path))
        self.merge_tools[event['data']].set_options({"blastout": self.work_dir + '/blast_tmp_' + event['data']})
        self.merge_tools[event['data']].on('start', self.set_step, {'start': self.merge_step[event['data']]})
        self.merge_tools[event['data']].on('end', self.set_step, {'end': self.merge_step[event['data']]})
        self.merge_tools[event['data']].run()

    def run_catcog(self):
        '''
        合并cog结果 不在此处实现
        '''
        pass

    def run_catgo(self):
        '''
        合并GO结果文件
        '''
        if os.path.exists(self.work_dir + '/go_tmp'):
            shutil.rmtree(self.work_dir + '/go_tmp')
        os.mkdir(self.work_dir + '/go_tmp')

        for k,i in enumerate(self.nr2go_tools):
            if len(os.listdir(i.output_dir)) != 1:
                self.logger.info(str(os.listdir(i.output_dir)))
                for f in os.listdir(i.output_dir):
                    if f.endswith("blast2go_annot.xls"):
                        file = f
            else:
                file = os.listdir(i.output_dir)[0]
            _path = os.path.join(i.output_dir, file)
            if os.path.exists(self.work_dir + '/go_tmp' + "/" + os.path.basename(_path) + str(k)):
                os.remove(self.work_dir + '/go_tmp' + "/" + os.path.basename(_path) + str(k))
            os.link(_path, self.work_dir + '/go_tmp' + "/" + os.path.basename(_path)+ str(k))
        options = {
            "tabledir": self.work_dir + '/go_tmp',
            "outname": "blast2go_merge.xls"
        }
        seq_dict = self.option("query").get_all_seq_name()
        with open(self.work_dir + '/trans_list', 'w') as f:
            for seq_id in seq_dict.keys():
                f.write("{}\n".format(seq_id))

        if self.option("known_go"):
            options.update({
                "known": self.option("known_go"),
                "choose_known": self.work_dir +  "/trans_list",
                "known_format": "{}\t{}\tknown\tknown\t0\t100\t100\n",
                "known_fill": "1,2"
            })
        self.merge_tools['GO'].set_options(options)
        # self.catblast.on('start', self.set_step, {'start': self.step.cat_blastout})
        self.merge_tools['GO'].run()

    def run_catkegg(self):
        '''
        合并kegg结果
        '''
        pass

    def set_output(self):
        '''
        设置输出结果
        '''
        self.logger.info("event name is {}".format(self._rely.keys()))
        for file in os.listdir(self.output_dir):
            if os.path.isfile(self.output_dir + "/" + file) and os.path.exists(self.output_dir + "/" + file):
                os.remove(self.output_dir + "/" + file)
            elif os.path.isdir(self.output_dir + "/" + file) and os.path.exists(self.output_dir + "/" + file):
                shutil.rmtree(self.output_dir + "/" + file)
        for key, tool in self.merge_tools.items():
            shutil.copytree(tool.output_dir, self.output_dir + "/" + key, True)
        self.end()

    def run(self):
        # 设置运行逻辑
        super(AnnotMapdbModule, self).run()
        for db in self.option("database").split(";"):
            if db == 'nr':
                self.merge_tools['nr'] = self.add_tool("ref_rna_v2.annotation.cat_blastout")
                self.merge_tools['GO'] = self.add_tool("ref_rna_v2.annotation.cat_table")
            elif db == 'string':
                self.merge_tools['string'] = self.add_tool("ref_rna_v2.annotation.cat_blastout")
            elif db == 'swissprot':
                self.merge_tools['swissprot'] = self.add_tool("ref_rna_v2.annotation.cat_blastout")
            elif db == 'kegg':
                self.merge_tools['kegg'] = self.add_tool("ref_rna_v2.annotation.cat_blastout")
            elif db == 'eggnog':
                self.merge_tools['eggnog'] = self.add_tool("ref_rna_v2.annotation.cat_blastout")
            else:
                self.logger.info("不支持该类型{}的数据库".format(db))

        if len(self.merge_tools.values()) > 1:
            self.logger.info("等待以下数据库完成{}".format(" ".join(self.merge_tools.keys())))
            self.on_rely(self.merge_tools.values(), self.set_output)
        elif len(self.merge_tools.values()) == 1:
            self.merge_tools.values()[0].on('end', self.set_output)
        else:
            pass


        self.run_splitfasta()

    def end(self):
        repaths = [
            [".", "", "blast输出目录"],
            ["blast.xml", "xml", "blast xml输出结果文件"],
            ["blast_table.xls", "xls", "blast xls输出结果文件"]
        ]
        sdir = self.add_upload_dir(self.output_dir)
        sdir.add_relpath_rules(repaths)
        super(AnnotMapdbModule, self).end()

class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        import datetime

        test_dir = '/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/test_ref_rna_v2'
        data = {
            "id": "annot_db" + datetime.datetime.now().strftime('%H-%M-%S'),
            "type": "module",
            "name": "ref_rna_v2.annot_mapdb",
            "instant": False,
            "options": dict(
                query = test_dir + "/" + "test.cds.fa",
                # database="nr",
                method = "diamond",
                nr_db = "fungi",
                lines = 50
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()

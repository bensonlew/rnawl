# -*- coding: utf-8 -*-
# __author__ = 'gao.hao'
# last_modify: 2019.01.10

import os,shutil
from biocluster.core.exceptions import OptionError
from biocluster.module import Module
import datetime

class AnnotationModule(Module):
    """
    细菌基因组基础注释模块
    分析内容：NR、swiss-prot、pfam、cog、go和kegg注释, cazy
    """

    def __init__(self, work_id):
        super(AnnotationModule, self).__init__(work_id)
        option = [
            {"name": "gene_seq", "type": "infile", "format": "sequence.fasta"},  # 序列
            {"name": "gene_gff", "type": "infile", "format": "gene_structure.gff3"},  # 基因预测gff文件
            {"name": "database_list", "type": "string", "default": "nr_v20200604,card_v3.0.9,eggnog,kegg_v94.2,cazy_v8"},
            {"name": "sample", "type": "string", "default": ""},  # 样品名称
            {"name": "tidy_nr", "type": "outfile", "format": "sequence.profile_table"},
            {"name": "tidy_cog", "type": "outfile", "format": "sequence.profile_table"},
            {"name": "tidy_kegg", "type": "outfile", "format": "sequence.profile_table"},
            {"name": "tidy_cazy", "type": "outfile", "format": "sequence.profile_table"},
            {"name": "tidy_summary", "type": "outfile", "format": "sequence.profile_table"},
            {"name": "category", "type": "string", "default": "drug_class"},
        ]
        self.add_option(option)
        self.diamond = self.add_module('align.diamond')  # nr
        self.hmmscan = self.add_module('align.hmmscan')  # cazy
        self.top_diamond = self.add_module('align.meta_diamond')  # cog、kegg
        self.cog_anno = self.add_tool('annotation.mg_cog_anno')
        self.kegg_anno = self.add_tool('metagbin.kegg_anno')
        self.cazy_anno = self.add_tool('annotation.cazy_anno')
        self.anno_tidy = self.add_tool('metagbin.anno_tidy')
        self.card = self.add_module("metagbin.card_dna_anno")
        self.step.add_steps( 'kegg', 'anno_tidy','anno_cazy','card')
        self.tools = []
        self.database = ""

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def check_options(self):
        """
    检查参数
        :return:
        """
        if not self.option("gene_seq").is_set:
            raise OptionError("必须提供基因序列文件")
        if self.option("database_list") == "":
            raise OptionError("必须提供比对的数据名称")
        if self.option("sample") == "":
            raise OptionError("必须提供样品名称")
        return True

    def run_diamond_blast(self):
        opts = {
            'query': self.option('gene_seq'),
            'query_type': 'prot',
            'database': self.database,
            "blast": 'blastp',
            'outfmt': 6,
            'sensitive':0, ##qingchen.zhang
            'lines' : 50000  ##qingchen.zhang
        }
        if self.database in ["nr", 'nr_v20200604']:
            self.tools.append(self.diamond)
            self.diamond.set_options(opts)
            self.diamond.run()

    def run_top_dimond(self):
        opts = {
            'query': self.option('gene_seq'),
            'query_type': 'prot',
            'database': self.database,
            'outfmt': 5
        }
        if self.database == "eggnog":
            self.top_diamond.set_options(opts)
            self.tools.append(self.top_diamond)
            self.tools.append(self.cog_anno)
            self.top_diamond.on("end", self.run_cog_anno)
            self.top_diamond.run()
        elif self.database in ["kegg", "kegg_v94.2"]:
            self.top_diamond2 = self.add_module('align.meta_diamond')
            self.top_diamond2.set_options(opts)
            self.tools.append(self.top_diamond2)
            self.tools.append(self.kegg_anno)
            self.top_diamond2.on("end", self.run_kegg_anno)
            self.top_diamond2.run()

    def run_hmmscm(self):
        opts = {
            'query': self.option('gene_seq'),
            'database': self.database,
        }
        if self.database in ["cazy", "cazy_v8"]:
            self.hmmscan.set_options(opts)
            self.tools.append(self.hmmscan)
            self.tools.append(self.cazy_anno)
            self.hmmscan.on("end", self.run_cazy_anno)
            self.hmmscan.run()

    def run_cog_anno(self):
        opts = {
            'cog_xml': self.top_diamond.option('outxml')
        }
        self.cog_anno.set_options(opts)
        self.tools.append(self.cog_anno)
        self.cog_anno.run()

    def run_kegg_anno(self):
        opts = {
            'kegg_xml': self.top_diamond2.option('outxml'),
            'database': self.top_diamond2.option("database") ## fix by qingchen.zhang 兼容新老版本
        }
        self.kegg_anno.set_options(opts)
        self.tools.append(self.kegg_anno)
        self.kegg_anno.on('end', self.set_output, "kegg")
        self.kegg_anno.run()

    def run_cazy_anno(self):
        opts = {
            'hmmscan_result': self.hmmscan.option('align_result'),
            'version': self.hmmscan.option("database"), ## fix by qingchen.zhang 兼容新老版本
        }
        self.cazy_anno.set_options(opts)
        self.tools.append(self.cazy_anno)
        self.cazy_anno.on('end', self.set_output, "anno_cazy")
        self.cazy_anno.run()

    def run_card(self):
        opts = {
            "query": self.option("gene_seq"),
            "sample": self.option("sample"),
            'database': self.database, ## 增加此参数是为了区分新老版本
            'category': self.option("category")
        }
        self.card.set_options(opts)
        self.card.on('end', self.set_output, "card")
        self.card.run()

    def run_anno_tidy(self):
        opts = {
            'gene_gff': self.option('gene_gff'),
            'anno_nr': self.diamond.option('outtable'),
            'anno_cog': self.cog_anno.output_dir + "/gene_cog_anno.xls",
            'kegg_xml':self.top_diamond2.option('outxml'),
            'anno_kegg': self.kegg_anno.option('anno_kegg'),
            'anno_cazy': self.cazy_anno.output_dir + "/gene_cazy_parse_anno.xls",
            'summary': 'T',
            'sample': self.option("sample")
        }
        if self.option("database_list") in ['nr_v20200604,card_v3.0.9,eggnog,kegg_v94.2,cazy_v8']:
            opts["version"] = "new" ## 用于兼容新老版本的问题
        else:
            opts["version"] = "old"
        if os.path.exists(self.card.output_dir + "/" + self.option("sample") + "_card_anno.xls"):
            opts['anno_card'] = self.card.output_dir + "/" + self.option("sample") + "_card_anno.xls"
        self.anno_tidy.set_options(opts)
        self.anno_tidy.on('end', self.set_output, "anno_tidy")
        self.anno_tidy.run()

    def run(self):
        super(AnnotationModule, self).run()
        if self.option('database_list') != "":
            base = self.option('database_list').split(",")
            for ba in base:
                self.database = ba
                if self.database in ["nr", "nr_v20200604"]:
                    self.run_diamond_blast()
                if self.database in ["eggnog", "kegg", "kegg_v94.2"]:
                    self.run_top_dimond()
                if self.database in ["cazy", "cazy_v8"]:
                    self.run_hmmscm()
                if self.database in ["card", "card_v3.0.9"]:
                    self.run_card()
        self.on_rely(self.tools, self.run_anno_tidy)
        self.anno_tidy.on('end',self.end)

    def set_output(self, event):
        obj = event['bind_object']
        if event['data'] == 'anno_tidy':
            if os.path.exists(obj.output_dir + '/NR/'):
                self.linkdir(obj.output_dir + '/NR', 'NR')
                self.option('tidy_nr', self.anno_tidy.option('tidy_nr'))
            if os.path.exists(obj.output_dir + '/COG/'):
                self.linkdir(obj.output_dir + '/COG', 'COG')
                self.option('tidy_cog', self.anno_tidy.option('tidy_cog'))
            if os.path.exists(obj.output_dir + '/KEGG/'):
                # self.linkdir(obj.output_dir + '/KEGG', 'KEGG')  其中有一个KEGG通路图的文件夹，不允许这样操作
                if os.path.exists(self.output_dir + "/KEGG/"):
                    shutil.rmtree(self.output_dir + "/KEGG/")
                os.rename(obj.output_dir + '/KEGG', self.output_dir + "/KEGG/")
                os.link(self.kegg_anno.output_dir + "/kegg_level_stat.xls",
                        self.output_dir + "/KEGG/" + self.option("sample") + "_kegg_level_stat.xls")
                self.option('tidy_kegg', self.anno_tidy.option('tidy_kegg'))
            if os.path.exists(self.output_dir + '/CAZY/' + self.option("sample") + "_anno_cazy.xls"):
                os.remove(self.output_dir + '/CAZY/' + self.option("sample") + "_anno_cazy.xls")
            os.link(obj.output_dir + '/CAZy/' + self.option("sample") + "_anno_cazy.xls", self.output_dir + '/CAZY/' + self.option("sample") + "_anno_cazy.xls")
            self.option('tidy_cazy', self.anno_tidy.option('tidy_cazy'))
            if os.path.exists(obj.output_dir + '/Summary/'):
                self.linkdir(obj.output_dir + '/Summary', 'Summary')
                self.option('tidy_summary', self.anno_tidy.option('tidy_summary'))
        if event['data'] == 'anno_cazy':
            if not os.path.exists(self.output_dir + '/CAZY'):
                os.mkdir(self.output_dir + '/CAZY')
            for i in ['gene_cazy_family_stat.xls','gene_cazy_class_stat.xls','gene_cazy_parse_anno.xls']:
                j = i.replace('gene', self.option("sample"))
                if os.path.exists(self.output_dir + '/CAZY/' + j):
                    os.remove(self.output_dir + '/CAZY/' + j)
                os.link(self.cazy_anno.output_dir + '/' + i,self.output_dir + '/CAZY/' + j)
        if event['data'] == 'kegg':
            if os.path.exists(self.output_dir + '/' + 'kegg_pathway_img.tar.gz'):
                os.remove(self.output_dir + '/' + 'kegg_pathway_img.tar.gz')
            os.link(self.kegg_anno.output_dir + '/kegg_pathway_img.tar.gz',self.output_dir + '/' + 'kegg_pathway_img.tar.gz')
        if event['data'] == 'card':
            if len(os.listdir(self.card.output_dir)) > 0:
                self.linkdir(self.card.output_dir, self.output_dir + '/CARD')
                # os.rename(self.output_dir + '/CARD/' + self.option("sample") + "_card_category.xls", self.output_dir + '/CARD/' + self.option("sample") + "_card_drug_class.xls")

    def linkdir(self, dirpath, dirname):
        allfiles = os.listdir(dirpath)
        newdir = os.path.join(self.output_dir, dirname)
        if not os.path.exists(newdir):
            os.mkdir(newdir)
        oldfiles = [os.path.join(dirpath, i) for i in allfiles]
        newfiles = [os.path.join(newdir, i) for i in allfiles]
        for newfile in newfiles:
            if os.path.exists(newfile):
                os.remove(newfile)
        for i in range(len(allfiles)):
            os.link(oldfiles[i], newfiles[i])

    def end(self):
        super(AnnotationModule, self).end()

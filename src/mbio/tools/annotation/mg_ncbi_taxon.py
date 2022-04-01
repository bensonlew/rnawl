# -*- coding: utf-8 -*-
# __author__ = 'shenghe'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from biocluster.config import Config
import threading
import re, os
import shutil
from mbio.packages.annotation.mg_annotation.mg_taxon import mg_taxon
from pymongo.errors import ServerSelectionTimeoutError
from pymongo.errors import NetworkTimeout
import time


class MgNcbiTaxonAgent(Agent):
    """
    gi2taxon.py v1.0
    author: shenghe
    last_modify: 2016.6.21
    """

    def __init__(self, parent):
        super(MgNcbiTaxonAgent, self).__init__(parent)
        options = [
            {"name": "blastout", "type": "infile", "format": "align.blast.blast_xml, align.blast.blast_table"},  # 输入文件
            {"name": "taxon_out", "type": "outfile", "format": "taxon.seq_taxon"},  # 输出结果文件
            {"name": "blastdb", 'type': 'string', 'default': 'None'},  # 输入文件的blast比对类型，必须为nr或者nt
            {"name": "nr_method", "type": "string", "default": "best_hit"},  # 新增NR不同筛选结果方法,best_hit,lca,deunclassied  @20190326qingchen.zhang
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("blastout").is_set:
            raise OptionError("必须设置输入文件")
        if self.option('blastdb') == 'None':
            raise OptionError("必须设置输入文件的blast比对类型")
        else:
            if self.option('blastdb') not in ['nr', 'nt', "nt_v20200604", 'nr_v20200604','nt_v20210917']:
                raise OptionError('blastdb必须为nr或者nt:%s', variables=(self.option('blastdb')))
        return True

    def set_resource(self):
        self._cpu = 5
        self._memory = '20G'
        self._memory_increase_step = 1

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
            ['query_taxons_detail.xls', 'xls', '序列详细物种分类文件']
        ])
        super(MgNcbiTaxonAgent, self).end()


class MgNcbiTaxonTool(Tool):
    def __init__(self, config):
        super(MgNcbiTaxonTool, self).__init__(config)
        self._version = "1.0"
        self.client = Config().get_mongo_client(mtype="metagenomic", ref=True)
        self.mongodb = self.client[Config().get_mongo_dbname("metagenomic", ref=True)]
        self.gi_tax = self.mongodb.NR_sequence_20200604
        self.gene_tax = {}
        self.process_rerun = 0

    def run(self):
        """
        运行
        :return:
        """
        super(MgNcbiTaxonTool, self).run()
        if self.option("blastout").format == 'align.blast.blast_xml':
            if self.option('blastdb') in ['nt_v20200604','nt_v20210917']:
                from mbio.packages.align.blast.xml2table import xml2table
                table = xml2table(self.option('blastout').path, self.work_dir + '/temp_blastable.xls', hit_id="yes")
            else:
                from mbio.packages.align.blast.xml2table import xml2table
                table = xml2table(self.option('blastout').path, self.work_dir + '/temp_blastable.xls')
        else:
            table = self.option('blastout').path
        if self.option('blastdb') in ['nr', 'nr_v20200604']:
            db = 'prot'
        elif self.option('blastdb') in ['nt', "nt_v20200604","nt_v20210917"]:
            db = 'nucl'
        if self.process_run(table, db):
            if self.new_taxon():
                if self.option('blastdb') in ['nr', 'nr_v20200604']:
                    # gene_tax = self.get_taxon_id()
                    gene_tax = self.gene_tax
                    data_length, data_identity = self.get_anno_info(table)
                    self.get_anno_nr(gene_tax, data_length, data_identity)
                elif self.option('blastdb') in ['nt', "nt_v20200604","nt_v20210917"]:
                    if os.path.exists(self.work_dir + '/temp_blastable.xls'):
                        os.link(self.work_dir + '/temp_blastable.xls', self.output_dir + '/blast_table.xls')
                self.set_output()
                self.end()
            else:
                self.set_error('注释结果转换出错！')
        else:
            self.set_error('注释查询出错！')

    def process_run(self, fp, db):
        """单独的线程运行查询分类注释
        modified by xieshichang 20200623
        """
        try:
            if db == 'prot':
                if self.option('blastdb') in ['nr']:
                    from mbio.packages.taxon.gi2taxinfo import gi_taxon
                    query = self.filter_query(fp, db)
                    hitnames = set([q[1] for q in query])
                    id_taxons = gi_taxon(hitnames)
                elif self.option('blastdb') in ['nr_v20200604']:
                    from mbio.packages.taxon.name2taxinfo import name_taxon
                    query = self.filter_query(fp, db)
                    hitnames = set([q[1] for q in query])
                    id_taxons = name_taxon(hitnames)
            else:  # fix by qingchen.zhang @20200820 add nt database
                if self.option('blastdb') in ['nt']:
                    from mbio.packages.taxon.accession2taxon import accession_taxon  # import进入了sqlite3的对象，这个对象不可以跨线程使用
                elif self.option('blastdb') in ['nt_v20210917']:
                    from mbio.packages.taxon.accession2taxon_v20210917 import accession_taxon
                else:
                    from mbio.packages.taxon.accession2taxon_update import accession_taxon
                query = self.filter_query(fp, db)
                hitnames = set([q[1] for q in query])
                id_taxons = accession_taxon(hitnames)
        except (ServerSelectionTimeoutError, NetworkTimeout): # 捕获因为mongo服务器问题导致的异常后重运行此方法
            if self.process_rerun < 5:
                self.process_rerun += 1
                self.logger.info("检测到TimeoutError, 第{}次重运行方法".format(self.process_rerun))
                time.sleep(5)
                self.process_run(fp, db)
            else:
                self.add_state('memory_limit', '检测到TimeoutError, 重运行tool')
        for i in range(len(query)):
            if query[i][1] in id_taxons:
                self.logger.info(query[i][1])
                if self.option('blastdb') in ['nr', 'nr_v20200604']:
                    query[i].append(id_taxons[query[i][1]][1])
                    genes = query[i][0].split("_")
                    newgene = "_".join(genes[0:len(genes) - 1])
                    self.gene_tax[newgene] = id_taxons[query[i][1]][0]
                else:
                    query[i].append(id_taxons[query[i][1]])
        with open(self.work_dir + '/query_taxons_detail.xls', 'w') as w:
            for item in query:
                if len(item) == 3:
                    if item[2]:
                        w.write(item[0] + '\t' + item[1] + '\t' + item[2] + '\n')
                elif len(item) == 2:
                    if item[1]:
                        w.write(item[0] + '\t' + item[1] + '\t' + "root{no rank}" + '\n')
        return True

    def filter_query(self, fp, db):
        """生成query对gi的
        modified by xieshichang 20200623
        以列表保存
        """
        query = list()
        openfp = open(fp)
        openfp.readline()
        for i in openfp:
            line_sp = i.split('\t')
            hitname = line_sp[10].strip()
            # hitname = line_sp[10].split('|')
            # if len(hitname) != 1:  #新版本只有accession号
            #     if len(hitname) >= 3 and (hitname[0] != 'gi' or (not re.match(r'^\d+$', hitname[1]))):
            #         self.logger.info("hitname: {}".format(hitname))
            #         self.set_error('输入文件中不是nr库比对结果,不含有gi信息:%s')
            #     elif len(hitname) < 3:#新版本只有accession号
            #         self.logger.info('新版本只用accession_id')
            if db == 'prot':
                if self.option('blastdb') in ['nr_v20200604']:
                    query.append([line_sp[5], hitname])
                else:
                    hitname = hitname.split('|')
                    query.append([line_sp[5], hitname[1]])
            else:
                hitname = line_sp[10].split('|')
                query.append([line_sp[5], hitname[0]])
        return query

    def new_taxon(self):
        '''
         注释结果转换
        :return:
        '''
        tax = mg_taxon()
        query_index = {}
        self.logger.info("start new_taxon(detail_to_level)")
        try:
            if "nr_method" in self.get_option_object().keys() and self.option("nr_method") == "best_hit":
                self.query_index = tax.detail_to_level(self.work_dir + '/query_taxons_detail.xls', self.work_dir)
            elif "nr_method" in self.get_option_object().keys() and self.option("nr_method") == "deunclassied":
                self.query_index = tax.detail_to_level_deunclassified(self.work_dir + '/query_taxons_detail.xls', self.work_dir)
            elif "nr_method" in self.get_option_object().keys() and self.option("nr_method") == "lca":
                self.query_index = tax.detail_to_level_lca(self.work_dir + '/query_taxons_detail.xls', self.work_dir)
            else:
                self.query_index = tax.detail_to_level(self.work_dir + '/query_taxons_detail.xls', self.work_dir)
        except Exception as e:
            self.set_error("new_taxon(detail_to_level) failed")
            raise Exception("new_taxon(detail_to_level) failed {}".format(e))
        return True

    def get_taxon_id(self):
        """
        根据NR注释结果得到注释结果和id关联
        :return:
        """
        gene_tax = {}
        with open(self.work_dir + '/query_taxons_detail.xls', "r") as f1:
            for line in f1:
                line = line.strip().split("\t")
                gene = line[0]
                genes = gene.split("_")
                newgene = "_".join(genes[0:len(genes) - 1])
                gi = str(line[1])
                detail = self.gi_tax.find_one({"origin_sequence_id": gi})
                if detail:
                    taxid = detail["taxid"]
                    gene_tax[newgene] = taxid
        return gene_tax

    def get_anno_info(self, table):
        """
        根据注释的identity和length
        :return:
        """
        data_length = {}
        data_identity = {}
        tmp_index = {}
        with open(table, "r") as f:
            head = f.next().strip()
            for line in f:
                line = line.strip()
                line1 = line.split("\t")
                if line1[0] != "Score":
                    identity = line1[3]
                    length = line1[2]
                    gene = line1[5]
                    new = gene.split("_")
                    newgene = "_".join(new[0:len(new) - 1])
                    if newgene not in data_length:
                        data_identity[newgene] = 0.0
                        data_length[newgene] = 0.0
                        tmp_index[gene] = 0
                    if "nr_method" in self.get_option_object().keys() and self.option("nr_method") == "lca":
                        data_length[newgene] += float(length) / self.query_index[gene]
                        data_identity[newgene] += float(identity) / self.query_index[gene]
                    elif tmp_index[gene] == self.query_index[gene]:
                        data_length[newgene] = length
                        data_identity[newgene] = identity
                    tmp_index[gene] += 1
            del tmp_index
        return data_length, data_identity

    def get_anno_nr(self, gene_tax, data_length, data_identity):
        """
        根据比对结果和获得的分类地位得到注释结果表
        :return:
        """
        with open(self.work_dir + "/query_taxons.xls", "r") as f, open(self.work_dir + "/gene_nr_anno.xls", "w") as outfile:
            #outfile.write("#Query\tTaxid\tDomain\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\tIdentity(%)\tAlign_len\n")
            for line in f:
                line = line.strip().split("\t")
                gene = line[0]
                tax = "\t".join(line[1].replace(" ","_").split(";"))
                if gene_tax.has_key(gene):
                    outfile.write(gene + "\t" + str(gene_tax[gene]) + "\t" + tax + "\t" + str(data_length[gene]) + "\t" + str(data_identity[gene]) + "\n")

    def set_output(self):
        """
        设置结果目录
        :return:
        """
        if os.path.exists(self.output_dir + "/query_taxons.xls"):
            os.remove(self.output_dir + "/query_taxons.xls")
        os.link(self.work_dir + "/query_taxons.xls", self.output_dir + "/query_taxons.xls")
        self.option('taxon_out', self.output_dir + '/query_taxons.xls')
        if os.path.exists(self.output_dir + "/gene_nr_anno.xls"):
            os.remove(self.output_dir + "/gene_nr_anno.xls")
        if os.path.exists(self.work_dir + "/gene_nr_anno.xls"):
            os.link(self.work_dir + "/gene_nr_anno.xls", self.output_dir + "/gene_nr_anno.xls")
        self.logger.info("Have got correct results")

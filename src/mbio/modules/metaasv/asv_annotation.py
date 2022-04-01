# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

from biocluster.module import Module
import os
import shutil
from biocluster.core.exceptions import OptionError
from mbio.packages.metaasv.common_function import link_dir, link_file
from biocluster.config import Config
from mbio.packages.annotation.mg_annotation.mg_taxon import mg_taxon



class AsvAnnotationModule(Module):
    """
    metaasv 注释统计
    方法有五种:blast、multi_blast、vsearch、bayes、rdp
    比对+注释+统计
    """
    def __init__(self, work_id):
        super(AsvAnnotationModule, self).__init__(work_id)
        options = [
            {'name': 'fasta', 'type': 'infile', 'format': 'sequence.fasta'},  # 输入fasta文件
            {'name': 'in_otu_table', 'type': 'infile', 'format': 'meta.otu.otu_table'},  # 输入的asv表
            # {'name': 'qza_otu_table', 'type': 'infile', 'format': 'metaasv.qza'},  # 输入asv文件的qza格式
            {'name': 'qza_fasta', 'type': 'infile', 'format': 'metaasv.qza'},  # 输入fasta文件的qza格式
            {'name': 'confidence', 'type': 'float', 'default': 0.7},  # 置信度值只针对用qiime比对的序列
            {'name': 'database', 'type': 'string'},  # 数据库选择
            {'name': 'ref_fasta', 'type': 'infile', 'format': 'sequence.fasta'},  # 参考fasta序列
            {'name': 'ref_taxon', 'type': 'infile', 'format': 'taxon.seq_taxon'},  # 参考taxon文件
            {'name': 'anno_method', 'type': 'string'},  # 注释方法
            {"name": "query_type", "type": "string", "default": "nucl"},# 输入的查询序列的格式，为nucl或者prot
            {'name': 'blast', 'type': 'string', "default": "blastn"},  # blast选择的方式blastn等方法
            {'name': 'identity', 'type': 'int', "default": 80},  # blast的比对identity
            {'name': 'coverage', 'type': 'int', "default": 80},  # blast的比对coverage
            {'name': 'taxon_file', 'type': 'outfile', 'format': 'taxon.seq_taxon'},  # 输出序列的分类信息文件
            {'name': 'asv_md5', 'type': 'infile', 'format': 'sequence.profile_table'},  # 二维表格
            {"name": "evalue", "type": "float", "default": 1e-5},
            {"name": "meta_pipeline", "type": "string", "default": ""},
        ]
        self.add_option(options)
        self.qiime2 = self.add_module("metaasv.qiime2_annotation")
        self.blast = self.add_module("metaasv.blast_annotation")
        self.blastp = self.add_module("align.blast")
        self.rdp = self.add_module("annotation.meta_tax")
        self.stat = self.add_tool("metaasv.asv_taxon_stat")
        self.diamond = self.add_module("metaasv.diamond2")

    def check_options(self):
        """
        参数二次检查
        :return:
        """
        if not self.option("fasta").is_set and not self.option("qza_fasta"):
            raise OptionError("必须设置参数fasta")
        if self.option('database') == "custom_mode":
            if not self.option("ref_fasta").is_set or not self.option("ref_taxon").is_set:
                raise OptionError("数据库自定义模式必须设置参考fasta序列和参考taxon文件")
        else:
            if self.option("database") not in ['unite8.0/its_fungi', 'unite7.2/its_fungi', 'unite7.0/its_fungi',
                                               'fgr/amoA', 'fgr/nosZ', 'fgr/nirK', 'fgr/nirS',
                                               'fgr/nifH', 'fgr/pmoA', 'fgr/mmoX', 'fgr/mcrA',
                                               'fgr/amoA_archaea', 'fgr/amoA_bacteria',
                                               'maarjam081/AM', 'Human_HOMD', 'Human_HPB', 'Protist_PR2_v4.5',
                                               'silva132/16s_archaea', 'silva132/16s_bacteria',
                                               'silva132/18s_eukaryota', 'silva132/16s',
                                               'silva138/16s_archaea', 'silva138/16s_bacteria',
                                               'silva138/18s_eukaryota', 'silva138/16s',
                                               'greengenes135/16s', 'greengenes135/16s_archaea',
                                               'greengenes135/16s_bacteria',
                                               'rdp11.5/16s', 'rdp11.5/16s_bacteria', 'rdp11.5/16s_archaea',
                                               'nt_v20200327/16s_archaea', 'nt_v20200327/16s_bacteria','nt_v20200327/16s',
                                               'Human_HOMD_v15.2', 'nt_v20200327/18s_eukaryota', 'nt_v20200327/its_fungi',
                                               'nt_v20210917/16s_archaea', 'nt_v20210917/16s_bacteria', 'nt_v20210917/16s',
                                               'nt_v20210917/18s_eukaryota', 'nt_v20210917/its_fungi','nt_v20210917',
                                               'nt', "nt_v20200604",'fgr/amoA_archaea_202012', 'fgr/amoA_bacteria_202012',
                                               'fgr/amoA_AOB_like_202012', 'fgr/amoA_comammox_202012', 'fgr/nosZ_202012',
                                               'fgr/nosZ_atypical_1_202012', 'fgr/nosZ_atypical_2_202012', 'fgr/nirK_202012',
                                               'fgr/nirS_202012', 'fgr/mcrA_202012', 'fgr/nifH_202012', 'fgr/pmoA_202012','fgr/mmoX_202012',
                                               "nr_v20210917"]:
                # 多样性目前的database
                raise OptionError("数据库%s不被支持", variables=(self.option("database")))
        if self.option("database") in ['nt', "nt_v20200604"]:
            if self.option('query_type') not in ['nucl', 'prot']:
                raise OptionError('query_type查询序列的类型为nucl(核酸)或者prot(蛋白):%s', variables=(self.option('query_type')))
            #else:
            #    if self._fasta_type[self.option('fasta').prop['seq_type']] != self.option('query_type'):
            #        raise OptionError('文件检查发现查询序列为: %s,而需要的文件类型为:%s', variables=(
            #            self._fasta_type[self.option('fasta').prop['seq_type']],self.option('query_type')))

                #raise OptionError('文件类型错误，请检查文件是否是核酸或者蛋白文件！')
        return True

    def run_rdp(self):
        """
        用rdp方法注释，不做nt库
        :return:
        """
        options = {
            "fasta": self.option("fasta"),
            "revcomp": "False",
            "confidence": self.option("confidence"),
            "database": self.option("database")
        }
        if self.option("database") == "custom_mode":
            options.update({
                "ref_fasta": self.option("ref_fasta"),
                "ref_taxon": self.option("ref_taxon")
            })
        self.rdp.set_options(options)
        self.rdp.run()

    def run_qiime2(self):
        """
        用blast、vsearch、bayes方法计算
        :return:
        """
        options = {
            "database": self.option("database"),
            "anno_method": self.option("anno_method")
        }
        if self.option("qza_fasta").is_set:
            options["qza_fasta"] = self.option("qza_fasta")
        if self.option("fasta").is_set:
            options["fasta"] = self.option("fasta")
        if self.option("anno_method") in ["bayes"]:
            options["confidence"] = self.option("confidence")
        else:
            options["identity"] = float(self.option("identity")) / 100
            options["coverage"] = float(self.option("coverage")) / 100
        if self.option("database") == "custom_mode":
            options.update({
                "ref_fasta": self.option("ref_fasta"),
                "ref_taxon": self.option("ref_taxon")
            })
        self.qiime2.set_options(options)
        self.qiime2.run()

    def run_blast(self):
        """
        用multi_blast方法注释，用nt和所有数据库做；
        :return:
        """
        options = {
            "query": self.option("fasta"),
            "database": self.option("database"),
            "identity": self.option("identity"),
            "coverage": self.option("coverage"),
            "blast": self.option("blast"),
            "query_type": self.option("query_type"),

        }
        if self.option("database") == "custom_mode":
            options.update({
                "ref_fasta": self.option("ref_fasta"),
                "ref_taxon": self.option("ref_taxon"),
                "reference_type": "nucl"
            })
        self.blast.set_options(options)
        self.blast.run()

    def run_blastp(self):
        """
        功能基因blastp注释；
        :return:
        """
        options = {
            "query": self.option("fasta"),
            "database": self.option("database"),
            "evalue": self.option("evalue"),
            "query_type": self.option("query_type"),
            "reference_type": "prot",
            "blast": "blastp"
        }
        if self.option("database") == "custom_mode":
            options.update({
                "ref_fasta": self.option("ref_fasta"),
                "ref_taxon": self.option("ref_taxon"),
                "evalue": self.option("evalue"),
                "reference_type": "prot"
            })
        self.blast.set_options(options)
        self.blast.run()

    def run_blastp2(self):
        """
        功能基因流程的blastp注释；
        :return:
        """

        options = {
            "query": self.option("fasta"),
            "database": self.option("database"),
            "evalue": self.option("evalue"),
            "query_type": self.option("query_type"),
            "reference_type": "prot",
            "blast": "blastp",
            "outfmt": 5
        }
        if self.option("database") == "custom_mode":
            options.update({
                "reference": self.option("ref_fasta"),
                "database": "customer_mode",
                #"ref_taxon": self.option("ref_taxon"),
                "evalue": self.option("evalue"),
                "reference_type": "prot"
            })
        self.blastp.set_options(options)
        self.blastp.run()

    def run_diamond2(self):
        self.diamond.set_options({
            "query": self.option("fasta"),
            "query_type": self.option("query_type"),
            "database": self.option("database"),
            "blast": self.option("blast"),
            "evalue": self.option("evalue")
        })
        self.diamond.run()

    def run_taxon(self):
        """
        nt注释比对结果
        :return:
        """
        self.ncbi_taxon = self.add_tool("annotation.mg_ncbi_taxon")
        self.ncbi_taxon.set_options({
            "blastout": self.diamond.option("outxml"),
            #"out_type": 1,
            "blastdb": self.option("database")
        })
        self.ncbi_taxon.on("end", self.set_output)
        self.ncbi_taxon.run()

    def run_taxon_stat(self):
        """
        对注释结果进行统计
        :return:
        """
        os.system("dos2unix {}".format(self.option("in_otu_table").prop["path"]))
        if self.option("anno_method") in ["rdp"]:
            taxon_file = self.rdp.option("taxon_file").prop['path']
        elif self.option("anno_method") in ["blast", "vsearch", "bayes"]:
            with open(os.path.join(self.work_dir, "ASV_tax_assignments.txt"), 'r') as f, open(os.path.join(self.work_dir, "ASV_tax_assignments2.txt"), 'w') as w:
                for line in f:
                    if line.strip().split("\t")[0] in ["Feature ID"]:
                        pass
                    else:
                        w.write(line)
            taxon_file = os.path.join(self.work_dir, "ASV_tax_assignments2.txt")
        elif self.option("database") in ["nr_v20210917"]:
            taxon_file = os.path.join(self.work_dir, "ASV_tax_assignments.txt")
            from mbio.packages.align.blast.xml2table import xml2table
            table = xml2table(self.diamond.output_dir + "/blast.xml",
                              self.work_dir + '/temp_blastable.xls')
            query = self.filter_query(table)
            from mbio.packages.taxon.nr2taxon_v20210917 import accession_taxon
            id_taxons = accession_taxon(set(query.values()))
            for i in query:
                query[i] = (query[i], id_taxons[query[i]])
            with open(self.work_dir + '/query_taxons_detail.xls', 'w') as w:
                for item in query.iteritems():
                    if item[1][1]:
                        w.write(item[0] + '\t' + item[1][0] + '\t' + item[1][1] + '\n')
            tax = mg_taxon()
            tax.detail_to_level(self.work_dir + '/query_taxons_detail.xls', self.work_dir)
            os.rename(self.work_dir + '/query_taxons.xls',taxon_file)
        else:
            if self.option("meta_pipeline") == "metaasv":
                taxon_file = os.path.join(self.work_dir, "ASV_tax_assignments.txt")
                if self.option("database") == "custom_mode":
                    db_path = self.option("ref_taxon").prop["path"]
                else:
                    db_path = Config().SOFTWARE_DIR + "/database/Framebot/tax/" + self.option("database").split("/")[1] + ".tax"
                from mbio.packages.align.blast.xml2table import xml2table
                self.logger.info(self.blastp.option("outxml"))
                table = xml2table(self.blastp.option("outxml").prop["path"], self.work_dir + '/temp_blastable.xls')
                with open(table, "r") as f, open(db_path, "r") as v, open(taxon_file,"w") as t:
                    data1 = f.readlines()
                    data2 = v.readlines()
                    tax_dict = {}
                    for x in data2:
                        tax_dict[x.strip().split("\t")[0]] = x.strip().split("\t")[1]
                    for i in data1[1:]:
                        if i.strip().split()[10] in tax_dict:
                            t.write(i.strip().split()[5] + "\t" + tax_dict[i.strip().split()[10]] + "\n")
            else:
                taxon_file = os.path.join(self.blast.output_dir, "ASV_tax_assignments.txt")
        self.stat.set_options({
            "in_otu_table": self.option("in_otu_table"),
            "taxon_file": taxon_file
        })
        self.stat.on("end", self.set_output)
        self.stat.run()

    def replace_name(self):
        """
        替换qiime2的md5值的结果，将未注释到的ASV归为Unassigned的情况，并将其过滤掉进行统计
        :return:
        """
        md5_dict = {}
        md5_dict2 = {}
        if self.option("asv_md5").is_set:
            self.logger.info("有md5文件")
            with open(self.option("asv_md5").prop["path"], 'r') as f:
                lines = f.readlines()
                for line in lines[1:]:
                    line = line.strip().split("\t")
                    if len(line[1]) <= 10:
                        md5_dict[line[0]] = line[1]
                    else:
                        md5_dict[line[1]] = line[0]
                        md5_dict2[line[0]] = line[1]
            self.logger.info(md5_dict)
        else:
            self.logger.info("没有md5文件")
            with open(self.option("in_otu_table").prop["path"], 'r') as f:
                lines = f.readlines()
                for line in lines[1:]:
                    line = line.strip().split("\t")
                    md5_dict[line[0]] = line[0]
        qiime2_outfile = os.path.join(self.qiime2.output_dir, "ASV_tax_assignments.txt")
        outfile = os.path.join(self.work_dir, "ASV_tax_assignments.txt")
        with open(qiime2_outfile, 'r') as m, open(outfile, 'w') as w:
            lines = m.readlines()
            for line in lines[1:]:
                line = line.strip().split("\t")
                # line[1] = "; ".join(line[1].strip().split(";"))
                md5_name = line[0]
                if md5_name in md5_dict:
                    asv_name = md5_dict[md5_name]
                    if line[1] not in ["Unassigned", "Unassigned;"]:
                        w.write("{}\t{}\n".format(asv_name, "\t".join(line[1:])))
                    else:
                        w.write("{}\t{}\n".format(md5_name, "d__Unassigned"))
                elif md5_name in md5_dict2:
                    asv_name = md5_dict2[md5_name]
                    if line[1] not in ["Unassigned", "Unassigned;"]:
                        w.write("{}\t{}\n".format(asv_name, "\t".join(line[1:])))
                    else:
                        w.write("{}\t{}\n".format(md5_name, "d__Unassigned"))
                else:
                    if not self.option("anno_method") in ["blast", "multi_blast"]:
                        raise Exception("请检查上传的fasta文件中的序列和丰度表中的序列ID是否一致！")
                    else:
                        w.write("{}\t{}\n".format(md5_name, "d__Unassigned"))
        self.run_taxon_stat()

    def set_output(self):
        """
        设置结果文件目录
        :return:
        """
        taxon_dir = os.path.join(self.output_dir, "Tax_assign")
        if os.path.exists(taxon_dir):
            shutil.rmtree(taxon_dir)
        os.mkdir(taxon_dir)
        result = os.path.join(taxon_dir, 'ASV_tax_assignments.txt')
        if self.option("anno_method") in ["rdp"]:
            link_file(self.rdp.option("taxon_file").prop['path'], result)
            self.option("taxon_file", result)
        else:
            if self.option("anno_method") in ["blast", "vsearch", "bayes"]:
                link_dir(self.qiime2.output_dir, taxon_dir)
                link_file(os.path.join(self.work_dir, "ASV_tax_assignments.txt"), result)
                self.option("taxon_file", result)
            else:
                if self.option("meta_pipeline") == "metaasv":
                    link_file(os.path.join(self.work_dir, "ASV_tax_assignments.txt"), result)
                else:
                    link_file(os.path.join(self.blast.output_dir, "ASV_tax_assignments.txt"), result)

        result_dir = os.path.join(self.output_dir, 'ASVTaxon_summary')
        link_dir(self.stat.output_dir, result_dir)
        self.end()

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

    def run(self):
        """
        运行和设置逻辑
        :return:
        """
        super(AsvAnnotationModule, self).run()
        if self.option("anno_method") in ["blast", "vsearch", "bayes"]:
            self.qiime2.on("end", self.replace_name)
            self.run_qiime2()
        elif self.option("anno_method") in ["rdp"]:
            self.rdp.on("end", self.run_taxon_stat)
            self.run_rdp()
        elif self.option("anno_method") in ["blastp"]:
            if self.option("meta_pipeline") == "metaasv":
                self.blastp.on("end", self.run_taxon_stat)
                self.run_blastp2()
            else:
                self.blast.on("end", self.run_taxon_stat)
                self.run_blastp()
        elif self.option("anno_method") in ["diamond2"]:
            self.diamond.on("end", self.run_taxon_stat)
            self.run_diamond2()
        else:##multi_blast
            self.blast.on("end", self.run_taxon_stat)
            self.run_blast()

    def end(self):
        """
        结束
        :return:
        """
        super(AsvAnnotationModule, self).end()

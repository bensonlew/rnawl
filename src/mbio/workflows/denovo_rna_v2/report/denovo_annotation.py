# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'
# last_modifiy = modified 2017.12.23

from biocluster.workflow import Workflow
from biocluster.config import Config
from bson.son import SON
from bson.objectid import ObjectId
from mbio.workflows.denovo_rna_v2.denovo_test_api import DenovoTestApiWorkflow
from biocluster.wsheet import Sheet
import os
import re
import shutil
import json
import glob
from biocluster.config import Config
from mbio.packages.denovo_rna_v2.copy_file import CopyFile
from mbio.packages.project_demo.run_log.get_run_log import GetRunLog
from biocluster.core.function import filter_error_info, link, CJsonEncoder
from mbio.packages.denovo_rna_v2.chart import Chart


class DenovoAnnotationWorkflow(Workflow):
    """
    交互分析进行筛选blast参数重注释时使用
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        #self.rpc = False
        super(DenovoAnnotationWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "nr_evalue", "type": "float", "default": 1e-5},
            {"name": "nr_similarity", "type": "float", "default": 0},
            {"name": "nr_identity", "type": "float", "default": 0},
            {"name": "swissprot_evalue", "type": "float", "default": 1e-5},
            {"name": "swissprot_similarity", "type": "float", "default": 0},
            {"name": "swissprot_identity", "type": "float", "default": 0},
            {"name": "cog_evalue", "type": "float", "default": 1e-5},
            {"name": "cog_similarity", "type": "float", "default": 0},
            {"name": "cog_identity", "type": "float", "default": 0},
            {"name": "kegg_evalue", "type": "float", "default": 1e-5},
            {"name": "kegg_similarity", "type": "float", "default": 0},
            {"name": "kegg_identity", "type": "float", "default": 0},
            {"name": "pfam_evalue", "type": "float", "default": 1e-5},
            {"name": "stat_id", "type": "string"},
            {"name": "task_id", "type": "string"},
            {"name": "origin_result", "type": "infile", "format": "denovo_rna_v2.common_dir"},
            {"name": "last_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "trans2gene", "type": "infile", "format": "denovo_rna_v2.common"},
            {"name": "gene_exp", "type": "infile", "format": "denovo_rna_v2.common"},
            {"name": "trans_exp", "type": "infile", "format": "denovo_rna_v2.common"},
            {"name": "group_dict", "type": "string"},
            {"name": "taxonomy", "type": "string", "default": None},
            {"name": "exclude_taxon", "type": "string", "default": None},
            {"name": "origin_param", "type": "string", "default": None},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.nr_filter = self.add_tool("denovo_rna_v2.filter_annot")
        self.swissprot_filter = self.add_tool("denovo_rna_v2.filter_annot")
        self.cog_filter = self.add_tool("denovo_rna_v2.filter_annot")
        self.kegg_filter = self.add_tool("denovo_rna_v2.filter_annot")
        self.pfam_filter = self.add_tool("denovo_rna_v2.filter_annot")
        # self.annotation = self.add_module("denovo_rna_v2.denovo_annotation")
        self.annotation = self.add_module("denovo_rna_v2.annot_class_beta")
        self.nr2go = self.add_tool('ref_rna_v2.annotation.nr2go')
        self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/01 Annotation')
        self.inter_dirs = []

    def send_log(self, data):
        # 中间目录修改
        m = re.match("^([\w\-]+)://(.*)interaction_result.*$", self._sheet.output)
        region = m.group(1)
        inter_dir = m.group(2)
        self.logger.info("更新结果目录")

        if "dirs" in data["data"]["sync_task_log"]:
            for dir_path in self.inter_dirs:
                dir_dict = {
                    "path": os.path.join(inter_dir, "interaction_results", dir_path[0]),
                    "size": "",
                    "format": dir_path[1],
                    "description": dir_path[2],
                    "region": region,
                }
                if len(dir_path) >= 5:
                    dir_dict.update({"code": "D" + dir_path[5]})

                data["data"]["sync_task_log"]["dirs"].append(dir_dict)
        with open(self.work_dir + "/post.changed.json", "w") as f:
            json.dump(data, f, indent=4, cls=CJsonEncoder)
        super(DenovoAnnotationWorkflow, self).send_log(data)

    def run_nr_blast_filter(self):
        xml = glob.glob(os.path.join(self.option("origin_result").prop['path'], "blast_xml/*_vs_nr.xml"))
        options = {
            'xml': xml[0],
            'types': "xml",
            'evalue': self.option('nr_evalue'),
            'identity': self.option('nr_identity'),
            'similarity': self.option('nr_similarity')
        }
        if self.option('exclude_taxon'):
            options.update({
                'exclude_taxon': self.option('exclude_taxon')
            })
        self.nr_filter.on('end', self.run_nr2go)
        self.nr_filter.set_options(options)
        self.nr_filter.run()

    def run_swissprot_blast_filter(self):
        xml = glob.glob(os.path.join(self.option("origin_result").prop['path'], "blast_xml/*_vs_swissprot.xml"))
        options = {
            'xml': xml[0],
            'types': "xml",
            'evalue': self.option('swissprot_evalue'),
            'identity': self.option('swissprot_identity'),
            'similarity': self.option('swissprot_similarity')
        }
        self.swissprot_filter.set_options(options)
        self.swissprot_filter.run()

    def run_cog_blast_filter(self):
        xml = glob.glob(os.path.join(self.option("origin_result").prop['path'], "blast_xml/*_vs_string.xml"))
        if len(xml) >= 1:
            pass
        else:
            xml = glob.glob(os.path.join(self.option("origin_result").prop['path'], "blast_xml/*_vs_eggnog.xml"))
        options = {
            'xml': xml[0],
            'types': "xml",
            'evalue': self.option('cog_evalue'),
            'identity': self.option('cog_identity'),
            'similarity': self.option('cog_similarity')
        }
        self.cog_filter.set_options(options)
        self.cog_filter.run()

    def run_kegg_blast_filter(self):
        xml = glob.glob(os.path.join(self.option("origin_result").prop['path'], "blast_xml/*_vs_kegg.xml"))
        options = {
            'xml': xml[0],
            'types': "xml",
            'evalue': self.option('kegg_evalue'),
            'identity': self.option('kegg_identity'),
            'similarity': self.option('kegg_similarity')
        }
        self.kegg_filter.set_options(options)
        self.kegg_filter.run()

    def run_pfam_filter(self):

        options = {
            'hmm': os.path.join(self.option("origin_result").prop['path'], "blast_xml/pfam_domain"),
            'types': "hmm",
            'evalue': self.option('pfam_evalue'),
        }
        self.pfam_filter.set_options(options)
        self.pfam_filter.run()

    def run_nr2go(self):
        options = {
            'blastout': self.nr_filter.option('outxml').prop['path'],
        }
        self.nr2go.set_options(options)
        self.nr2go.run()

    def run_annotation_v2(self):
        if os.path.exists(self.option("trans2gene").prop['path']):
            trans2gene = self.option("trans2gene").prop['path']
        else:
            trans2gene = self.option("origin_result").prop['path']  + '/all_tran2gene.txt'
        options = {
            'db': 'nr,swissprot,kegg,eggnog,pfam,go',
            'type': 'new',
            'gene2trans':  trans2gene,
            'blast_nr_xml': self.nr_filter.option('outxml').prop['path'],
            'blast_swissprot_xml': self.swissprot_filter.option('outxml').prop['path'],
            'blast_eggnog_xml': self.cog_filter.option('outxml').prop['path'],
            'blast_kegg_xml': self.kegg_filter.option('outxml').prop['path'],
            'taxonomy': self.option('taxonomy'),
            'link_bgcolor': 'green',
            'png_bgcolor': '00CD00',
            'blast2go_annot': self.nr2go.output_dir + "/blast2go_annot.xls",
            'pfam_domain': self.pfam_filter.option('outtable').prop['path'],
        }

        self.annotation.set_options(options)
        self.annotation.on('end', self.set_output)
        self.annotation.run()

    def run_annotation(self):
        options = {
            'gene2trans': self.option('trans2gene'),
            'go_annot': True,
            'blast_nr_xml': self.nr_filter.option('outxml').prop['path'],
            'nr_annot': True,
            'blast_kegg_xml': self.kegg_filter.option('outxml').prop['path'],
            'taxonomy': self.option('taxonomy'),
            'blast_string_xml': self.cog_filter.option('outxml').prop['path'],
            'blast_swissprot_xml': self.swissprot_filter.option('outxml').prop['path'],
            'pfam_domain': self.pfam_filter.option('outtable').prop['path']
        }
        self.annotation.set_options(options)
        self.annotation.on('end', self.set_output)
        self.annotation.run()

    def run(self):
        self.on_rely([self.nr_filter, self.swissprot_filter, self.cog_filter, self.kegg_filter, self.pfam_filter, self.nr2go], self.run_annotation_v2)
        self.get_run_log()
        self.run_nr_blast_filter()
        self.run_cog_blast_filter()
        self.run_kegg_blast_filter()
        self.run_swissprot_blast_filter()
        self.run_pfam_filter()

        # self.annotation.output_dir = "/mnt/ilustre/users/sanger-dev/workspace/20171224/DenovoAnnotation_test_anno_web_1224_113531/DenovoAnnotation/output"
        # self.set_output()
        super(DenovoAnnotationWorkflow, self).run()

    def get_run_log(self):
        get_run_log = GetRunLog("denovo_rna_v2", table="sg_annotation_stat", main_id=self.option('stat_id'),
                                dir_path=self.work_dir)
        self.run_log = get_run_log.run()

    def set_output(self):
        self.logger.info("结果路径为{}".format(self.annotation.output_dir))
        #output_dir = self.annotation.output_dir
        self.logger.info("结果路径为{}".format(self.annotation.output_dir))
        self.set_db()

    def linkdir(self, olddir, newname, mode='link'):
        """
        移动目录下的输出文件/文件夹到输出文件夹下
        """
        if not os.path.isdir(olddir):
            self.set_error('需要移动到output目录的文件夹不存在。', code = "12000501")
        newdir = os.path.join(self.output_dir, newname)
        if os.path.exists(newdir):
            shutil.rmtree(newdir)
        os.mkdir(newdir)
        allfiles = os.listdir(olddir)
        oldfiles = [os.path.join(olddir, i) for i in allfiles]
        newfiles = [os.path.join(newdir, i) for i in allfiles]
        for i in range(len(allfiles)):
            if os.path.isfile(oldfiles[i]):
                os.link(oldfiles[i], newfiles[i])
            else:
                new1 = os.path.join(newdir, os.path.basename(oldfiles[i]))
                os.system("mv {} {}".format(oldfiles[i], new1))


    def set_db(self):
        self.logger.info("保存结果到mongo")
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_AFTER_END = False
        self.test_api = self.api.api("denovo_rna_v2.denovo_annotation")
        params = {
            "nr_evalue": self.option("nr_evalue"),
            "nr_similarity": self.option("nr_similarity"),
            "nr_identity": self.option("nr_identity"),
            "swissprot_evalue":self.option("swissprot_evalue"),
            "swissprot_similarity": self.option("swissprot_similarity"),
            "swissprot_identity": self.option("swissprot_identity"),
            "cog_evalue": self.option("cog_evalue"),
            "cog_similarity": self.option("cog_similarity"),
            "cog_identity": self.option("cog_identity"),
            "kegg_evalue": self.option("kegg_evalue"),
            "kegg_similarity": self.option("kegg_similarity"),
            "kegg_identity": self.option("kegg_identity"),
            "pfam_evalue": self.option("pfam_evalue"),
        }
        if self.option("exclude_taxon"):
            params.update({
                "exclude_taxon": self.option("exclude_taxon")
            })
        else:
            params.update({
                "exclude_taxon": ""
            })
        self.test_api.anno_type = 'latest'
        self.test_api.run_webroot(self.annotation.output_dir, self.option("trans2gene").prop['path'],
                                  params,
                                  task_id=self.option("task_id"),
                                  stat_id=self.option("stat_id"),
                                  last_id=self.option("last_id"),
                                  taxonomy=self.option("taxonomy"),
                                  version="v2",
                                  gene_exp=self.option("gene_exp").prop['path'],
                                  trans_exp=self.option("trans_exp").prop['path'],
        )
        self.end()

    def update_task_id(self, stat_id, blast_id, nr_id, swissprot_id):
        """更新主表task_id"""
        self.logger.info("更新主表task_id")
        db = Config().get_mongo_client(mtype="denovo_rna_v2")[Config().get_mongo_dbname("denovo_rna_v2")]
        #client = Config().mongo_client
        #db_name = Config().MONGODB + '_ref_rna'
        stat_coll = db['sg_annotation_stat']
        results = stat_coll.find_one({'_id': ObjectId(stat_id)})
        task_id = results['task_id']
        blast_coll = db['sg_annotation_blast']
        nr_coll = db['sg_annotation_nr']
        sw_coll = db['sg_annotation_swissprot']
        blast_coll.update({'_id': ObjectId(blast_id)}, {'$set': {'task_id': task_id}})
        nr_coll.update({'_id': ObjectId(nr_id)}, {'$set': {'task_id': task_id}})
        sw_coll.update({'_id': ObjectId(swissprot_id)}, {'$set': {'task_id': task_id}})
        self.logger.info("更新主表ID成功")

    def chart(self):
        chart = Chart()
        chart.work_dir = self.work_dir + "/"
        annot_stat = self.annotation.output_dir + '/all_stat.xls'
        venn_dir = os.path.dirname(annot_stat)
        annot_nr_species_pie_gene = os.path.join(self.annotation.output_dir, 'nr', "gene_nr_species_stat.xls")
        annot_nr_evalue_pie_gene = os.path.join(self.annotation.output_dir, 'nr', "gene_nr_evalue.xls")
        annot_nr_similar_pie_gene = os.path.join(self.annotation.output_dir, 'nr', "gene_nr_similar.xls")
        annot_nr_species_pie_trans = os.path.join(self.annotation.output_dir, 'nr', "tran_nr_species_stat.xls")
        annot_nr_evalue_pie_trans = os.path.join(self.annotation.output_dir, 'nr', "trans_nr_evalue.xls")
        annot_nr_similar_pie_trans = os.path.join(self.annotation.output_dir, 'nr', "trans_nr_similar.xls")
        print annot_nr_species_pie_gene
        print annot_nr_evalue_pie_gene
        print annot_nr_similar_pie_gene
        print annot_nr_species_pie_trans
        print annot_nr_evalue_pie_trans
        print annot_nr_similar_pie_trans
        annot_swissprot_evalue_pie_gene = os.path.join(self.annotation.output_dir, 'swissprot', "gene_swissprot_evalue.xls")
        annot_swissprot_similar_pie_gene = os.path.join(self.annotation.output_dir, 'swissprot', "gene_swissprot_similar.xls")
        annot_swissprot_evalue_pie_trans = os.path.join(self.annotation.output_dir, 'swissprot', "trans_swissprot_evalue.xls")
        annot_swissprot_similar_pie_trans = os.path.join(self.annotation.output_dir, 'swissprot', "trans_swissprot_similar.xls")
        # annotation pfam bar
        annot_pfam_bar_gene = os.path.join(self.annotation.output_dir, 'pfam', "pfam_domain_gene.xls")
        annot_pfam_bar_trans = os.path.join(self.annotation.output_dir, 'pfam', "pfam_domain_tran.xls")
        # annotation cog bar
        annot_cog_bar_gene = os.path.join(self.annotation.output_dir, 'cog', "summary.G.tsv")
        annot_cog_bar_trans = os.path.join(self.annotation.output_dir, 'cog', "summary.T.tsv")
        # annotation go pie level
        annot_go2_pie_gene = os.path.join(self.annotation.output_dir, 'go', "go_lev2_gene.stat.xls")
        annot_go3_pie_gene = os.path.join(self.annotation.output_dir, 'go', "go_lev3_gene.stat.xls")
        annot_go4_pie_gene = os.path.join(self.annotation.output_dir, 'go', "go_lev4_gene.stat.xls")
        annot_go2_pie_trans = os.path.join(self.annotation.output_dir, 'go', "go_lev2_tran.stat.xls")
        annot_go3_pie_trans = os.path.join(self.annotation.output_dir, 'go', "go_lev3_tran.stat.xls")
        annot_go4_pie_trans = os.path.join(self.annotation.output_dir, 'go', 'go_lev4_tran.stat.xls')
        # annotation kegg layer
        annot_kegg_layer_gene = os.path.join(self.annotation.output_dir, 'kegg', 'kegg_layer_gene.xls')
        annot_kegg_layer_trans = os.path.join(self.annotation.output_dir, 'kegg', 'kegg_layer_tran.xls')
        chart.denovo_chart_annotation_stat("", self.option("gene_exp").prop['path'], self.option("trans_exp").prop['path'], annot_stat, venn_dir)
        chart.denovo_chart_annotation_nr_species(annot_nr_species_pie_gene, 'gene')
        chart.denovo_chart_annotation_nr_evalue(annot_nr_evalue_pie_gene, 'gene')
        chart.denovo_chart_annotation_nr_similary(annot_nr_similar_pie_gene, 'gene')
        chart.denovo_chart_annotation_nr_species(annot_nr_species_pie_trans, 'tran')
        chart.denovo_chart_annotation_nr_evalue(annot_nr_evalue_pie_trans, 'tran')
        chart.denovo_chart_annotation_nr_similary(annot_nr_similar_pie_trans, 'tran')
        chart.denovo_chart_annotation_swissprot_evalue(annot_swissprot_evalue_pie_gene, 'gene')
        chart.denovo_chart_annotation_swissprot_similary(annot_swissprot_similar_pie_gene, 'gene')
        chart.denovo_chart_annotation_swissprot_evalue(annot_swissprot_evalue_pie_trans, 'tran')
        chart.denovo_chart_annotation_swissprot_similary(annot_swissprot_similar_pie_trans, 'tran')
        chart.denovo_chart_annotation_pfam(annot_pfam_bar_gene, 'gene')
        chart.denovo_chart_annotation_pfam(annot_pfam_bar_trans, 'tran')
        chart.denovo_chart_annotation_cog(annot_cog_bar_gene, 'gene')
        chart.denovo_chart_annotation_cog(annot_cog_bar_trans, 'tran')
        chart.denovo_chart_annotation_go(annot_go2_pie_gene, 2, 'gene')
        chart.denovo_chart_annotation_go(annot_go3_pie_gene, 3, 'gene')
        chart.denovo_chart_annotation_go(annot_go4_pie_gene, 4, 'gene')
        chart.denovo_chart_annotation_go(annot_go2_pie_trans, 2, 'tran')
        chart.denovo_chart_annotation_go(annot_go3_pie_trans, 3, 'tran')
        chart.denovo_chart_annotation_go(annot_go4_pie_trans, 4, 'tran')
        chart.denovo_chart_annotation_kegg(annot_kegg_layer_gene, 'gene')
        chart.denovo_chart_annotation_kegg(annot_kegg_layer_trans, 'tran')
        chart.to_pdf()

        # move pdf to result dir
        pdf_file = glob.glob(self.work_dir + '/*annot*.pdf')
        for i in pdf_file:
            os.link(i, self.annotation.output_dir + '/' + os.path.basename(i))

        pdf_file = glob.glob(self.work_dir + "/*nr*.pdf")
        for i in pdf_file:
            os.link(i, self.annotation.output_dir + "/nr/" + os.path.basename(i))

        pdf_file = glob.glob(self.work_dir + "/*swissprot*.pdf")
        for i in pdf_file:
            os.link(i, self.annotation.output_dir + "/swissprot/" + os.path.basename(i))

        pdf_file = glob.glob(self.work_dir + "/*Pfam*.pdf")
        for i in pdf_file:
            os.link(i, self.annotation.output_dir + "/pfam/" + os.path.basename(i))

        pdf_file = glob.glob(self.work_dir + "/*COG*.pdf")
        for i in pdf_file:
            os.link(i, self.annotation.output_dir + "/cog/" + os.path.basename(i))

        pdf_file = glob.glob(self.work_dir + "/*GO*.pdf")
        for i in pdf_file:
            os.link(i, self.annotation.output_dir + "/go/" + os.path.basename(i))

        pdf_file = glob.glob(self.work_dir + "/*KEGG*.pdf")
        for i in pdf_file:
            os.link(i, self.annotation.output_dir + "/kegg/" + os.path.basename(i))


        # pdf_file = glob.glob(self.work_dir + "/Trans*nr*.pdf")
        # for i in pdf_file:
        #     os.link(i, self.annotation.output_dir + "/Transcript_Anno/nr/" + os.path.basename(i))
        # pdf_file = glob.glob(self.work_dir + "/Trans*swissprot*.pdf")
        # for i in pdf_file:
        #     os.link(i, self.annotation.output_dir + "/Transcript_Anno/swissprot/" + os.path.basename(i))
        # pdf_file = glob.glob(self.work_dir + "/Trans*Pfam*.pdf")
        # for i in pdf_file:
        #     os.link(i, self.annotation.output_dir + "/Transcript_Anno/pfam/" + os.path.basename(i))
        # pdf_file = glob.glob(self.work_dir + "/Trans*COG*.pdf")
        # for i in pdf_file:
        #     os.link(i, self.annotation.output_dir + "/Transcript_Anno/cog/" + os.path.basename(i))
        # pdf_file = glob.glob(self.work_dir + "/Trans*GO*.pdf")
        # for i in pdf_file:
        #     os.link(i, self.annotation.output_dir + "/Transcript_Anno/go/" + os.path.basename(i))
        # pdf_file = glob.glob(self.work_dir + "/Trans*KEGG*.pdf")
        # for i in pdf_file:
        #     os.link(i, self.annotation.output_dir + "/Transcript_Anno/kegg/" + os.path.basename(i))



    def end(self):
        self.chart()
        origin_dir = self.annotation.output_dir
        target_dir = self.output_dir
        
        os.system('cp {0}/all* {1}'.format(origin_dir,target_dir))
        os.system('cp {0}/annot*.pdf {1}'.format(origin_dir, target_dir))
        for d in ['cog','go','kegg','nr','pfam','swissprot']:
            dirct = '{}/{}/'.format(origin_dir,d)
            for f in os.listdir(dirct):
                if os.path.isfile(os.path.join(dirct,f)):
                    if '.T.' in str(f) or 'tran' in str(f):
                        os.system('install -D {0}/{1}/{2} {3}/Transcript_Anno/{1}/{2}'.format(origin_dir,d,f,target_dir))
                    elif '.G.' in str(f) or 'gene' in str(f):
                        os.system('install -D {0}/{1}/{2} {3}/Unigene_Anno/{1}/{2}'.format(origin_dir,d,f,target_dir))
                elif os.path.isdir(os.path.join(dirct,f)):
                    if '.T.' in str(f) or 'tran' in str(f):
                        if not os.path.exists("{1}/Transcript_Anno/{0}".format(d,target_dir)):
                            os.makedirs( "{1}/Transcript_Anno/{0}".format(d,target_dir))
                        os.system('cp -r {0}/{1}/{2} {3}/Transcript_Anno/{1}/{2}'.format(origin_dir,d,f,target_dir))
                    elif '.G.' in str(f) or 'gene' in str(f):
                        if not os.path.exists("{1}/Unigene_Anno/{0}".format(d,target_dir)):
                            os.makedirs("{1}/Unigene_Anno/{0}".format(d,target_dir))
                        os.system('cp -r {0}/{1}/{2} {3}/Unigene_Anno/{1}/{2}'.format(origin_dir,d,f,target_dir))


        # CopyFile().linkdir(self.annotation.output_dir,
        #                    os.path.join(self.output_dir, 'Annotation'))
        # os.mkdir(target_dir + "/Annotation")
        '''
        os.mkdir(target_dir + "/Annotation/blast_xml")
        blast_nr = os.path.join(origin_dir + "/blast_xml/Trinity_vs_nr.xml")
        os.link(blast_nr, target_dir + "/Annotation/blast_xml/Trinity_vs_nr.xml")
        blast_swiss = os.path.join(origin_dir + "/blast_xml/Trinity_vs_swissprot.xml")
        os.link(blast_swiss, target_dir + "/Annotation/blast_xml/Trinity_vs_swissprot.xml")
        blast_string = os.path.join(origin_dir + "/blast_xml/Trinity_vs_string.xml")
        os.link(blast_string, target_dir + "/Annotation/blast_xml/Trinity_vs_string.xml")
        blast_kegg = os.path.join(origin_dir + "/blast_xml/Trinity_vs_kegg.xml")
        os.link(blast_kegg, target_dir + "/Annotation/blast_xml/Trinity_vs_kegg.xml")
        blast_pfam = os.path.join(origin_dir + "/blast_xml/pfam_domain")
        os.link(blast_pfam, target_dir + "/Annotation/blast_xml/pfam_domain")
        '''
        '''
        os.mkdir(target_dir + "/Annotation/AnnoStat")
        gene_anno_stat = os.path.join(origin_dir + "/anno_stat/all_annotation_statistics.xls")
        os.link(gene_anno_stat, target_dir + "/Annotation/AnnoStat/all_annotation_statistics.xls")
        os.mkdir(target_dir + "/Annotation/AnnoStat/NR")
        gene_nr = os.path.join(origin_dir + "/anno_stat/blast/gene_nr.xls")
        os.link(gene_nr, target_dir + "/Annotation/AnnoStat/NR/unigene_nr_anno_detail.xls")
        trans_nr = os.path.join(origin_dir + "/anno_stat/blast/nr.xls")
        os.link(trans_nr, target_dir + "/Annotation/AnnoStat/NR/transcript_nr_anno_detail.xls")
        taxon_nr = os.path.join(origin_dir + "/blast_nr_taxon/query_taxons_detail.xls.pdf")
        os.link(taxon_nr, target_dir + "/Annotation/AnnoStat/NR/query_taxons_detail.xls.pdf")
        os.mkdir(target_dir + "/Annotation/AnnoStat/Swiss-Prot")
        gene_swiss = os.path.join(origin_dir + "/anno_stat/blast/gene_swissprot.xls")
        os.link(gene_swiss, target_dir + "/Annotation/AnnoStat/Swiss-Prot/unigene_swissprot_anno_detail.xls")
        trans_swiss = os.path.join(origin_dir + "/anno_stat/blast/swissprot.xls")
        os.link(trans_swiss, target_dir + "/Annotation/AnnoStat/Swiss-Prot/transcript_swissprot_anno_detail.xls")
        os.mkdir(target_dir + "/Annotation/AnnoStat/Pfam")
        gene_pfam = os.path.join(origin_dir + "/anno_stat/pfam_stat/gene_pfam_domain")
        os.link(gene_pfam, target_dir + "/Annotation/AnnoStat/Pfam/unigene_pfam_anno_detail.xls")
        trans_pfam = os.path.join(origin_dir + "/blast_xml/pfam_domain")
        os.link(trans_pfam, target_dir + "/Annotation/AnnoStat/Pfam/transcript_pfam_anno_detail.xls")
        os.mkdir(target_dir + "/Annotation/AnnoStat/COG")
        gene_coglist = os.path.join(origin_dir + "/anno_stat/cog_stat/gene_cog_list.xls")
        os.link(gene_coglist, target_dir + "/Annotation/AnnoStat/COG/unigene_cog_list.xls")
        gene_coganno = os.path.join(origin_dir + "/anno_stat/cog_stat/gene_cog_table.xls")
        os.link(gene_coganno, target_dir + "/Annotation/AnnoStat/COG/unigene_cog_anno_detail.xls")
        gene_cogsum = os.path.join(origin_dir + "/anno_stat/cog_stat/gene_cog_summary.xls")
        os.link(gene_cogsum, target_dir + "/Annotation/AnnoStat/COG/unigene_cog_summary.xls")
        trans_coglist = os.path.join(origin_dir + "/cog/cog_list.xls")
        os.link(trans_coglist, target_dir + "/Annotation/AnnoStat/COG/transcript_cog_list.xls")
        trans_coganno = os.path.join(origin_dir + "/cog/cog_table.xls")
        os.link(trans_coganno, target_dir + "/Annotation/AnnoStat/COG/transcript_cog_anno_detail.xls")
        trans_cogsum = os.path.join(origin_dir + "/cog/cog_summary.xls")
        os.link(trans_cogsum, target_dir + "/Annotation/AnnoStat/COG/transcript_cog_summary.xls")
        os.mkdir(target_dir + "/Annotation/AnnoStat/GO")
        gene_golist = os.path.join(origin_dir + "/anno_stat/go_stat/gene_gos.list")
        os.link(gene_golist, target_dir + "/Annotation/AnnoStat/GO/unigene_gos.list")
        gene_goanno = os.path.join(origin_dir + "/anno_stat/go_stat/gene_blast2go.annot")
        os.link(gene_goanno, target_dir + "/Annotation/AnnoStat/GO/unigene_blast2go.annot")
        gene_golevel12 = os.path.join(origin_dir + "/anno_stat/go_stat/gene_go12level_statistics.xls")
        os.link(gene_golevel12, target_dir + "/Annotation/AnnoStat/GO/unigene_go12level_statistics.xls")
        gene_golevel123 = os.path.join(origin_dir + "/anno_stat/go_stat/gene_go123level_statistics.xls")
        os.link(gene_golevel123, target_dir + "/Annotation/AnnoStat/GO/unigene_go123level_statistics.xls")
        gene_golevel1234 = os.path.join(origin_dir + "/anno_stat/go_stat/gene_go1234level_statistics.xls")
        os.link(gene_golevel1234, target_dir + "/Annotation/AnnoStat/GO/unigene_go1234level_statistics.xls")
        trans_golist = os.path.join(origin_dir + "/go/query_gos.list")
        os.link(trans_golist, target_dir + "/Annotation/AnnoStat/GO/transcript_gos.list")
        trans_goanno = os.path.join(origin_dir + "/go/blast2go.annot")
        os.link(trans_goanno, target_dir + "/Annotation/AnnoStat/GO/transcript_blast2go.annot")
        trans_golevel12 = os.path.join(origin_dir + "/go/go12level_statistics.xls")
        os.link(trans_golevel12, target_dir + "/Annotation/AnnoStat/GO/transcript_go12level_statistics.xls")
        trans_golevel123 = os.path.join(origin_dir + "/go/go123level_statistics.xls")
        os.link(trans_golevel123, target_dir + "/Annotation/AnnoStat/GO/transcript_go123level_statistics.xls")
        trans_golevel1234 = os.path.join(origin_dir + "/go/go1234level_statistics.xls")
        os.link(trans_golevel1234, target_dir + "/Annotation/AnnoStat/GO/transcript_go1234level_statistics.xls")
        os.mkdir(target_dir + "/Annotation/AnnoStat/KEGG")
        os.mkdir(target_dir + "/Annotation/AnnoStat/KEGG/unigene_pathways")
        os.mkdir(target_dir + "/Annotation/AnnoStat/KEGG/transcript_pathways")
        gene_kegg_table = os.path.join(origin_dir + "/anno_stat/kegg_stat/gene_kegg_table.xls")
        os.link(gene_kegg_table, target_dir + "/Annotation/AnnoStat/KEGG/unigene_kegg_table.xls")
        gene_pathway_table = os.path.join(origin_dir + "/anno_stat/kegg_stat/gene_pathway_table.xls")
        os.link(gene_pathway_table, target_dir + "/Annotation/AnnoStat/KEGG/unigene_pathway_table.xls")
        gene_kegg_layer = os.path.join(origin_dir + "/anno_stat/kegg_stat/gene_kegg_layer.xls")
        os.link(gene_kegg_layer, target_dir + "/Annotation/AnnoStat/KEGG/unigene_kegg_layer.xls")
        trans_kegg_table = os.path.join(origin_dir + "/kegg/kegg_table.xls")
        os.link(trans_kegg_table, target_dir + "/Annotation/AnnoStat/KEGG/transcript_kegg_table.xls")
        trans_pathway_table = os.path.join(origin_dir + "/kegg/pathway_table.xls")
        os.link(trans_pathway_table, target_dir + "/Annotation/AnnoStat/KEGG/transcript_pathway_table.xls")
        trans_kegg_layer = os.path.join(origin_dir + "/kegg/kegg_layer.xls")
        os.link(trans_kegg_layer, target_dir + "/Annotation/AnnoStat/KEGG/transcript_kegg_layer.xls")
        for file in os.listdir(origin_dir + "/anno_stat/kegg_stat/gene_pathway"):
            gene_pathway_path = os.path.join(origin_dir + "/anno_stat/kegg_stat/gene_pathway", file)
            os.link(gene_pathway_path, target_dir + "/Annotation/AnnoStat/KEGG/unigene_pathways/" + file)
        for file in os.listdir(origin_dir + "/kegg/pathways"):
            trans_pathway_path = os.path.join(origin_dir + "/kegg/pathways", file)
            os.link(trans_pathway_path, target_dir + "/Annotation/AnnoStat/KEGG/transcript_pathways/" + file)
        # AnnoQuery
        os.mkdir(target_dir + "/AnnoQuery")
        gene_anno_detail = os.path.join(origin_dir + "/anno_stat/gene_anno_detail.xls")
        os.link(gene_anno_detail, target_dir + "/AnnoQuery/unigene_anno_detail.xls")
        trans_anno_detail = os.path.join(origin_dir + "/anno_stat/trans_anno_detail.xls")
        os.link(trans_anno_detail, target_dir + "/AnnoQuery/transcript_anno_detail.xls")
        '''

        repaths = [
                    [".", "", "转录组功能注释文件",0],
                    ["all_annot.xls", "", " ", 0 ],
                    ["all_stat.xls", "", "注释概况统计表", 0 ],
                    ["all_stat_detail.xls", "", " ", 0 ],
                    ["all_tran2gene.txt", "", " ", 0 ],
                    ["Transcript_Anno", "", "Transcript注释结果", 0 ],
                    ["Transcript_Anno/cog", "", "COG注释结果", 0 ],
                    ["Transcript_Anno/cog/cog_list_tran.xls", "", " ", 0 ],
                    ["Transcript_Anno/cog/cog_venn_tran.txt", "", " ", 0 ],
                    ["Transcript_Anno/cog/summary.T.tsv", "", "COG分类统计表", 0 ],
                    ["Transcript_Anno/go", "", "GO分类统计结果", 0 ],
                    ["Transcript_Anno/go/go_lev2_tran.stat.xls", "", "GO二级分类统计表", 0 ],
                    ["Transcript_Anno/go/go_lev3_tran.stat.xls", "", "GO三级分类统计表", 0 ],
                    ["Transcript_Anno/go/go_lev4_tran.stat.xls", "", "GO四级分类统计表", 0 ],
                    ["Transcript_Anno/go/go_list_tran.xls", "", "序列对应GO分类列表", 0 ],
                    ["Transcript_Anno/go/go_venn_tran.txt", "", " ", 0 ],
                    ["Transcript_Anno/kegg", "", "KEGG注释结果", 0 ],
                    ["Transcript_Anno/kegg/kegg_gene_tran.xls", "", "Pathway分类统计表", 0 ],
                    ["Transcript_Anno/kegg/kegg_layer_tran.xls", "", " ", 0 ],
                    ["Transcript_Anno/kegg/kegg_pathway_tran.xls", "", "Pathway对应的序列统计表", 0 ],
                    ["Transcript_Anno/kegg/kegg_venn_tran.txt", "", " ", 0 ],
                    ["Transcript_Anno/kegg/kegg_pathway_tran_dir", "", "KEGG分析结果通路图", 0 ],
                    ["Transcript_Anno/kegg/kegg_pathway_tran_dir/*.html", "", "KEGG通路html文件", 0 ],
                    ["Transcript_Anno/kegg/kegg_pathway_tran_dir/*.html.mark", "", " ", 0 ],
                    ["Transcript_Anno/kegg/kegg_pathway_tran_dir/*.png", "", "KEGG通路png图片", 0 ],
                    ["Transcript_Anno/nr", "", "NR注释结果", 0 ],
                    ["Transcript_Anno/nr/nr.T.tsvspecies.txt", "", " ", 0 ],
                    ["Transcript_Anno/nr/nr.T.tsvspecies.txt.exp.xls", "", " ", 0 ],
                    ["Transcript_Anno/nr/nr.T.tsvspecies.txt.pdf", "", " ", 0 ],
                    ["Transcript_Anno/nr/nr_blast_tran.xls", "", " ", 0 ],
                    ["Transcript_Anno/nr/nr_venn_tran.txt", "", " ", 0 ],
                    ["Transcript_Anno/nr/tran_nr_species_stat.xls", "", " ", 0 ],
                    ["Transcript_Anno/nr/trans_nr_evalue.xls", "", " ", 0 ],
                    ["Transcript_Anno/nr/trans_nr_similar.xls", "", " ", 0 ],
                    ["Transcript_Anno/pfam", "", "Pfam注释结果", 0 ],
                    ["Transcript_Anno/pfam/pfam_domain_tran.xls", "", "Pfam注释详情表", 0 ],
                    ["Transcript_Anno/pfam/pfam_venn_tran.txt", "", " ", 0 ],
                    ["Transcript_Anno/swissprot", "", "Swiss-Prot注释结果", 0 ],
                    ["Transcript_Anno/swissprot/swissprot_blast_tran.xls", "", " ", 0 ],
                    ["Transcript_Anno/swissprot/swissprot_venn_tran.txt", "", " ", 0 ],
                    ["Transcript_Anno/swissprot/trans_swissprot_evalue.xls", "", " ", 0 ],
                    ["Transcript_Anno/swissprot/trans_swissprot_similar.xls", "", " ", 0 ],
                    ["Unigene_Anno", "", "Unigene注释结果", 0 ],
                    ["Unigene_Anno/cog", "", "COG注释结果", 0 ],
                    ["Unigene_Anno/cog/cog_venn_gene.txt", "", " ", 0 ],
                    ["Unigene_Anno/cog/summary.G.tsv", "", "COG分类统计表", 0 ],
                    ["Unigene_Anno/go", "", "GO分类统计结果", 0 ],
                    ["Unigene_Anno/go/go_lev2_gene.stat.xls", "", "GO二级分类统计表", 0 ],
                    ["Unigene_Anno/go/go_lev3_gene.stat.xls", "", "GO三级分类统计表", 0 ],
                    ["Unigene_Anno/go/go_lev4_gene.stat.xls", "", "GO四级分类统计表", 0 ],
                    ["Unigene_Anno/go/go_list_gene.xls", "", "序列对应GO分类列表", 0 ],
                    ["Unigene_Anno/go/go_venn_gene.txt", "", " ", 0 ],
                    ["Unigene_Anno/kegg", "", "KEGG注释结果", 0 ],
                    ["Unigene_Anno/kegg/kegg_gene_gene.xls", "", "Pathway分类统计表", 0 ],
                    ["Unigene_Anno/kegg/kegg_layer_gene.xls", "", " ", 0 ],
                    ["Unigene_Anno/kegg/kegg_pathway_gene.xls", "", "Pathway对应的序列统计表", 0 ],
                    ["Unigene_Anno/kegg/kegg_venn_gene.txt", "", " ", 0 ],
                    ["Unigene_Anno/kegg/kegg_pathway_gene_dir", "", "KEGG分析结果通路图", 0 ],
                    ["Unigene_Anno/kegg/kegg_pathway_gene_dir/*.html", "", "KEGG通路html文件", 0 ],
                    ["Unigene_Anno/kegg/kegg_pathway_gene_dir/*.html.mark", "", " ", 0 ],
                    ["Unigene_Anno/kegg/kegg_pathway_gene_dir/*.png", "", "KEGG通路png图片", 0 ],
                    ["Unigene_Anno/nr", "", "NR注释结果", 0 ],
                    ["Unigene_Anno/nr/gene_nr_evalue.xls", "", " ", 0 ],
                    ["Unigene_Anno/nr/gene_nr_similar.xls", "", " ", 0 ],
                    ["Unigene_Anno/nr/gene_nr_species_stat.xls", "", " ", 0 ],
                    ["Unigene_Anno/nr/nr.G.tsvspecies.txt", "", " ", 0 ],
                    ["Unigene_Anno/nr/nr.G.tsvspecies.txt.exp.xls", "", " ", 0 ],
                    ["Unigene_Anno/nr/nr.G.tsvspecies.txt.pdf", "", " ", 0 ],
                    ["Unigene_Anno/nr/nr_blast_gene.xls", "", " ", 0 ],
                    ["Unigene_Anno/nr/nr_venn_gene.txt", "", " ", 0 ],
                    ["Unigene_Anno/pfam", "", "Pfam注释结果", 0 ],
                    ["Unigene_Anno/pfam/pfam_domain_gene.xls", "", "Pfam注释详情表", 0 ],
                    ["Unigene_Anno/pfam/pfam_venn_gene.txt", "", " ", 0 ],
                    ["Unigene_Anno/swissprot", "", "Swiss-Prot注释结果", 0 ],
                    ["Unigene_Anno/swissprot/gene_swissprot_evalue.xls", "", " ", 0 ],
                    ["Unigene_Anno/swissprot/gene_swissprot_similar.xls", "", " ", 0 ],
                    ["Unigene_Anno/swissprot/swissprot_blast_gene.xls", "", "", 0 ],
                    ["Unigene_Anno/swissprot/swissprot_venn_gene.txt", "", " ", 0 ],
                    ['run_parameter.txt', 'txt', '运行参数日志', 0],
        ]

        self.workflow_output_tmp = self._sheet.output
        if re.match(r'tsanger:', self.workflow_output_tmp):
            self.workflow_output = self.workflow_output_tmp.replace('tsanger:', '/mnt/ilustre/tsanger-data/')
        else:
            self.workflow_output = self.workflow_output_tmp.replace('sanger:', '/mnt/ilustre/data/')

        if os.path.exists(os.path.join(target_dir, os.path.basename(self.run_log))):
            os.remove(os.path.join(target_dir, os.path.basename(self.run_log)))
        os.link(self.run_log, os.path.join(target_dir, os.path.basename(self.run_log)))

        result_dir = self.add_upload_dir(target_dir)
        self.inter_dirs = [
            ["01 Annotation", "", "转录组功能注释结果目录",0],
        ]
        result_dir.add_regexp_rules(repaths)
        db = Config().get_mongo_client(mtype="denovo_rna_v2")[Config().get_mongo_dbname("denovo_rna_v2")]
        col1 = db["sg_annotation_stat"]
        col1.update({"_id": ObjectId(self.option("stat_id"))}, {"$set": {"result_dir": self.workflow_output}}, upsert=True)

        super(DenovoAnnotationWorkflow, self).end()

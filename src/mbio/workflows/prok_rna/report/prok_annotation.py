# -*- coding: utf-8 -*-
# __author__ = 'liubinxu'
# last_modifiy = modified 2018.03.27

from biocluster.workflow import Workflow
from biocluster.config import Config
from bson.son import SON
from bson.objectid import ObjectId
from biocluster.wsheet import Sheet
import os
import re
import shutil
import json
from biocluster.config import Config
from mbio.packages.rna.annot_config import AnnotConfig
from mbio.packages.prok_rna.chart import Chart
import tarfile


class ProkAnnotationWorkflow(Workflow):
    """
    交互分析进行筛选blast参数重注释时使用
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        #self.rpc = False
        super(ProkAnnotationWorkflow, self).__init__(wsheet_object)
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
            {"name": "kegg_species", "type": "string", "default": ""},
            {"name": "pfam_evalue", "type": "float", "default": 1e-5},
            {"name": "stat_id", "type": "string"},
            {"name": "task_id", "type": "string"},
            {"name": "origin_result", "type": "infile", "format": "prok_rna.common_dir"},
            {"name": "last_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "trans2gene", "type": "string"},
            {"name": "taxonomy", "type": "string", "default": None},
            {"name": "origin_param", "type": "string", "default": None},
            {"name": "annot_group", "type": "string", "default": "REFRNA_GROUP_202007"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.ref_filter = self.add_module("ref_rna_v2.annot_filter")

        self.annotation = self.add_module("prok_rna.annot_class")
        self.annot_config_dict = AnnotConfig().get_group_option_detail(section=self.option("annot_group"))
        self.old_param = json.loads(self.option("origin_param"))
        self.step.add_steps("ref_filter", "ref_class")

    def set_step(self, event):
        if 'start' in event['data'].keys():
            event['data']['start'].start()
        if 'end' in event['data'].keys():
            event['data']['end'].finish()
        self.step.update()

    def run_ref_filter(self):
        '''
        根据筛选参数过滤参考注释
        '''

        if self.option("origin_result").prop['path'].endswith("origin_result/"):
            ref_annot_dir = self.option("origin_result").prop['path'] + "Annotation/anvermapdb"
        else:
            ref_annot_dir = self.option("origin_result").prop['path'] + "/annot_mapdb"
        if self.option("annot_group") in ["REFRNA_GROUP_202110"]:
            cog_xml = ref_annot_dir + "/cog/blast.xml"
        else:
            cog_xml = ref_annot_dir + "/eggnog/blast.xml"
        options = {
            "blast_nr_xml" : ref_annot_dir + "/nr/blast.xml",
            "blast_eggnog_xml" : cog_xml,
            "blast_kegg_xml": ref_annot_dir + "/kegg/blast.xml",
            "blast_swissprot_xml" : ref_annot_dir + "/swissprot/blast.xml",
            "blast2go_annot" : ref_annot_dir + "/GO/go_annot.xls",
            "pfam_domain" : ref_annot_dir + "/pfam/pfam_domain",
            'nr_evalue': self.option('nr_evalue'),
            'nr_identity': self.option('nr_identity'),
            'nr_similarity': self.option('nr_similarity'),
            'swissprot_evalue': self.option('swissprot_evalue'),
            'swissprot_identity': self.option('swissprot_identity'),
            'swissprot_similarity': self.option('swissprot_similarity'),
            'eggnog_evalue': self.option('cog_evalue'),
            'eggnog_identity': self.option('cog_identity'),
            'eggnog_similarity': self.option('cog_similarity'),
            'kegg_evalue': self.option('kegg_evalue'),
            'kegg_identity': self.option('kegg_identity'),
            'kegg_similarity': self.option('kegg_similarity'),
            'pfam_evalue': self.option('pfam_evalue'),
        }
        self.ref_filter.set_options(options)
        self.ref_filter.on('start', self.set_step, {'start': self.step.ref_filter})
        self.ref_filter.on('end', self.set_step, {'end': self.step.ref_filter})
        self.ref_filter.on('end', self.run_ref_class)
        self.ref_filter.run()

    def run_ref_class(self):
        ref_filter_dir = self.ref_filter.output_dir
        options = {
            'taxonomy': self.option('taxonomy'),
            'blast_nr_xml': ref_filter_dir + "/nr/blast.xml.filter.xml",
            'blast_kegg_xml' :ref_filter_dir + "/kegg/blast.xml.filter.xml",
            'blast_swissprot_xml': ref_filter_dir + "/swissprot/blast.xml.filter.xml",
            'pfam_domain': ref_filter_dir + "/pfam/pfam_domain.filter.xls",
            "blast2go_annot" : ref_filter_dir + "/go/go_annot.xls.filter.xls",
            "des" : self.option("origin_result").prop['path'] + "/ref.txt",
            "kegg_version" : self.annot_config_dict['kegg']['version'],
            "nr_version" : self.annot_config_dict['nr']['version'],
            "eggnog_version" : self.annot_config_dict['eggnog']['version'],
            "string_version" : self.annot_config_dict['string']['version'],
            "go_version" : self.annot_config_dict['go']['version'],
            "pir_version" : self.annot_config_dict['pir']['version'],
            "swissprot_version" : self.annot_config_dict['swissprot']['version'],

        }
        if self.option("annot_group") in ["REFRNA_GROUP_202110"]:
            options.update({
                'blast_ncbicog_xml': ref_filter_dir + "/eggnog/blast.xml.filter.xml",
                'cog_verion': self.annot_config_dict['cog']['version']
            })
        else:
            options.update({
                'blast_eggnog_xml': ref_filter_dir + "/eggnog/blast.xml.filter.xml"
            })
        self.annotation.set_options(options)
        self.annotation.on('end', self.set_output)
        self.annotation.run()

    def run(self):

        self.run_ref_filter()

        super(ProkAnnotationWorkflow, self).run()

    def chart(self):
        chart = Chart()
        chart.work_dir = self.work_dir + '/'
        annot_venn = os.path.join(self.annotation.output_dir, 'anno_stat/venn')
        annot_stats = os.path.join(self.annotation.output_dir, 'anno_stat/all_annotation_statistics.xls')
        cog_annot = os.path.join(self.annotation.output_dir, 'cog/cog_summary.xls')
        go_annot_level2 = os.path.join(self.annotation.output_dir, 'go/go_detail_stat.tsv')
        kegg_annot = os.path.join(self.annotation.output_dir, 'kegg/kegg_layer.xls')
        if os.path.exists(annot_venn):
            chart.prok_annot_venn(annot_venn)
        if os.path.exists(annot_stats):
            chart.prok_annot_bar(annot_stats)
        if os.path.exists(cog_annot):
            chart.prok_annot_cog_bar(cog_annot)
        if os.path.exists(go_annot_level2):
            chart.prok_annot_go_bar_pie(go_annot_level2)
        if os.path.exists(kegg_annot):
            chart.prok_annot_kegg_bar(kegg_file=kegg_annot)
        chart.to_pdf()

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
            self.set_error("需要移动到output目录的文件夹不存在", code = "15000601")
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
        self.test_api = self.api.api("prok_rna.prokrna_annotation")
        new_params = self.old_param
        new_params.update({
            "nr_evalue": str(self.option("nr_evalue")),
            "cog_evalue": str(self.option("cog_evalue")),
            "kegg_evalue": str(self.option("kegg_evalue")),
            "pfam_evalue": str(self.option("pfam_evalue")),
            "swissprot_evalue": str(self.option("swissprot_evalue")),
            "task_type": 2
        })
        self.test_api.anno_type = 'latest'
        genome_stat = self.option("origin_result").prop['path'] + "/genome_stat.xls"
        self.test_api.run_prok_webroot(genome_stat, self.annotation.output_dir, new_params, self.option("task_id"), self.option("stat_id"), self.option("last_id"))

        pathway_file = os.path.join(self.annotation.output_dir, 'kegg', 'pathways')
        pngs = os.listdir(pathway_file)
        tar_file = pathway_file + ".tar.gz"
        with tarfile.open(tar_file, mode='w:gz') as f:
            for png in pngs:
                f.add(pathway_file + "/" + png, arcname=png)
        shutil.rmtree(pathway_file)
        self.end()

    def move_file(self, old_file, new_file):
        if os.path.isfile(old_file):
            os.link(old_file, new_file)
        else:
            os.mkdir(new_file)
            for file in os.listdir(old_file):
                file_path = os.path.join(old_file, file)
                new_path = os.path.join(new_file, file)
                self.move_file(file_path, new_path)

    def move2outputdir(self, olddir, newname, mode='link'):
        """
        移动一个目录下的所有文件/文件夹到workflow输出文件夹下
        """
        # start = time.time()
        if not os.path.isdir(olddir):
            self.set_error("需要移动到output目录的文件夹不存在", code = "15000602")
        newdir = os.path.join(self.output_dir, newname)
        if not os.path.exists(newdir):
            os.makedirs(newdir)
        allfiles = os.listdir(olddir)
        oldfiles = [os.path.join(olddir, i) for i in allfiles]
        newfiles = [os.path.join(newdir, i) for i in allfiles]
        self.logger.info(newfiles)
        for newfile in newfiles:
            if os.path.isfile(newfile) and os.path.exists(newfile):
                os.remove(newfile)
            elif os.path.isdir(newfile) and os.path.exists(newfile):
                shutil.rmtree(newfile)
        for i in range(len(allfiles)):
            self.move_file(oldfiles[i], newfiles[i])
        # end = time.time()
        # duration = end - start
        self.logger.info("文件夹{}到{}".format(olddir, newdir))

    def move_pdf(self, origin, new):
        if os.path.exists(new):
            os.remove(new)
        if not os.path.exists(os.path.dirname(new)):
            os.mkdir(os.path.dirname(new))
        if os.path.exists(origin):
            os.link(origin, new)

    def end(self):
        self.chart()

        origin_dir = self.annotation.output_dir
        target_dir = self.output_dir
        self.move2outputdir(origin_dir, target_dir)

        # Annotation Chart
        cog_pdf = os.path.join(self.work_dir, 'annot_cog.cog_bar.pdf')
        new_path = os.path.join(target_dir, 'cog', 'cog_level.pdf')
        self.move_pdf(cog_pdf, new_path)

        go_bar_pdf = os.path.join(self.work_dir, 'annot_go_bar.multi_bar.pdf')
        new_path = os.path.join(target_dir, 'go', 'go_level_bar.pdf')
        self.move_pdf(go_bar_pdf, new_path)

        go_pie_pdf = os.path.join(self.work_dir, 'annot_go_pie.multi_pie.pdf')
        new_path = os.path.join(target_dir, 'go', 'go_level_pie.pdf')
        self.move_pdf(go_pie_pdf, new_path)

        kegg_pdf = os.path.join(self.work_dir, 'annot_kegg_bar.kegg_bar.pdf')
        new_path = os.path.join(target_dir, 'kegg', 'kegg_level.pdf')
        self.move_pdf(kegg_pdf, new_path)

        for i in [['stats_annot_venn.venn.pdf', 'all_anno_venn.pdf'],
                  ['stats_annot_num_bar.bar.pdf', 'all_anno_number.pdf'],
                  ['stats_annot_percent_bar.bar.pdf', 'all_anno_percent.pdf']]:
            pdf = os.path.join(self.work_dir, i[0])
            new_path = os.path.join(target_dir, 'summary', i[1])
            self.move_pdf(pdf, new_path)

        '''
        '''
        repaths = [
            [".", "", "基础注释重运行统计结果"],
            ["pfam_domain", "", "pfam注释详情表"],
            ["anno_stat", "", "注释汇总目录"],
            ["anno_stat/blast", "", "blast比对结果目录"],
            ["anno_stat/blast/swissprot.xls", "", "swissprot注释详情表"],
            ["anno_stat/blast/nr.xls", "", "nr注释详情表"],
            ["anno_stat/all_anno_detail.xls", "", "注释结果汇总详情表"],
            ["anno_stat/all_annotation_statistics.xls", "", "注释结果统计表"],
            ["anno_stat/venn", "", "基础注释统计Venn图结果目录"],
            ["anno_stat/venn/pfam_venn.txt", "", "存在pfam注释的基因列表"],
            ["anno_stat/venn/kegg_venn.txt", "", "存在kegg注释的基因列表"],
            ["anno_stat/venn/go_venn.txt", "", "存在go注释的基因列表"],
            ["anno_stat/venn/cog_venn.txt", "", "存在cog注释的基因列表"],
            ["anno_stat/venn/nr_venn.txt", "", "存在nr注释的基因列表"],
            ["anno_stat/venn/swissprot_venn.txt", "", "存在swissprot注释的基因列表"],
            ["go", "", "go注释结果目录"],
            ["go/GO.stat.xls", "", "go注释统计表"],
            ["go/query_gos.list", "", "go注释表"],
            ["go/GO.list", "", "go注释表"],
            ["go/go_detail.xls", "", "基因与go对应关系详情表"],
            ["go/go12level_statistics.xls", "", "go二级分类统计表"],
            ["go/go123level_statistics.xls", "", "鉴定CDS的GO1、2、3层级的详细信息"],
            ["go/go1234level_statistics.xls", "", "go四级分类统计表"],
            ["go/go_level_bar.pdf", "pdf", "GO注释分类统计柱形图"],
            ["go/go_level_pie.pdf", "pdf", "GO注释分类统计饼图"],
            ["cog", "", "cog注释结果目录"],
            ["cog/cog_stat.xls", "", "cog注释详情表"],
            ["cog/cog_table.xls", "", "鉴定CDS序列文件与COG比对结果表"],
            ["cog/cog_summary.xls", "", "cog注释统计表"],
            ["cog/cog_list.xls", "", "鉴定CDS与COG的对应关系"],
            ["cog/cog.xls", "", "基因与cog对应关系详情表"],
            ["cog/cog_level.pdf", "pdf", "COG分类统计柱状图"],
            ["kegg", "", "kegg注释结果目录"],
            ["kegg/pid.txt", "", "kegg注释pathway统计表"],
            ["kegg/kegg_layer.xls", "", "Pathway分类统计表"],
            ["kegg/kegg_table.xls", "", "基因与Pathway对应关系详情表"],
            ["kegg/pathways.tar.gz", "", "KEGG通路注释结果图片目录"],
            ["kegg/pathway_table.xls", "", "pathway信息详情表"],
            ["kegg/kegg_level.pdf", "", "pathway分类统计柱状图"],
            ["kegg/kegg_table.xls", "", "每行以CDS为单位显示通路注释详情表"],
            ["kegg/pid.txt", "", "每个注释通路和映射上去的CDS联系表"],
            ["summary", "", "注释结果汇总目录"],
            ["summary/all_anno_venn.pdf", "pdf", "基础注释统计Venn图"],
            ["summary/all_anno_percent.pdf", "pdf", "基础注释统计柱状占比图"],
            ["summary/all_anno_number.pdf", "pdf", "基础注释统计柱状图"],
            ["blast_xml", "", "blast比对结果目录", 1],
            ["blast_xml/vs_kegg.xml", "", "kegg数据库比对结果"],
            ["blast_xml/vs_eggnog.xml", "", "eggnog数据库比对结果式"],
            ["blast_xml/pfam_domain", "", "pfam数据库比对结果"],
            ["blast_xml/vs_nr.xml", "", "nr数据库比对结果"],
            ["blast_xml/vs_swissprot.xml", "", "swissprot数据库比对结果"],
            ["annot_mapdb", "", "数据库原始比对结果(evalue 1e-3)"]
        ]

        self.workflow_output_tmp = self._sheet.output
        if re.match(r'tsanger:', self.workflow_output_tmp):
            self.workflow_output = self.workflow_output_tmp.replace('tsanger:', '/mnt/ilustre/tsanger-data/')
        else:
            self.workflow_output = self.workflow_output_tmp.replace('sanger:', '/mnt/ilustre/data/')

        result_dir = self.add_upload_dir(target_dir)
        result_dir.add_regexp_rules(repaths)
        db = Config().get_mongo_client(mtype="prok_rna")[Config().get_mongo_dbname("prok_rna")]
        col1 = db["sg_annotation_stat"]
        col1.update({"_id": ObjectId(self.option("stat_id"))}, {"$set": {"result_dir": self.workflow_output}}, upsert=True)

        super(ProkAnnotationWorkflow, self).end()

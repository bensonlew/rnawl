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
from mbio.packages.dia_v3.chart import Chart
from biocluster.core.function import filter_error_info, link, CJsonEncoder
import glob


class ProteinAnnotationWorkflow(Workflow):
    """
    交互分析进行筛选blast参数重注释时使用
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        #self.rpc = False
        super(ProteinAnnotationWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "go_evalue", "type": "float", "default": 1e-5},
            {"name": "go_similarity", "type": "float", "default": 0},
            {"name": "go_identity", "type": "float", "default": 0},
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
            {"name": "origin_result", "type": "infile", "format": "labelfree.common_dir"},
            {"name": "last_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "trans2gene", "type": "string"},
            {"name": "taxonomy", "type": "string", "default": None},
            {"name": "origin_param", "type": "string", "default": None},
            {"name": "annot_group", "type": "string", "default": "REFRNA_GROUP_202007"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.nr_filter = self.add_tool("labelfree.annotation.filter_annot")
        self.swissprot_filter = self.add_tool("labelfree.annotation.filter_annot")
        self.cog_filter = self.add_tool("labelfree.annotation.filter_annot")
        self.eggnog_filter = self.add_tool("labelfree.annotation.filter_annot")
        self.kegg_filter = self.add_tool("labelfree.annotation.filter_annot")
        self.pfam_filter = self.add_tool("labelfree.annotation.filter_annot")
        self.annotation = self.add_module("labelfree.protein_annotation")
        self.old_param = json.loads(self.option("origin_param"))
        self.annot_config_dict = AnnotConfig().get_group_option_detail(section=self.option("annot_group"))
        self._sheet.output = self._sheet.output.replace('interaction_results', 'interaction_results/3_Anno')
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
        super(ProteinAnnotationWorkflow, self).send_log(data)

    def run_nr_blast_filter(self):
        options = {
            'xml': os.path.join(self.option("origin_result").prop['path'], "blast_xml/nr.xml"),
            'types': "xml",
            'evalue': self.option('go_evalue'),
            'identity': self.option('go_identity'),
            'similarity': self.option('go_similarity')
        }
        self.nr_filter.set_options(options)
        self.nr_filter.run()

    def run_swissprot_blast_filter(self):
        options = {
            'xml': os.path.join(self.option("origin_result").prop['path'], "blast_xml/swissprot.xml"),
            'types': "xml",
            'evalue': self.option('swissprot_evalue'),
            'identity': self.option('swissprot_identity'),
            'similarity': self.option('swissprot_similarity')
        }
        self.swissprot_filter.set_options(options)
        self.swissprot_filter.run()

    def run_cog_blast_filter(self):
        options = {
            'xml': os.path.join(self.option("origin_result").prop['path'], "blast_xml/string.xml"),
            'types': "xml",
            'evalue': self.option('cog_evalue'),
            'identity': self.option('cog_identity'),
            'similarity': self.option('cog_similarity')
        }
        self.cog_filter.set_options(options)
        self.cog_filter.run()

    def run_eggnog_blast_filter(self):
        options = {
            'xml': os.path.join(self.option("origin_result").prop['path'], "blast_xml/eggnog.xml"),
            'types': "xml",
            'evalue': self.option('cog_evalue'),
            'identity': self.option('cog_identity'),
            'similarity': self.option('cog_similarity')
        }
        self.eggnog_filter.set_options(options)
        self.eggnog_filter.run()

    def run_kegg_blast_filter(self):
        options = {
            'xml': os.path.join(self.option("origin_result").prop['path'], "blast_xml/kegg.xml"),
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

    def run_annotation(self):
        options = {
            'des': os.path.join(self.option("origin_result").prop['path'], "blast_xml/protein.xls"),
            'fa': os.path.join(self.option("origin_result").prop['path'], "blast_xml/protein.fa"),
            'go_annot': True,
            'blast_nr_xml': self.nr_filter.option('outxml').prop['path'],
            'nr_annot': True,
            'blast_kegg_xml': self.kegg_filter.option('outxml').prop['path'],
            'taxonomy': self.old_param['taxon'],
            'kegg_species': self.old_param['kegg_species'],
            'blast_string_xml': self.cog_filter.option('outxml').prop['path'],
            'pfam_domain': self.pfam_filter.option('outtable').prop['path'],
            "kegg_version" : self.annot_config_dict['kegg']['version'],
            "nr_version" : self.annot_config_dict['nr']['version'],
            "eggnog_version" : self.annot_config_dict['eggnog']['version'],
            "string_version" : self.annot_config_dict['string']['version'],
            "go_version" : self.annot_config_dict['go']['version'],
            "pir_version" : self.annot_config_dict['pir']['version'],
            "swissprot_version" : self.annot_config_dict['swissprot']['version'],
        }
        if self.option("annot_group") in ["REFRNA_GROUP_202007"]:
            options.update({
                'blast_eggnog_xml': self.cog_filter.option('outxml').prop['path']
            })
        self.annotation.set_options(options)
        self.annotation.on('end', self.set_output)
        self.annotation.run()

    def run(self):
        if self.option("annot_group") in ["REFRNA_GROUP_202007"]:
            self.on_rely([self.nr_filter, self.cog_filter, self.kegg_filter, self.pfam_filter, self.eggnog_filter], self.run_annotation)
            self.run_nr_blast_filter()
            self.run_cog_blast_filter()
            self.run_kegg_blast_filter()
            self.run_eggnog_blast_filter()
            # self.run_swissprot_blast_filter()
            self.run_pfam_filter()
        else:
            # 旧版
            self.on_rely([self.nr_filter, self.cog_filter, self.kegg_filter, self.pfam_filter], self.run_annotation)

            self.run_nr_blast_filter()
            self.run_cog_blast_filter()
            self.run_kegg_blast_filter()
            # self.run_swissprot_blast_filter()
            self.run_pfam_filter()

        # self.annotation.output_dir = "/mnt/ilustre/users/sanger-dev/workspace/20171224/ProteinAnnotation_test_anno_web_1224_113531/ProteinAnnotation/output"
        # self.set_output()
        super(ProteinAnnotationWorkflow, self).run()

    def chart(self):
        chart = Chart()
        chart.work_dir = self.work_dir + '/'
        chart.output_dir = self.work_dir + '/'
        all_annotation_stat = os.path.join(self.annotation.output_dir, "anno_stat", "all_annotation_statistics.xls")
        go12level_statistics = os.path.join(self.annotation.output_dir, "go", "go12level_statistics.xls")
        go123level_statistics = os.path.join(self.annotation.output_dir, "go", "go123level_statistics.xls")
        go1234level_statistics = os.path.join(self.annotation.output_dir, "go", "go1234level_statistics.xls")
        kegg_layer = os.path.join(self.annotation.output_dir, "kegg", "kegg_layer.xls")
        kegg_pathway = os.path.join(self.annotation.output_dir, "kegg", "pathway_table.xls")
        cog_summary = os.path.join(self.annotation.output_dir, "cog", "cog_summary.xls")
        pfam_domain = os.path.join(self.annotation.output_dir, "blast_xml", "pfam_domain")
        subloc_stat = os.path.join(self.annotation.output_dir, "subloc", "multiloc_stat.xls")
        if os.path.exists(all_annotation_stat):
            chart.chart_all_annotation_stat(all_annotation_stat)
        if os.path.exists(go12level_statistics):
            chart.chart_go12level_statistics(go12level_statistics)
        if os.path.exists(go123level_statistics):
            chart.chart_go123level_statistics(go123level_statistics)
        if os.path.exists(go1234level_statistics):
            chart.chart_go1234level_statistics(go1234level_statistics)
        if os.path.exists(kegg_layer):
            chart.chart_kegg_layer(kegg_layer)
        if os.path.exists(kegg_pathway):
            chart.chart_kegg_pathway(kegg_pathway)
        if os.path.exists(cog_summary):
            chart.chart_cog_summary(cog_summary)
        if os.path.exists(pfam_domain):
            chart.chart_pfam_domain(pfam_domain)
        if os.path.exists(subloc_stat):
            chart.chart_subloc_stat(subloc_stat)
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
            self.set_error('需要移动到output目录的文件夹不存在。', code = "12500801")
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
        self.test_api = self.api.api("dia.protein_annotation")
        new_params = self.old_param
        new_params.update({
            "go_evalue": str('%.0e'%(self.option("go_evalue"))),
            "cog_evalue": str('%.0e'%(self.option("cog_evalue"))),
            "kegg_evalue": str('%.0e'%(self.option("kegg_evalue"))),
            "pfam_evalue": str('%.0e'%(self.option("pfam_evalue"))),
            "go_identity": str(self.option("go_identity")),
            "cog_identity": str(self.option("cog_identity")),
            "kegg_identity": str(self.option("kegg_identity")),
            "task_type": 2
        })
        self.test_api.anno_type = 'latest'
        self.test_api.run_protein_webroot(self.annotation.output_dir, new_params, task_id=self.option("task_id"), stat_id=self.option("stat_id"), last_id=self.option("last_id"))
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
            self.set_error('需要移动到output目录的文件夹不存在。', code = "12500802")
        newdir = os.path.join(self.output_dir, newname)
        if not os.path.exists(newdir):
            os.makedirs(newdir)
        allfiles = os.listdir(olddir)

        newfiles = list()
        for each in allfiles:
            if each == 'anno_stat':
                newfiles.append(os.path.join(newdir, '01_anno_Stat'))
            elif each == 'go':
                newfiles.append(os.path.join(newdir, '02_anno_GO'))
            elif each == 'kegg':
                newfiles.append(os.path.join(newdir, '03_anno_KEGG'))
            elif each == 'cog':
                newfiles.append(os.path.join(newdir, '04_anno_COG'))
            elif each == 'pfam':
                newfiles.append(os.path.join(newdir, '05_anno_Pfam'))
            elif each == 'subloc':
                newfiles.append(os.path.join(newdir, '06_anno_SubLoc'))
            else:
                newfiles.append(os.path.join(newdir, each))

        oldfiles = [os.path.join(olddir, i) for i in allfiles]
        # newfiles = [os.path.join(newdir, i) for i in allfiles]
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

    def end(self):
        self.chart()

        origin_dir = self.annotation.output_dir
        target_dir = self.output_dir
        self.move2outputdir(origin_dir, target_dir)

        # move pdf to result dir
        def movefile(old, new):
            if not os.path.exists(os.path.dirname(new)):
                os.mkdir(os.path.dirname(new))
            if os.path.exists(new):
                os.remove(new)
            os.link(old, new)

        pdf_files = glob.glob(os.path.join(self.work_dir, 'go_lev*pdf'))
        for i in pdf_files:
            movefile(i, os.path.join(self.output_dir, '02_anno_GO', os.path.basename(i)))

        pdf_files = glob.glob(os.path.join(self.work_dir, '*path*pdf'))
        for i in pdf_files:
            movefile(i, os.path.join(self.output_dir, '03_anno_KEGG', os.path.basename(i)))

        pdf_files = glob.glob(os.path.join(self.work_dir, 'anno_stat*pdf'))
        for i in pdf_files:
            movefile(i, os.path.join(self.output_dir, '01_anno_Stat', os.path.basename(i)))

        pdf_files = glob.glob(os.path.join(self.work_dir, 'cog_bar*pdf'))
        for i in pdf_files:
            movefile(i, os.path.join(self.output_dir, '04_anno_COG', os.path.basename(i)))

        pdf_files = glob.glob(os.path.join(self.work_dir, 'pfam_bar*pdf'))
        for i in pdf_files:
            movefile(i, os.path.join(self.output_dir, '05_anno_Pfam', os.path.basename(i)))

        pdf_files = glob.glob(os.path.join(self.work_dir, 'subloc_bar*pdf'))
        for i in pdf_files:
            movefile(i, os.path.join(self.output_dir, '06_anno_SubLoc', os.path.basename(i)))

        '''
        '''
        repaths = [
            ['.', '', '蛋白注释结果目录', 0],
            ["./01_anno_Stat", "", "蛋白注释统计结果目录", 0, "220066"],
            ["./02_anno_GO", "", "蛋白GO注释结果目录", 0, "220067"],
            ["./04_anno_COG", "", "蛋白COG注释结果目录", 0, "220068"],
            ["./03_anno_KEGG", "", "蛋白KEGG注释结果目录", 0, "220069"],
            ["./05_anno_Pfam", "", "蛋白Pfam注释结果目录", 0],
            ["./06_anno_SubLoc", "", "蛋白亚细胞定位注释结果目录", 0, "220070"],
            ["./blast_xml", "", "蛋白比对注释结果目录", 1, "220071"]
        ]

        self.workflow_output_tmp = self._sheet.output
        if re.match(r'tsanger:', self.workflow_output_tmp):
            self.workflow_output = self.workflow_output_tmp.replace('tsanger:', '/mnt/ilustre/tsanger-data/')
        else:
            self.workflow_output = self.workflow_output_tmp.replace('sanger:', '/mnt/ilustre/data/')

        result_dir = self.add_upload_dir(target_dir)
        self.inter_dirs = [
            ["3_Anno", "", "蛋白注释结果目录", 0]
        ]
        result_dir.add_relpath_rules(repaths)
        result_dir.add_regexp_rules([
            [r".*anno_stat.pdf", "", "注释结果柱状图"],
            [r".*go_lev.*_pie.pdf", "", "饼图"],
            [r".*go_lev.*_bar.pdf", "", "直方图"],
            [r".*cog_bar.*pdf", "", "COG分类统计柱状图"],
            [r".*path_class.*pdf", "", "Pathway分类统计柱状图"],
            [r".*key_path.*pdf", "", "重要通路统计图"],
            [r".*subloc_bar.*pdf", "", "Subloc注释柱状图"],
            [r".*pfam_bar.*pdf", "", "Pfam注释柱状图"],
        ])
        db = Config().get_mongo_client(mtype="dia")[Config().get_mongo_dbname("dia")]
        col1 = db["sg_annotation_stat"]
        col1.update({"_id": ObjectId(self.option("stat_id"))}, {"$set": {"result_dir": self.workflow_output}}, upsert=True)

        super(ProteinAnnotationWorkflow, self).end()

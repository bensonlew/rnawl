# -*- coding: utf-8 -*-
# __author__ = 'gaohao'

"""RandomForest分析"""

from biocluster.workflow import Workflow
import os,re
from bson.objectid import ObjectId
import types
from mbio.packages.metagenomic.common import get_save_pdf_status, get_name, get_submit_loc
import json

class RandomforestWorkflow(Workflow):
    """
    报告中调用RandomForest分析时使用
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(RandomforestWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "update_info", "type": "string"},
            {"name": "group", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "pre_group", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "anno_table", "type": "infile", "format": "meta.profile"},
            {"name": "anno_type", "type": "string"},
            {"name": "geneset_id", "type": "string"},
            {"name": "anno_id", "type": "string"},
            {"name": "level_id", "type": "string"},
            {"name": "geneset_table", "type": "infile", "format": "meta.otu.otu_table"},
            {"name": "abund_file", "type": "infile", "format": "meta.otu.otu_table"},
            {"name": "gene_list", "type": "infile", "format": "meta.profile"},
            {"name": "level_type", "type": "string"},
            {"name": "level_type_name", "type": "string", "default": ""},
            {"name": "lowest_level", "type": "string", "default": ""},
            {"name": "group_detail", "type": "string"},
            {"name": "pre_group_detail", "type": "string"},
            {"name": "group_id", "type": "string"},
            {"name": "pre_group_id", "type": "string"},
            {"name": "auth_method", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "tree_num", "type": "int"},
            {"name": "problem_type", "type": "int", "default": 1},
            # {'name': 'save_pdf', 'type': 'int', 'default': 1}
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.abundance = self.add_tool("meta.create_abund_table")
        self.sort_samples = self.add_tool("meta.otu.sort_samples_mg")
        self.sort_samples2 = self.add_tool("meta.otu.sort_samples_mg")
        self.rf = self.add_tool("statistical.randomforest")
        self.list = [self.sort_samples, self.sort_samples2]

    def run_abundance(self):
        options = {
            'anno_table': self.option('anno_table'),
            'geneset_table': self.option("geneset_table"),
            'level_type': self.option('level_type'),
            'level_type_name': self.option('level_type_name'),
            'gene_list': self.option('gene_list'),
            'lowest_level': self.option('lowest_level')
        }
        self.abundance.set_options(options)
        if self.option("pre_group").is_set:
            self.abundance.on("end", self.run_sort_samples)
            self.abundance.on("end", self.run_tiqu_sample)
        else:
            self.abundance.on("end", self.run_sort_samples)
        self.abundance.run( )

    def run_sort_samples(self):
        if self.option("abund_file").is_set:
            abund_table = self.option("abund_file")
        else:
            abund_table = self.abundance.option("out_table").prop['path']
        self.sort_samples.set_options({
            "in_otu_table": abund_table,
            "group_table": self.option("group"),
            "sample_del_warn": "T"
        })
        self.sort_samples.run()

    def run_tiqu_sample(self):
        if self.option("abund_file").is_set:
            abund_table = self.option("abund_file")
        else:
            abund_table = self.abundance.option("out_table").prop['path']
        self.sort_samples2.set_options({
            "in_otu_table": abund_table,
            "group_table": self.option("pre_group"),
            "sample_del_warn": "F",
            "variable_del":'F',
        })
        self.sort_samples2.run()

    def run_rf(self):
        abund_table = self.sort_samples.option("out_otu_table").prop['path']
        self.rename_file(abund_table,self.work_dir + "/abund_table.xls",self.work_dir + "/names.xls")
        options={
            "abutable": self.work_dir + "/abund_table.xls",
            "grouptable": self.option("group"),
            "method": self.option("auth_method"),
            "ntree":self.option("tree_num"),
        }
        if self.option("pre_group").is_set:
            pre_table = self.sort_samples2.option("out_otu_table").prop['path']
            self.rename_file(pre_table, self.work_dir + "/pre_table.xls", self.work_dir + "/pre_names.xls")
            options['predict_sample']=self.work_dir + "/pre_table.xls"
        self.logger.info(options)
        self.rf.set_options(options)
        self.rf.on('end', self.set_db)
        self.rf.run()

    def end(self):
        # if self.option("save_pdf"):
        if self.pdf_status:
            if self.pdf_status:
                pdf_outs = self.figsave.output_dir
                os.system("cp -r {}/* {}/".format(pdf_outs, self.rf.output_dir))
        result_dir = self.add_upload_dir(self.rf.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "RandomForest分析结果目录", 0, "121008"],
            ["./randomForest_confusion_subRF_table.xls", "xls", "随机森林计算出的分类结果(根据挑选的重要性变量)", 0, "121009"],
            ["./randomForest_confusion_table.xls", "xls", "随机森林计算出的分类结果(所有变量)", 0, "121010"],
            ["./randomForest_imptance_table.xls", "xls", "所有物种（变量）的重要性衡量值", 0, "121011"],
            ["./randomForest_10-fold_CV.xls", "xls", "不同物种(变量)数下十折交叉验证的分类错误率表格", 0, "121012"],
            ["./randomForest_AUC.xls", "xls", "不同物种(变量)数下随机森林的AUC表", 0, "121013"],
            ["./randomForest_predict.xls", "xls", "随机森林预测样品结果表", 0, "121014"],
            ["./randomForest_subRF_pcoa_sites.xls", "xls", "根据挑选出的重要物种(变量)构建随机森林的Pcoa坐标", 0, "121015"],
            ["randomForest_top_vim.pdf", "pdf", "重要性柱形图"],
            ["AUC_evaluation.pdf", "pdf", "AUC验证结果图"],
            ["10fold_evaluation.pdf", "pdf", "十折交叉验证结果图"],
        ])
        result_dir.add_regexp_rules([[r'randomForest_top.*_vimp.xls$', 'xls', '重要物种（变量）丰度表格', 0, "121016"]])
        super(RandomforestWorkflow, self).end()

    def set_db(self):
        """
        保存两组比较分析的结果表保存到mongo数据库中
        """
        self.logger.info("正在写入mongo数据库")
        api_rf = self.api.api("metagenomic.randomforest")
        main_id = self.option("main_id")
        if not isinstance(main_id, ObjectId):
            if isinstance(main_id, types.StringTypes):
                main_id = ObjectId(main_id)
            else:
                self.set_error('main_id必须为ObjectId对象或其对应的字符串！', code="12802901")
        self.logger.info(main_id)
        files =os.listdir(self.rf.output_dir)
        for file in files:
            if re.search('randomForest_top[0-9]*_vimp.xls',file):
                self.rename_result(self.work_dir + "/names.xls",self.rf.output_dir + '/' + file,self.rf.output_dir + '/result.xls')
                os.remove(self.rf.output_dir + '/' + file)
                os.renames(self.rf.output_dir + '/result.xls',self.rf.output_dir + '/' + file)
        if os.path.exists(self.rf.output_dir + '/randomForest_subRF_pcoa_sites.xls'):
            api_rf.add_rf_scatter(main_id,self.rf.output_dir + '/randomForest_subRF_pcoa_sites.xls')
        else:
            self.set_error('randomForest_subRF_pcoa_sites.xls file is not exists！', code="12802902")
        if os.path.exists(self.rf.output_dir + '/randomForest_imptance_table.xls'):
            self.rename_result(self.work_dir + "/names.xls",self.rf.output_dir + '/randomForest_imptance_table.xls',self.rf.output_dir + '/randomForest_imptance_table.xls2')
            os.remove(self.rf.output_dir + '/randomForest_imptance_table.xls')
            os.renames(self.rf.output_dir + '/randomForest_imptance_table.xls2',self.rf.output_dir + '/randomForest_imptance_table.xls')
            api_rf.add_rf_bar(main_id,self.rf.output_dir + '/randomForest_imptance_table.xls')
        else:
            self.set_error('randomForest_imptance_table.xls file is not exists！', code="12802903")
        if self.option('auth_method') == "CV":
            if os.path.exists(self.rf.output_dir + '/randomForest_10-fold_CV.xls'):
                api_rf.add_rf_evaluate(main_id,self.rf.output_dir + '/randomForest_10-fold_CV.xls')
            else:
                self.set_error('randomForest_10-fold_CV.xls file is not exists！', code="12802904")
        if self.option('auth_method') == "AUC":
            if os.path.exists(self.rf.output_dir + '/randomForest_AUC.xls'):
                api_rf.add_rf_evaluate(main_id,self.rf.output_dir + '/randomForest_AUC.xls')
            else:
                self.set_error('randomForest_AUC.xls file is not exists！', code="12802905")
        if self.option('pre_group').is_set:
            if os.path.exists(self.rf.output_dir + '/randomForest_predict.xls'):
                api_rf.add_rf_evaluate(main_id, self.rf.output_dir + '/randomForest_predict.xls')
            else:
                self.set_error('randomForest_predict.xls file is not exists！', code="12802906")
        # if self.option("save_pdf"):
        self.pdf_status = get_save_pdf_status("_".join(self.sheet.id.split('_')[:2]))
        if self.pdf_status:
            name = get_name(self.option("main_id"), "randomforest")
            self.figsave = self.add_tool("metagenomic.fig_save")
            self.figsave.on('end', self.end)
            submit_loc = get_submit_loc(self.option("main_id"), "randomforest")
            self.figsave.set_options({
                "task_id": "_".join(self.sheet.id.split('_')[:2]),
                "table_id": self.option("main_id"),
                "table_name": name,
                "project": "metagenomic",
                "submit_loc": submit_loc,
                "interaction": 1,
            })
            self.figsave.run()
        else:
            self.end()

    def run(self):
        if self.option("abund_file").is_set:
            if self.option("pre_group").is_set:
                self.on_rely(self.list,self.run_rf)
                self.run_sort_samples()
                self.run_tiqu_sample()
            else:
                self.sort_samples.on("end",self.run_rf)
                self.run_sort_samples()
        else:
            if self.option("pre_group").is_set:
                self.on_rely(self.list,self.run_rf)
                self.run_abundance()
            else:
                self.sort_samples.on("end",self.run_rf)
                self.run_abundance()
        super(RandomforestWorkflow, self).run()

    def rename_file(self,input,output,output2):
        with open (input,"r") as f,open (output,"w") as p,open (output2,"w") as s:
            lines =f.readlines()
            p.write(lines[0])
            s.write("name" + '\t' + 'new_name' + '\n')
            n=1
            for lin in lines[1:]:
                du =lin.rstrip('\n\r').split('\t')
                des ='tax' + str(n)
                s.write(du[0]  + '\t' + des + '\n')
                du[0]=des
                p.write("\t".join(du) + "\n")
                n +=1

    def rename_result(self,input,input2,output):
        with open (input,"r") as f,open (input2,"r") as p,open (output,"w") as s:
            dict1={}
            lines =f.readlines()
            for lin in lines[1:]:
                du =lin.rstrip('\n\r').split('\t')
                dict1[du[1]]=du[0]
            lines2 = p.readlines()
            s.write(lines2[0])
            for li in lines2[1:]:
                duu = li.rstrip('\n\r').split('\t')
                if duu[0] in dict1:
                    duu[0]=dict1[duu[0]]
                des="\t".join(duu)
                s.write(des + "\n")
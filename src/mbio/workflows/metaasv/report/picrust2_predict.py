#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# @Author  :qingchen.zhang @20191225

import os,re
from biocluster.core.exceptions import OptionError
from biocluster.workflow import Workflow
from mbio.packages.metagbin.common_function import link_dir


class Picrust2PredictWorkflow(Workflow):
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(Picrust2PredictWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "otu_fasta", "type": "infile", "format": "sequence.fasta"}, ##代表序列OTU或者ASV序列
            {"name": "otu_table", "type": "infile", "format": "meta.otu.otu_table"},##代表序列的otu表
            {"name": "thread", "type": "string", "default": "8"}, ##线程数
            {"name": "database", "type": "string", "default": "COG,EC,KO"},##注释参考库
            {"name": "hsp", "type": "string", "default": "mp"}, ##HSP方法选择，主要有5中方法{mp, emp_prob, subtree_average,pic,scp}
            {"name": "analysis_type", "type": "string", "default": "16S"},##选择方法16S、18S、ITS
            {"name": "group_method", "type": "string", "default": ""}, ##是否进行分组合并
            {"name": "group", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "main_id", "type": "string"},##主表id
            {"name": "update_info", "type": "string"}, #用于框架更新状态
            {"name": "group_detail", "type": "string", "default": "none"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_AFTER_END = False
        self.sort_samples = self.add_tool("meta.otu.sort_samples_mg")
        self.predict = self.add_tool("meta.picrust2_predict")
        self.api_picrust2 = self.api.api("metaasv.picrust2_predict")

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option("otu_fasta").is_set:
            raise OptionError("请传入otu_reps序列！")
        if not self.option("otu_table").is_set:
            raise OptionError("请传入otu的丰度表文件！")

    def run_sort_samples(self):
        """
        对丰度表是否进行合并计算
        :return:
        """
        self.sort_samples.set_options({
            "in_otu_table": self.option("otu_table"),
            "group_table": self.option("group"),
            "method": self.option("group_method"),
        })
        self.sort_samples.run()

    def run_predict(self):
        """
        用picrust2进行预测计算
        :return:
        """
        self.check_params()
        if self.option("group_method") != "":
            otu_table = self.sort_samples.option("out_otu_table")
        else:
            otu_table = self.option("otu_table")
        opts = {
            "otu_fasta": self.option("otu_fasta"),
            "otu_table": otu_table,
            "thread": self.option("thread"),
            "database": self.option("database"),
            "hsp": self.option("hsp"),
            "analysis_type": self.option("analysis_type")
        }
        self.predict.set_options(opts)
        self.predict.run()

    def check_params(self):
        """
        对数据进行检查，具体逻辑为：
        如果工作流跑的数据为16s数据，则进行16S的运行，否则报错；
        如果工作流跑的数据为18s数据，则进行18S的运行，否则报错；
        如果工作流跑的数据为its数据，则进行ITS的运行，否则报错；
        :return:
        """
        task_id = self._sheet.id
        task_id_list = task_id.split("_")
        new_task_id = "_".join(task_id_list[0:2])
        info = self.api_picrust2.get_database(new_task_id)
        database = info['database']
        self.logger.info("正在对数据库进行检查")
        if (self.option("analysis_type") in ["16S"]) and (database not in ['silva123/16s_bacteria', 'silva123/16s_archaea', 'silva123/16s','silva119/16s_bacteria', 'silva119/16s_archaea', 'silva119/16s','silva128/16s_archaea', 'silva128/16s_bacteria', 'silva128/16s','silva132/16s_archaea', 'silva132/16s_bacteria', 'silva132/16s','greengenes135/16s', 'greengenes135/16s_archaea', 'greengenes135/16s_bacteria','rdp11.5/16s', 'rdp11.5/16s_bacteria', 'rdp11.5/16s_archaea', 'fgr/amoA', 'fgr/nosZ', 'fgr/nirK', 'fgr/nirS','fgr/nifH', 'fgr/pmoA', 'fgr/mmoX', 'fgr/mcrA','fgr/amoA_archaea', 'fgr/amoA_bacteria','maarjam081/AM', 'Human_HOMD', 'Human_HPB', 'Protist_PR2_v4.5', "nt", 'nt/16s','silva138/16s_archaea', 'silva138/16s_bacteria', 'silva138/16s']):
            raise OptionError("注释数据库%s不是16s数据库，不能进行Picrust2分析！", variables=(database))
        elif (self.option("analysis_type") in ["18S"]) and (database not in ["silva123/18s_eukaryota", 'silva119/18s_eukaryota', 'silva128/18s_eukaryota', 'silva132/18s_eukaryota','fgr/amoA', 'fgr/nosZ', 'fgr/nirK', 'fgr/nirS','fgr/nifH', 'fgr/pmoA', 'fgr/mmoX', 'fgr/mcrA','fgr/amoA_archaea', 'fgr/amoA_bacteria','maarjam081/AM', 'Human_HOMD', 'Human_HPB', 'Protist_PR2_v4.5', "nt",'nt/18S', 'silva138/18s_eukaryota']):
            raise OptionError("注释数据库%s不是18s数据库，不能进行Picrust2分析！", variables=(database))
        elif (self.option("analysis_type") in ["ITS"]) and (database not in ['unite8.0/its_fungi', 'unite7.2/its_fungi', 'unite7.0/its_fungi', 'fgr/amoA', 'fgr/nosZ', 'fgr/nirK', 'fgr/nirS','fgr/nifH', 'fgr/pmoA', 'fgr/mmoX', 'fgr/mcrA','fgr/amoA_archaea', 'fgr/amoA_bacteria','maarjam081/AM', 'Human_HOMD', 'Human_HPB', 'Protist_PR2_v4.5', "nt", 'nt/its']):
            raise OptionError("注释数据库%s不是its数据库，不能进行Picrust2分析！", variables=(database))
        else:
            self.logger.info("开始运行picrust2预测啦")

    def replace_dir(self,kegg_dir_old, metacyc_dir_old=None):
        """
        将丰度表中的所有丰度为0的在此基础上加上整张丰度表最小值的十分之一
        :return:
        """
        self.logger.info("正在将结果文件中的丰度值0替换")
        kegg_dir = kegg_dir_old
        if metacyc_dir_old:
            metacyc_dir = metacyc_dir_old
        files = ["prediction_enzyme.xls", "prediction_KO.xls", "prediction_module.xls", "prediction_pathway.L1.xls", "prediction_pathway.L2.xls", "prediction_pathway.L3.xls", "MetaCyc_pathway_pred.xls"]
        out_kegg_dir = os.path.join(self.output_dir, "KEGG")
        out_metacyc_dir = os.path.join(self.output_dir, "MetaCyc")
        for file in os.listdir(kegg_dir):
            if file in files:
                file_path = os.path.join(kegg_dir, file)
                out_file = os.path.join(out_kegg_dir, file)
                if os.path.exists(out_file):
                    os.remove(out_file)
                self.replace_file(file_path, out_file)
        if os.path.exists(metacyc_dir_old):
            for file in os.listdir(metacyc_dir):
                if file in files:
                    file_path = os.path.join(metacyc_dir, file)
                    out_file = os.path.join(out_metacyc_dir, file)
                    if os.path.exists(out_file):
                        os.remove(out_file)
                    self.replace_file(file_path, out_file)

    def replace_file(self, infile, outfile):
        """
        对文件进行替换0
        :param infile: 输入的kegg要做heatmap图的文件，取最小值绝对丰度的十分之一
        :param outfile: 输出的output文件
        :return:
        """
        second_list = []
        with open(infile, 'r') as f, open(outfile, 'w') as w:
            lines = f.readlines()
            w.write("{}".format(lines[0]))
            head_list = lines[0].split("\t")
            if head_list[0] in ["Enzyme_ID", "KO_ID", "Module_ID", "Level2"]:
                head_num = 2
            elif head_list[0] in ["level3"]:
                head_num = 4
            elif head_list[0] in ["Level1"]:
                head_num = 1
            else: #酶的表
                head_num = 2
            for line in lines[1:]:##求最小值
                line = line.strip().split('\t')
                min_list = []
                for i in range(head_num, len(line)):
                    if float(line[i]) != 0.0:
                        line_min = float(line[i])
                        min_list.append(line_min)
                min_line = min(min_list)
                second_list.append(min_line)
            min_table = float(min(second_list))##得到全表最小值
            for line in lines[1:]:
                line = line.strip().split('\t')
                data_list = []
                for i in range(head_num, len(line)):
                    if float(line[i]) == 0.0:
                        line[i] = float(line[i]) + float(min_table) / 10
                    else:
                        line[i] = float(line[i])
                    data_list.append(line[i])
                w.write("{}\t{}\n".format("\t".join(line[0:head_num]), "\t".join(str(i) for i in data_list)))

    def set_db(self):
        """
        将数据导入MongoDB
        :return:
        """
        self.logger.info("正在准备将picrust2预测结果进行替换")
        link_dir(self.predict.output_dir, self.output_dir)
        kegg_dir = os.path.join(self.predict.output_dir, 'KEGG')
        metacyc_dir = os.path.join(self.predict.output_dir, 'MetaCyc')
        self.replace_dir(kegg_dir, metacyc_dir)
        self.logger.info("正在写入mongo数据库")
        api_picrust2 = self.api_picrust2
        main_id = self.option("main_id")
        # main_id = api_picrust2.add_function_prediction(params=params, otu_id="5db24ebf17b2bf304ed8877e")
        ##导入COG注释结果
        cog_path = os.path.join(self.output_dir, "COG/prediction_cog.xls")
        function_path = os.path.join(self.output_dir, "COG/prediction_function.xls")
        if os.path.exists(cog_path):
            self.logger.info("开始导入cog预测结果")
            api_picrust2.add_cog_annotation(prediction_id=main_id, annotation_path=cog_path, type='cog')
        else:
            if self.option("analysis_type") in ['16s', "16S"]:
                self.set_error("未能成功正确的生成cog结果文件或结果文件为空，请检查！")
            else:
                self.logger.info("18S和ITS没有cog功能分析")
        if os.path.exists(function_path):
            self.logger.info("开始导入cog的function功能结果")
            api_picrust2.add_cog_annotation(prediction_id=main_id, annotation_path=function_path, type="function")
        else:
            if self.option("analysis_type") in ['16s', "16S"]:
                self.set_error("未能成功正确的生成cog结果文件或结果文件为空，请检查！",)
            else:
                self.logger.info("18S和ITS没有cog功能分析")

        ##导入KEGG注释结果
        ko_path = os.path.join(self.output_dir, "KEGG/prediction_KO.xls")
        module_path = os.path.join(self.output_dir, "KEGG/prediction_module.xls")
        enzyme_path = os.path.join(self.output_dir, "KEGG/prediction_enzyme.xls")
        if os.path.exists(ko_path):
            self.logger.info("正在导入ko表")
            api_picrust2.add_kegg_abu(prediction_id=main_id, abu_path=ko_path, type='ko')
            self.logger.info("开始更新ko表")
            api_picrust2.update_main_type(main_id=main_id, type='ko')
        else:
            if self.option("analysis_type") in ['16s', "16S"]:
                self.set_error("未能成功正确的生成ko结果文件或结果文件为空，请检查！")
            else:
                self.logger.info("18S和ITS没有ko功能分析")
        if os.path.exists(module_path):
            self.logger.info("正在导入module表")
            api_picrust2.add_kegg_abu(prediction_id=main_id, abu_path=module_path, type='module')
            api_picrust2.update_main_type(main_id=main_id, type='module')
        else:
            if self.option("analysis_type") in ['16s', "16S"]:
                self.set_error("未能成功正确的生成module结果文件或结果文件为空，请检查！")
            else:
                self.logger.info("18S和ITS没有module功能分析")
        if os.path.exists(enzyme_path):
            self.logger.info("正在导入enzyme表")
            api_picrust2.add_kegg_abu(prediction_id=main_id, abu_path=enzyme_path, type='enzyme')
            api_picrust2.update_main_type(main_id=main_id, type='enzyme')
            api_picrust2.update_specimen(sample_path=enzyme_path, prediction_id=main_id)
            self.logger.info("更新主表specimen成功！")
        else:
            self.set_error("%s未能成功正确的生成enzyme结果文件或结果文件为空，请检查！" , variables=(self.option("analysis_type")), code="12705105")

        pathway_path = os.path.join(self.output_dir, "KEGG")
        if os.path.exists(pathway_path):
            self.logger.info("正在导入pathway表")
            list_dirs = os.listdir(pathway_path)
            for file in list_dirs:
                file_path = os.path.join(pathway_path, file)
                file_name = os.path.basename(file_path)
                if re.search(r"pathway", file_name):
                    api_picrust2.add_kegg_level(prediction_id=main_id, kegg_path=pathway_path)
                    break
                else:
                    self.logger.info("18S和ITS没有pathway功能分析")

        ##导入Metacyc注释结果
        metacyc_path = os.path.join(self.output_dir, "MetaCyc/MetaCyc_pathway_pred.xls")
        if os.path.exists(metacyc_path):
            self.logger.info("正在导入metacyc表")
            api_picrust2.add_metacyc_abu(prediction_id=main_id, abu_path=metacyc_path)
            api_picrust2.update_main_type(main_id=main_id, type='metacyc')
        self.end()

    def end(self):
        """
        结束和上传结果文件
        :return:
        """
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "dir", "PICRUSt2功能预测结果目录",0,""],
            ["./COG/", "dir", "COG功能预测结果目录",0,""],
            ["./KEGG/", "dir", "KEGG功能预测结果目录",0,""],
            ["./MetaCyc/", "dir", "MetaCyc功能预测结果目录",0,""],
            [r"./COG/prediction_cog.xls", "xls", "COG样品丰度表",0,""],
            [r"./COG/prediction_function.xls", "xls", "Function样品丰度表",0,""],
            [r"./COG/weighted_nsti.xls", "xls", "COG功能预测的NSTI值文件",0,""],
            [r"./KEGG/prediction_KO.xls", "xls", "KO样品丰度表",0,""],
            [r"./KEGG/prediction_pathway.L1.xls", "xls", "代谢通路level 1丰度表",0,""],
            [r"./KEGG/prediction_pathway.L2.xls", "xls", "代谢通路level 2丰度表",0,""],
            [r"./KEGG/prediction_pathway.L3.xls", "xls", "代谢通路level 3丰度表",0,""],
            [r"./KEGG/kegg.pathway.profile.xls", "xls", "各个样品的pathway丰度统计",0,""],
            [r"./KEGG/prediction_enzyme.xls", "xls", "enzyme样品丰度表",0,""],
            [r"./KEGG/weighted_nsti.xls", "xls", "KEGG功能预测的NSTI值文件",0,""],
            [r"./KEGG/prediction_module.xls", "xls", "module样品丰度表",0,""],
            [r"./MetaCyc/MetaCyc_pathway_pred.xls", "xls", "MetaCyc pathway样品丰度表",0,""],
            [r"./COG_predicted.xls", "xls", "预测的OTU中COG的数量",0,""],
            [r"./EC_predicted.xls", "xls", "预测的OTU中enzyme的数量",0,""],
            [r"./KO_predicted.xls", "xls", "预测的OTU中KO的数量",0,""],
            [r"./marker_predicted_and_nsti.xls", "xls", "预测的16S拷贝数和NSTI值",0,""],
            [r"./out.tre", "tre", "参考序列树文件",0,""],
        ])
        super(Picrust2PredictWorkflow, self).end()

    def run(self):
        """
        运行
        :return:
        """
        if self.option("group_method") != "":
            self.sort_samples.on("end", self.run_predict)
            self.predict.on("end", self.set_db)
            self.run_sort_samples()
        else:
            self.predict.on("end", self.set_db)
            self.run_predict()
        super(Picrust2PredictWorkflow, self).run()

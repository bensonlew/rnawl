# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

from mbio.packages.metaasv.common_function import link_dir,normalize_data
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError
import pandas as pd
import re
import os


class RandomforestWorkflow(Workflow):
    """
    metaasv 随机森林分析
    特点：有环境因子
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(RandomforestWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "otutable", "type": "infile", "format": "meta.otu.otu_table, meta.otu.tax_summary_dir"},
            {"name": "level", "type": "int", "default": 9},
            {"name": "grouptable", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "predict_sample", "type": "infile", "format": "meta.otu.otu_table"},  # 用于预测丰度表样品分类文件
            {"name": "ntree", "type": "int", "default": 500},
            {"name": "problem_type", "type": "int", "default": 1},
            {"name": "method", "type": "string", "default": "CV"},
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "group_id", "type": 'string'},
            {"name": "group_id2", "type": 'string'},  # 用于“选择预测样本”
            {"name": "asv_id", "type": 'string'},
            {"name": "group_detail", "type": "string"},
            {"name": "group_detail2", "type": "string"},  # 用于“选择预测样本”
            {"name": "envtable", "type": "infile", "format": "meta.otu.otu_table"}, ##输入环境因子丰度表
            {"name": "env_labs", "type": "string", "default": ""},#选择的环境因子
            {"name": "env_id","type":"string"},## 输入环境因子的主表id
            {"name": "data_standard","type":"string", "default": ""},##数据归一化标准化方法
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_AFTER_END = False
        self.randomforest = self.add_tool("meta.beta_diversity.randomforest")
        # if self.option("data_standard") != "":
        #     self.otu_table = os.path.join(self.work_dir, "abundance.xls")
        # else:
        self.otu_table = self.option("otutable").prop["path"]

    def check_option(self):
        """
        参数检查
        :return:
        """
        if self.option('envtable').is_set:
            self.option('envtable').get_info()
            if self.option('envlabs'):
                labs = self.option('envlabs').split(',')
                for lab in labs:
                    if lab not in self.option('envtable').prop['group_scheme']:
                        raise OptionError('该envlabs中的因子不存在于环境因子表：%s', variables=(lab))
            else:
                pass
            env_table = self.option('envtable').prop['path']
            sample_e = []
            with open(env_table, 'r') as e:
                e.next()
                for line in e:
                    line = line.strip().split('\t')
                    sample_e.append(line[0])
                    for i in range(1, len(line)):
                        if float(line[i]) or line[i] == '0' or float(line[i]) == 0.0:
                            continue
                        else:
                            raise OptionError('环境因子表中存在分类型环境因子')

            otu_path = self.option("otutable").prop['path']
            with open(otu_path, 'r') as o:
                line = o.readline()
                line = line.strip('\n').split('\t')
                sample_o = line[1:]
                self.logger.info(sample_e)
                self.logger.info(sample_o)
            for i in sample_o:
                if i in str(sample_e):
                    continue
                else:
                    sample_name = i
                    self.logger.info(i)
                    raise OptionError('ASV表中的样本%s和环境因子表中的样本不一致，请剔除ASV表中非法样本！', variables=(sample_name))
            otu_table = pd.DataFrame(pd.read_table(otu_path, sep='\t', index_col=0))
            env_abund = pd.DataFrame(pd.read_table(env_table, sep='\t', index_col=0))
            name = env_abund.index.name
            env_abund.index = [str(i) for i in env_abund.index]
            new_env_abund = env_abund.ix[otu_table.columns]
            new_env_abund.index.name = name
            a = (new_env_abund != 0).any(axis=0)
            env_empt = []
            for i in range(len(new_env_abund.columns)):
                if a[i]:
                    pass
                else:
                    sample_name = new_env_abund.columns[i]
                    env_empt.append(sample_name)
            if env_empt:
                self.set_error('环境因子：%s在所选样品中均为0，请剔除该环境因子!', variables=(env_empt))
        if self.option("data_standard") not in ["absolute", "Relative", "Min-Max", "log10", "Z-score", ""]:
            self.set_error('数据归一化/标准化方法不存在，请检查方法名称是否正确!')

    def run_data_standard(self, otutable, abundance_path):
        """
        功能：对数据进行标准化处理
        :return:
        """
        self.logger.info("开始对数据进行标准化处理！")
        if os.path.exists(abundance_path):
            os.remove(abundance_path)
        try:
            normalize_data(otutable, abundance_path, abundance_method=self.option("data_standard"), type="asv")
            self.logger.info("对数据进行标准化处理结束！")
        except:
            self.set_error("未能成功的对数据进行标准化处理！")
        env_path = os.path.join(self.work_dir, "normalized_env.xls")
        if self.option('envtable').is_set:
            if os.path.exists(env_path):
                os.remove(env_path)
            try:
                normalize_data(self.option("envtable").prop["path"], env_path, abundance_method=self.option("data_standard"), type='env')
                self.logger.info("对数据进行标准化处理结束！")
            except:
                self.set_error("未能成功的对数据进行标准化处理！")

    def cat_env_and_asv(self, asv_file, env_file, out_file):
        """
        功能：将env表和asv表进行合并
        将env表转置成为asv表的一行数据进行分析
        :return:
        """
        asv_file = pd.read_table(asv_file, sep="\t", header=0)
        env_file = pd.read_table(env_file, sep="\t", header=0)
        asv_columns = list(asv_file.columns)
        sample_list = asv_columns[1:]##固定样本顺序
        env_columns = list(env_file.columns)
        env_columns[0] = "sample"
        env_file.columns = env_columns
        env_file.set_index("sample", inplace=True)
        new_env_file = env_file.T
        new_env_file.reset_index(inplace=True)
        new_env_file_columns = list(new_env_file.columns)
        new_env_file_columns[0] = "ASV ID"
        new_env_file.columns = new_env_file_columns
        new_env_file = new_env_file[["ASV ID"] + sample_list]##调整后的顺序
        result = pd.concat([asv_file, new_env_file], axis=0, join='outer')
        result.to_csv(out_file, sep="\t", index=0)
        return out_file

    def change_otuname(self, tablepath, file_name):
        """
        改换asv名称
        :param tablepath:
        :param file_name:
        :return:
        """
        newtable = self.work_dir + "/" + file_name + "_input_abund.xls"
        with open(tablepath, "r") as f, open(newtable, "w") as g:
            head = f.readline()
            g.write(head)
            for line in f:
                lines = line.split("\t", 1)
                specimen = re.subn("^.*; ", "", lines[0])[0]
                g.write(specimen + "\t" + lines[1])
        return newtable

    def run_randomforest(self):
        """
        运行 randomforest tool
        :return:
        """
        newtable = self.change_otuname(self.otu_table, "otutable")
        if self.option("data_standard") not in ["", "none", None]:
            abundance_path = os.path.join(self.work_dir, "abundance.xls")
            self.run_data_standard(newtable, abundance_path)
            otu_table = os.path.join(self.work_dir, "abundance.xls")
        else:
            otu_table = newtable
        if self.option('envtable').is_set and self.option("data_standard") not in ["", "none", None]:
            otu_env_file = os.path.join(self.work_dir, "asv_env_abund.xls")
            env_path = os.path.join(self.work_dir, "normalized_env.xls")
            out_file = self.cat_env_and_asv(otu_table, env_path, otu_env_file)
        elif self.option('envtable').is_set and self.option("data_standard") in ["", "none", None]:
            env_path = self.option('envtable').prop['path']
            otu_env_file = os.path.join(self.work_dir, "asv_env_abund.xls")
            out_file = self.cat_env_and_asv(otu_table, env_path, otu_env_file)
        else:
            out_file = otu_table
        options = {
            'otutable': out_file,
            'method': self.option('method'),
            'level': self.option('level'),
            'grouptable': self.option('grouptable'),
            'ntree': self.option('ntree'),
            'problem_type': self.option('problem_type'),
        }
        if self.option('predict_sample').is_set:
            newtable2 = self.change_otuname(self.option('predict_sample').prop['path'], "predict_sample")
            if self.option("data_standard") not in ["", "none", None]:
                abundance_file = os.path.join(self.work_dir, "abundance2.xls")
                self.run_data_standard(newtable2, abundance_file)
                otu_table2 = os.path.join(self.work_dir, "abundance2.xls")
            else:
                otu_table2 = newtable2
            if self.option('envtable').is_set and self.option("data_standard") not in ["", "none", None]:
                otu_env_file = os.path.join(self.work_dir, "asv_env_abund2.xls")
                env_path = os.path.join(self.work_dir, "normalized_env.xls")
                out_file2 = self.cat_env_and_asv(otu_table2, env_path, otu_env_file)
            elif self.option('envtable').is_set and self.option("data_standard") in ["", "none", None]:
                env_path = self.option('envtable').prop['path']
                otu_env_file = os.path.join(self.work_dir, "asv_env_abund2.xls")
                out_file2 = self.cat_env_and_asv(otu_table2, env_path, otu_env_file)
            else:
                out_file2 = newtable2
            options["predict_sample"] = out_file2
        self.randomforest.set_options(options)
        self.randomforest.on('end', self.set_db)
        self.randomforest.run()

    def end(self):
        """
        上传结果文件和结束
        :return:
        """
        repaths = [
            [".", "", "RandomForest分析结果目录", 0, ""],
            ["./randomForest_confusion_subRF_table.xls", "xls", "随机森林计算出的分类结果(根据挑选的重要性变量)", 0, ""],
            ["./randomForest_confusion_table.xls", "xls", "随机森林计算出的分类结果(所有变量)", 0, ""],
            ["./randomForest_imptance_table.xls", "xls", "所有物种（变量）的重要性衡量值", 0, ""],
            ["./randomForest_10-fold_CV.xls", "xls", "不同物种(变量)数下十折交叉验证的分类错误率表格", 0, ""],
            ["./randomForest_AUC.xls", "xls", "不同物种(变量)数下随机森林的AUC表", 0, ""],
            ["./randomForest_predict.xls", "xls", "随机森林预测样品结果表", 0, ""],
            ["./randomForest_subRF_pcoa_sites.xls", "xls", "根据挑选出的重要物种(变量)构建随机森林的Pcoa坐标", 0, ""],
        ]
        regexps = [
            [r'randomForest_top.*_vimp.xls$', 'xls', '重要物种（变量）丰度表格', 0, ""],
        ]
        sdir = self.add_upload_dir(self.output_dir)
        sdir.add_relpath_rules(repaths)
        sdir.add_regexp_rules(regexps)
        super(RandomforestWorkflow, self).end()

    def set_db(self):
        """
        连接结果文件和导入MongoDB
        :return:
        """
        link_dir(self.randomforest.output_dir, self.output_dir)
        api_randomforest = self.api.api("metaasv.randomforest")

        datadim = self.output_dir + '/randomForest_subRF_pcoa_sites.xls'
        datavip = self.output_dir + '/randomForest_imptance_table.xls'
        if self.option('method') == "AUC":
            datamethod = self.output_dir + '/randomForest_AUC.xls'
            api_randomforest.add_randomforest_evaluate(file_path=datamethod, table_id=self.option("main_id"))
        elif self.option('method') == "CV":
            datamethod = self.output_dir + '/randomForest_10-fold_CV.xls'
            api_randomforest.add_randomforest_evaluate(file_path=datamethod, table_id=self.option("main_id"))
        if self.option('predict_sample').is_set:
            predict = self.output_dir + '/randomForest_predict.xls'
            api_randomforest.add_randomforest_predict(file_path=predict, table_id=self.option("main_id"))
        if not os.path.isfile(datadim):
            self.logger.error("找不到报告文件:{}".format(datadim))
            self.set_error("找不到报告文件")
        if not os.path.isfile(datavip):
            self.logger.error("找不到报告文件:{}".format(datavip))
            self.set_error("找不到报告文件")
        api_randomforest.add_randomforest_dim(file_path=datadim, table_id=self.option("main_id"))
        api_randomforest.add_randomforest_vip(file_path=datavip, table_id=self.option("main_id"))
        self.end()

    def run(self):
        """
        运行
        :return:
        """
        self.run_randomforest()
        super(RandomforestWorkflow, self).run()

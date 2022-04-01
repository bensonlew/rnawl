# -*- coding: utf-8 -*-
# __author__ = 'qingchen.zhang'

from mbio.packages.metaasv.common_function import filter_asv_set
from mbio.packages.metaasv.common_function import link_dir,normalize_data
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError
import pandas as pd
import os
import re


class RocWorkflow(Workflow):
    """
    metaasv ROC曲线分析
    特点：有环境因子 增加数据合并方法
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(RocWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "otu_table", "type": "infile", "format": "meta.otu.otu_table, meta.otu.tax_summary_dir"},
            {"name": "level", "type": "int", "default": 9},
            {"name": "group_table", "type": "infile", "format": "meta.otu.group_table"},  # 分组表，只接受两个分组
            {"name": "method", "type": "string", "default": "sum"},  # 选定的各物种指标在组内求和、均值、中位数后进行roc计算
            {"name": "confidence_interval", "type": "float", "default": 0.95},##置信区间的阈值
            {"name": "top_n", "type": "int", "default": 0},
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"},
            {"name": "group_id", "type": 'string'},
            {"name": "otu_id", "type": 'string'},
            {"name": "asv_id", "type": 'string'},
            {"name": "group_detail", "type": "string"},
            {"name": "envtable", "type": "infile", "format": "meta.otu.otu_table"}, ##输入环境因子丰度表
            {"name": "env_labs", "type": "string", "default": ""},#选择的环境因子
            {"name": "env_id","type":"string"},
            {"name": "data_standard","type":"string", "default": ""},##数据归一化标准化方法
            {"name": "set_id", "type": "string", "default": ""},  ##筛选的ASV集的set_id
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_AFTER_END = False
        self.roc = self.add_tool("statistical.roc")
        # if self.option("data_standard") != "":
        #     self.otu_table = os.path.join(self.work_dir, "abundance.xls")
        # else:

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

            otu_path = self.option("otu_table").prop['path']
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
        改换名称
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

    def run_roc(self):
        """
        运行tool roc
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

        if self.option('top_n') != 0:
            top_n = self.option('top_n')
        else:
            top_n = 0
        options = {
            'abu_table': out_file,
            'group_table': self.option('group_table'),
            'method': self.option('method'),
            'top_n': top_n,
            'confidence_interval': str(self.option('confidence_interval')),
        }
        self.roc.set_options(options)
        self.roc.on('end', self.set_db)
        self.roc.run()

    def end(self):
        '''
        结束和上传结果文件目录
        '''
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "ROC分析结果目录", 0, ""],
            ["./roc_curve.xls", "xls", "ROC曲线结果表", 0, ""],
            ["./roc_auc.xls", "xls", "AUC值计算结果表", 0, ""],
            ["./roc_interval.xls", "xls", "ROC曲线置信区间", 0, ""],
            ["./best_loc.xls", "xls", "ROC曲线临界值", 0, ""],
            ["./roc_curve_smooth.xls", "xls", "ROC曲线平滑处理结果表", 0, ""]
        ])
        super(RocWorkflow, self).end()

    def set_db(self):
        """
        链接结果文件和导入MongoDB
        :return:
        """
        link_dir(self.roc.output_dir, self.output_dir)
        api_roc = self.api.api("metaasv.roc")
        datacurve = self.output_dir + '/roc_curve.xls'
        api_roc.add_roc_curve(roc_id=self.option("main_id"), type="curve", file=datacurve)

        if os.path.exists(self.output_dir + '/roc_curve_smooth.xls'):
            datacurve_s = self.output_dir + '/roc_curve_smooth.xls'
            api_roc.add_roc_curve(roc_id=self.option("main_id"), type="smooth", file=datacurve_s)

        if os.path.exists(self.output_dir + '/roc_interval.xls'):
            datainterval = self.output_dir + '/roc_interval.xls'
            api_roc.add_roc_interval(roc_id=self.option("main_id"), file=datainterval)

        dataauc = self.output_dir + '/roc_auc.xls'
        api_roc.add_roc_auc(roc_id=self.option("main_id"), file=dataauc, type="curve")

        if os.path.exists(self.output_dir + '/roc_auc_smooth.xls'):
            dataauc_s = self.output_dir + '/roc_auc_smooth.xls'
            api_roc.add_roc_auc(roc_id=self.option("main_id"), file=dataauc_s, type="smooth")
        databestloc = self.output_dir + '/best_loc.xls'
        if os.path.exists(databestloc):
            dataloc = self.work_dir + '/youden.xls'
            with open(databestloc, 'r') as f, open(dataloc, 'w') as w:
                lines = f.readlines()
                line = lines[1].strip().split("\t")
                y = (float(line[5]))/100.0
                x = (100-float(line[2]))/100.0
                value = y - x
                w.write(str(value))
            api_roc.add_roc_best_loc(roc_id=self.option("main_id"), file=databestloc, type="best_loc")
            api_roc.add_roc_youden(roc_id=self.option("main_id"), file=dataloc, type="youden_index")

        report_files = [datacurve, dataauc, databestloc]
        for f in report_files:
            if not os.path.isfile(f):
                self.logger.error("找不到报告文件:{}".format(f))
                self.set_error("找不到报告文件")

        self.end()

    def run(self):
        """
        运行
        :return:
        """
        if self.option("set_id") != "":
            try:
                filter_asv_set(self.option("otu_table").prop['path'], self.option("set_id"),self.option("asv_id"), self.option("level") ,os.path.join(self.work_dir, "otu_table2.xls"))
            except:
                self.set_error("筛选asv集失败！")
            self.option("otu_table", os.path.join(self.work_dir, "otu_table2.xls"))
        else:
            self.option("otu_table", os.path.join(self.work_dir, "otu_table.xls"))
        self.otu_table = self.option("otu_table").prop["path"]
        self.run_roc()
        super(RocWorkflow, self).run()        


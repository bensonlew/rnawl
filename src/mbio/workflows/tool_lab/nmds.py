# -*- coding: utf-8 -*-
# __author__ = 'linna'
# from biocluster.agent import Agent
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError
import os
import pandas as pd


class NmdsWorkflow(Workflow):
    """
    NMDS小工具工作流
    """
    METHOD = ['abund_jaccard', 'binary_chisq', 'binary_chord',
              'binary_euclidean', 'binary_hamming', 'binary_jaccard',
              'binary_lennon', 'binary_ochiai',
              'binary_pearson', 'binary_sorensen_dice',
              'bray_curtis', 'bray_curtis_faith', 'bray_curtis_magurran',
              'canberra', 'chisq', 'chord', 'euclidean', 'gower',
              'hellinger', 'kulczynski', 'manhattan', 'morisita_horn',
              'pearson', 'soergel', 'spearman_approx', 'specprof']

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(NmdsWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "tooltable", "type": "infile", "format": "tool_lab.table"},
            {"name": "specimen_name", "type": "string", "default": "column"},
            {"name": "group_table", "type": "infile", "format": "toolapps.group_table"},
            {"name": "scale", "type": "string", "default": "False"},
            {"name": "dis_method", "type": "string", "default": "bray_curtis"},
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())
        self.nmds = self.add_tool("tool_lab.nmds")

    def check_options(self):
        """
        参数检查
        """
        if not self.option('tooltable').is_set:
            raise OptionError('必须提供数据表', code="32702903")
        self.option('tooltable').get_info()
        if self.option('tooltable').prop['sample_num'] < 3:
            raise OptionError('列数少于3，不可进行分析', code="32702904")
        if self.option('dis_method') not in NmdsWorkflow.METHOD:
            raise OptionError('错误或者不支持该距离矩阵计算方法', code="32701903")
        # if self.option('group_table').is_set:
        #     if self.option('specimen_name') == 'column':
        #         for i in self.option('group_table').prop['sample_name']:
        #             if i not in self.option('tooltable').prop['col_sample']:
        #                 raise OptionError('分组文件中的样本不存在于表格中，查看是否数据选择错误', code="32702909")
        #     else:
        #         for i in self.option('group_table').prop['sample_name']:
        #             if i not in self.option('tooltable').prop['row_sample']:
        #                 raise OptionError('分组文件中的样本不存在于表格中，查看是否数据选择错误', code="32702910")
        return True

    def run_nmds(self):
        if self.option('group_table').is_set:
            options = {
                "tooltable": self.option('tooltable'),
                "specimen_name": self.option('specimen_name'),
                "group_table": self.option('group_table'),
                "scale": self.option('scale'),
                "dis_method": self.option('dis_method'),
            }
        else:
            options = {
                "tooltable": self.option('tooltable'),
                "specimen_name": self.option('specimen_name'),
                "scale": self.option('scale'),
                "dis_method": self.option('dis_method'),
            }
        self.nmds.set_options(options)
        self.nmds.on("end", self.set_output, "nmds")
        self.nmds.run()

    def set_output(self, event):
        obj = event['bind_object']
        if event['data'] == 'nmds':
            self.linkdir(obj.output_dir, 'nmds')
        self.set_db()

    def linkdir(self, dirpath, dirname):
        allfiles = os.listdir(dirpath)
        newdir = os.path.join(self.output_dir, dirname)
        if not os.path.exists(newdir):
            os.mkdir(newdir)
        oldfiles = [os.path.join(dirpath, i) for i in allfiles]
        newfiles = [os.path.join(newdir, i) for i in allfiles]
        for newfile in newfiles:
            if os.path.exists(newfile):
                if os.path.isfile(newfile):
                    os.remove(newfile)
                else:
                    os.system('rm -r %s' % newfile)
                    # self.logger.info('rm -r %s' % newfile)
        for i in range(len(allfiles)):
            if os.path.isfile(oldfiles[i]):
                os.link(oldfiles[i], newfiles[i])
            elif os.path.isdir(oldfiles[i]):
                # self.logger.info('cp -r %s %s' % (oldfiles[i], newdir))
                os.system('cp -r %s %s' % (oldfiles[i], newdir))

    def set_db(self):
        self.logger.info("开始导表")
        api_nmds = self.api.api("tool_lab.nmds_api")
        nmds_dir = self.output_dir + '/nmds'
        api_nmds.add_nmds_detail(self.option('main_id'), nmds_dir)
        api_nmds.add_nmds_ellipse_detail(self.option('main_id'), nmds_dir)
        api_nmds.add_nmds_table_detail(self.option('main_id'), nmds_dir)
        api_nmds.add_nmds_stress_detail(self.option('main_id'), nmds_dir)
        self.logger.info("导表结束")
        self.end()

    def run(self):
        """
        运行
        """
        self.run_nmds()
        super(NmdsWorkflow, self).run()

    def end(self):
        if self.option('group_table').is_set:
            result_dir = self.add_upload_dir(self.output_dir)
            result_dir.add_relpath_rules([
                [".", "", "Nmds分析结果输出目录"],
                ["./nmds_sites.xls", "xls", "样本坐标表"],
                ["./nmds_stress.xls", "xls", "样本特征拟合度值"],
                ["./group.xls", "xls", "分组表"],
                ["./ellipse.xls", "xls", "置信椭圆"],
            ])
        else:
            result_dir = self.add_upload_dir(self.output_dir)
            result_dir.add_relpath_rules([
                [".", "", "Nmds分析结果输出目录"],
                ["./nmds_sites.xls", "xls", "样本坐标表"],
                ["./nmds_stress.xls", "xls", "样本特征拟合度值"],
                ["./ellipse.xls", "xls", "置信椭圆"],
            ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(NmdsWorkflow, self).end()

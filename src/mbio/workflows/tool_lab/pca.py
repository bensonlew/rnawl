# !usr/bin/python
# -*- coding: utf-8 -*-
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError
import pandas as pd
import os
import subprocess
import re


class PcaWorkflow(Workflow):
    """
    pca工作流
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(PcaWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "tooltable", "type": "infile", "format": "tool_lab.table"},
            {"name": "specimen_name", "type": "string", "default": "column"},
            {"name": "group_table", "type": "infile", "format": "tool_lab.group_table"},
            {"name": "scale", "type": "string", "default": "True"},
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())
        self.pca = self.add_tool("tool_lab.pca")

    def check_options(self):
        """
        检查参数
        """
        if not self.option('tooltable').is_set:
            raise OptionError('必须提供数据表', code="32702903")
        self.option('tooltable').get_info()
        if self.option('tooltable').prop['sample_num'] < 3:
            raise OptionError('列数少于3，不可进行分析', code="32702904")
        # if self.option('group_table').is_set:
        #     if self.option('specimen_name') == '列标签':
        #         for i in self.option('group_table').prop['sample_name']:
        #             if i not in self.option('tooltable').prop['col_sample']:
        #                 raise OptionError('分组文件中的样本不存在于表格中，查看是否数据选择错误', code="32702909")
        #     else:
        #         for i in self.option('group_table').prop['sample_name']:
        #             if i not in self.option('tooltable').prop['row_sample']:
        #                 raise OptionError('分组文件中的样本不存在于表格中，查看是否数据选择错误', code="32702910")
        return True

    def run_pca(self):
        if self.option('group_table').is_set:
            options = {
                "tooltable": self.option('tooltable'),
                "specimen_name": self.option('specimen_name'),
                "group_table": self.option('group_table'),
                "scale": self.option('scale')
            }
        else:
            options = {
                "tooltable": self.option('tooltable'),
                "specimen_name": self.option('specimen_name'),
                "scale": self.option('scale')
            }
        self.pca.set_options(options)
        self.pca.on("end", self.set_output, "pca")
        self.pca.run()

    def set_output(self, event):
        obj = event['bind_object']
        if event['data'] == 'pca':
            self.linkdir(obj.output_dir, 'pca')
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
        """
        保存结果到Mongo库
        """
        self.logger.info("保存结果到Mongo")
        main_id = self.option("main_id")
        pca_api = self.api.api("tool_lab.pca")
        pca_dir = self.output_dir + '/pca'
        pca_api.add_pca_detail(main_id, pca_dir)
        pca_api.add_pca_ellipse_detail(main_id, pca_dir)
        pca_api.add_pca_box_detail(main_id, pca_dir)
        pca_api.add_pca_result_detail(main_id, pca_dir)
        pca_api.add_pca_importance_detail(main_id, pca_dir)
        pca_api.add_pca_rotation_detail(main_id, pca_dir)
        self.logger.info("保存结果到Mongo结束")
        self.end()

    def run(self):
        """
        运行
        """
        self.run_pca()
        super(PcaWorkflow, self).run()

    def end(self):
        if self.option('group_table').is_set:
            result_dir = self.add_upload_dir(self.output_dir)
            result_dir.add_relpath_rules([
                [".", "", "PCA分析结果输出目录"],
                ["./pca_importance.xls", "xls", "主成分解释度表"],
                ["./pca_rotation.xls", "xls", "主成分贡献度表"],
                ["./pca_rotation_all.xls", "xls", "主成分贡献度表全部"],
                ["./pca_sites.xls", "xls", "样本坐标表"],
                ["./group.xls", "xls", "分组表"],
                ["./ellipse.xls", "xls", "置信椭圆"],
            ])
        else:
            result_dir = self.add_upload_dir(self.output_dir)
            result_dir.add_relpath_rules([
                [".", "", "PCA分析结果输出目录"],
                ["./pca_importance.xls", "xls", "主成分解释度表"],
                ["./pca_rotation.xls", "xls", "主成分贡献度表"],
                ["./pca_rotation_all.xls", "xls", "主成分贡献度表全部"],
                ["./pca_sites.xls", "xls", "样本坐标表"],
                ["./ellipse.xls", "xls", "置信椭圆"],
            ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(PcaWorkflow, self).end()

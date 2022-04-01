# -*- coding: utf-8 -*-
# __author__ = 'zhaobinbin'
# modified 20200427

import os
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError


class RegressionWorkflow(Workflow):
    """
    相关性小工具工作流
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(RegressionWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "tooltable", "type": "infile", "format": "toolapps.table"},
            {"name": "group_table", "type": "infile", "format": "toolapps.group_table"},
            {"name": "scale", "type": "string", "default": "F"},
            {"name": "themename", "type": "string", "default": "线性回归分析"},
            {"name": "main_id", "type": "string"},
            {"name": "specimen_name", "type": "string", "default": "列标签"},
            {'name': "update_info", 'type': 'string'},
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())
        self.regression = self.add_tool("tool_lab.regression")

    def check_options(self):
        if not self.option("tooltable"):
            raise OptionError("必须输入数据表")
        if not self.option("group_table"):
            raise OptionError("必须输入分组文件")
        if not self.option("scale"):
            raise OptionError("必须明确是否标准化")
        if not self.option("specimen_name"):
            raise OptionError("必须明确是列标签还是行标签")
        return True

    def run_regression(self):
        options = {
            "tooltable": self.option('tooltable'),
            "group_table": self.option('group_table'),
            "scale": self.option('scale'),
            'themename': self.option('themename'),
            "specimen_name": self.option('specimen_name'),
        }
        self.regression.set_options(options)
        self.regression.on("end", self.set_output, "regression")
        self.regression.run()

    def set_output(self, event):
        obj = event['bind_object']
        if event['data'] == 'regression':
            self.linkdir(obj.output_dir, 'regression')
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

    def formattable(self, tablepath):
        """
        转置表格
        """
        this_table = tablepath
        if self.option('specimen_name') != '列标签':
            return this_table
        else:
            newtable = this_table + '.T'
            self.t_table(this_table, newtable)
            return newtable

    def t_table(self, table_file, new_table):
        """
        转置表格实现
        """
        with open(table_file) as f, open(new_table, 'w') as w:
            table_list = [i.rstrip().split('\t') for i in f.readlines()]
            table_list = map(lambda *a: '\t'.join(a) + '\n', *table_list)
            w.writelines(table_list)

    def set_db(self):
        self.logger.info("保存结果到mongo")
        self.logger.info(self.option('main_id'))
        self.logger.info(os.path.join(self.output_dir, "regression"))
        self.logger.info(self.formattable(self.option('tooltable').path))
        api_regression = self.api.api('tool_lab.regression')
        api_regression.add_regression(self.option('main_id'), os.path.join(self.output_dir, "regression"),
                                      self.option('tooltable').path)

        self.logger.info("导表结束")
        self.end()

    def run(self):
        self.run_regression()
        super(RegressionWorkflow, self).run()

    def end(self):
        if self.option('group_table'):
            result_dir = self.add_upload_dir(self.output_dir)
            result_dir.add_relpath_rules([
                [".", "", "线性回归分析结果输出目录"],
                ["./regression_data.xls", "xls", "数据表"],
                ["./regression_message.xls", "xls", "信息表"],
                ["./group.xls", "xls", "分组表"],
            ])
        else:
            result_dir = self.add_upload_dir(self.output_dir)
            result_dir.add_relpath_rules([
                [".", "", "线性回归分析结果输出目录"],
                ["./regression_data.xls", "xls", "数据表"],
                ["./regression_message.xls", "xls", "信息表"],
            ])
        super(RegressionWorkflow, self).end()

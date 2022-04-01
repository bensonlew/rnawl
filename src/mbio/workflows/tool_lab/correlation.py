# -*- coding: utf-8 -*-
# __author__ = 'HD'
# modified 20200413

import os
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError


class CorrelationWorkflow(Workflow):
    """
    相关性小工具工作流
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(CorrelationWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "corr_table", "type": "infile", "format": "toolapps.table"},
            {"name": "strategy", "type": "string", "default": "sample"},
            {"name": "sample", "type": "string", "default": "col"},
            {"name": "group_table", "type": "infile", "format": "toolapps.group_table"},
            {"name": "group_method", "type": "string", "default": "mean"},
            {"name": "method", "type": "string", "default": "pearson"},
            {"name": "distance_method", "type": "string", "default": "euclidean"},
            {"name": "cluster", "type": "string", "default": "complete"},
            {"name": "cluster_tree", "type": "string", "default": "True"},
            {"name": "log", "type": "string", "default": "False"},  # 去掉了取log功能，默认不取log
            {'name': 'main_id', 'type': 'string'},
            {'name': "update_info", 'type': 'string'},
            # {'name': "title", "type": "string"}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())
        self.correlation = self.add_tool("tool_lab.correlation")

    def check_options(self):
        if self.option("strategy") not in ["group", "sample"]:
            raise OptionError("不支持该分析策略")
        if self.option("method") not in ["pearson", "spearman"]:
            raise OptionError("不支持该相关系数方法")
        if self.option("distance_method") not in ["euclidean", "maximum", "manhattan",
                                                  "canberra", "binary", "minkowski"]:
            raise OptionError("不支持该距离算法")
        if self.option("cluster") not in ["complete", "single", "average", 'no']:
            raise OptionError("不支持该聚类方式")
        return True

    def run_correlation(self):
        options = {
            "corr_table": self.option('corr_table'),
            "strategy": self.option('strategy'),
            "sample": self.option('sample'),
            'group_table': self.option('group_table'),
            "group_method": self.option('group_method'),
            "method": self.option('method'),
            "distance_method": self.option('distance_method'),
            "cluster": self.option('cluster'),
            "cluster_tree": self.option('cluster_tree'),
            "log": self.option('log')
        }
        self.correlation.set_options(options)
        self.correlation.on("end", self.set_output, "correlation")
        self.correlation.run()

    def set_output(self, event):
        obj = event['bind_object']
        if event['data'] == 'correlation':
            self.linkdir(obj.output_dir, 'correlation')
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
        self.logger.info("保存结果到mongo")
        api_correlation = self.api.api("tool_lab.correlation")
        if self.option('cluster_tree') in ["True"] and self.option('cluster') != 'no':
            tree_col = self.output_dir + "/correlation/corr_col.tre"
            tree_row = self.output_dir + "/correlation/corr_row.tre"
        else:
            tree_col = None
            tree_row = None
        api_correlation.insert_detail_table(self.option("main_id"),
                                            self.output_dir + "/correlation/correlation_matrix.xls",
                                            tree_col, tree_row, self.output_dir + "/correlation/pvalue_matrix.xls")
        self.end()

    def run(self):
        self.run_correlation()
        super(CorrelationWorkflow, self).run()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "Correlation计算结果输出目录"],
            ["./correlation_matrix.xls", "xls", "Correlation矩阵"],
            ["./pvalue_matrix.xls", "xls", "Pvalue矩阵"],
            ["./cluster.tre", "tre", "聚类结果"]
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(CorrelationWorkflow, self).end()

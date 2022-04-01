# -*- coding: utf-8 -*-

from biocluster.workflow import Workflow
import pandas as pd
import glob


class ChiSquareWorkflow(Workflow):


    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(ChiSquareWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "table", "type": "infile", "format": "tool_lab.table"},
            #{'name': 'side_type', 'type': 'string', 'default': 'two.side'},  # 单尾或双尾检验 two.side,less,greater
            {'name': 'correct', 'type': 'string', 'default':'bonferroni'},
            {"name": "update_info", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "sample1", "type": "string"},
            {"name": "sample2", "type": "string"},
            {"name": "side_type","type": "string"}  #没有用到该参数
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())


    def run(self):
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        self.logger.info("start run workflow")

        self.run_diff_sample()

        super(ChiSquareWorkflow, self).run()


    def run_diff_sample(self):
        self.two_sample_tool = self.add_tool("tool_lab.metastat2")
        correct = self.option("correct")
        if correct == 'FDR':
            correct = 'fdr'

        options = {
            "chi_input": self.option("table").path,
            "chi_sample1": self.option("sample1"),
            "chi_sample2": self.option("sample2"),
            "chi_correction": correct,
            "test": 'chi',
            "chi_methor": 'DiffBetweenPropAsymptotic',
            "chi_coverage": 0.95,
        }

        self.two_sample_tool.set_options(options)
        self.two_sample_tool.on('end',self.set_db)
        self.two_sample_tool.run()


    def set_db(self):
        """
        保存结果到mongo数据库中
        """
        result = glob.glob(self.two_sample_tool.output_dir+"/*_result.xls")[0]
        data = pd.read_table(result,sep='\t',index_col=0)
        col1 = data.columns[0]
        col2 = data.columns[1]
        data['odds_ratio'] = data[col1]/data[col2]
        data.to_csv(result, sep='\t')

        api_name = self.api.api("tool_lab.chi_square")
        self.logger.info("正在写入mongo数据库")
        main_id = self.option("main_id")
        api_name.add_two_sample_detail(self.two_sample_tool.output_dir,main_id, self.option("sample1"), self.option("sample2"))
        self.end()


    def end(self):
        result_dir = self.add_upload_dir(self.two_sample_tool.output_dir)
        relpath_rules = [
            [".", "", "比较结果文件夹", 0, ],
        ]
        regexps = [
            [r".*/.*_result.xls", "xls", "显著性比较结果表，包括均值，标准差，p值", 0],
            [r".*/.*_CI.xls", "xls", "两样本比较的置信区间值", 0]
        ]
        result_dir.add_relpath_rules(relpath_rules)
        result_dir.add_regexp_rules(regexps)
        super(ChiSquareWorkflow, self).end()


if __name__ == '__main__':
    from biocluster.wsheet import Sheet
    data = {
        'name': 'test_sample_diff',
        'id': 'tsg_36964',
        'type': 'workflow',
        'options': {
            "table" : "/mnt/ilustre/users/sanger-dev/workspace/20200228/Metabolome_tsg_36991/Preprocess/output/pos/metab_abund.txt",
            "main_id" : "5e9e6a6017b2bf2049a81be9",
            "sample1" : "NG_D1_A",
            "sample2" : "NG_D1_B"
        }
    }

    wsheet = Sheet(data=data)

    wf = ChiSquareWorkflow(wsheet)
    wf.run()
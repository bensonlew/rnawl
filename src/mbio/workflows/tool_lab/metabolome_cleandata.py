# -*- coding: utf-8 -*-


from biocluster.workflow import Workflow
import os
import shutil


class MetabolomeCleandataWorkflow(Workflow):
    """
    代谢质控
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(MetabolomeCleandataWorkflow, self).__init__(wsheet_object)
        options = [
            {"name" : "pos_table", "type":"infile","format":"tool_lab.simple"},
            {"name" : "neg_table", "type":"infile","format":"tool_lab.simple"},
            {"name" : "order_file", "type":"infile","format":"tool_lab.simple"},
            {"name": "column", "type": "string", "default":"id"},
            {"name": "qc_filter", "type": "float", "default":0.5},
            {"name": "sample_filter", "type": "float", "default":0.5},
            {"name": "method", "type": "string", "default":"svr"},
            {"name": "update_info", "type": "string"},
            {"name": "metab_name", "type": "string"},
            {"name": "main_id", "type": "string"}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())


    def run(self):
        self.clean_tool = self.add_tool("tool_lab.metabolome_cleandata")
        options = {
            "pos_table": self.option('pos_table'),
            "order_file" : self.option("order_file"),
            "column" : self.option("column"),
            "qc_filter" : self.option("qc_filter"),
            "sample_filter" : self.option("sample_filter"),
            "method" :self.option("method")
        }
        if self.option("neg_table").is_set:
            options['neg_table'] = self.option("neg_table")

        self.clean_tool.set_options(options)
        self.clean_tool.on('end',self.set_db)

        self.clean_tool.run()

        super(MetabolomeCleandataWorkflow, self).run()

    def set_db(self):
        if os.path.exists(self.output_dir):
            shutil.rmtree(self.output_dir)
        shutil.copytree(self.clean_tool.output_dir, self.output_dir)

        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "代谢数据校正", 0]

        ])

        result_dir.add_regexp_rules([
            [".*", "", "", 0],
        ])

        super(MetabolomeCleandataWorkflow, self).end()

if __name__ == "__main__":
    from biocluster.wsheet import Sheet
    data = {
        'name': 'tool_lab.metabolome_cleandata',
        'id': 'tsg_36964_3',
        'type': 'workflow',
        'options': {
            "pos_table" : "/mnt/ilustre/users/sanger-dev/sg-users/zouguanqing/small_tool/metabolome_cleandata/pos_measurement.xls",
            "neg_table" : "/mnt/ilustre/users/sanger-dev/sg-users/zouguanqing/small_tool/metabolome_cleandata/neg_measurement.xls",
            "order_file" : "/mnt/ilustre/users/sanger-dev/sg-users/zouguanqing/small_tool/metabolome_cleandata/order.txt",
            "main_id" : "5e4ce5e817b2bf4b326f4b20"
        }
    }

    wsheet = Sheet(data=data)

    wf = MetabolomeCleandataWorkflow(wsheet)
    wf.run()

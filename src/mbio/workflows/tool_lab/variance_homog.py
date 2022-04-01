# -*- coding: utf-8 -*-


from biocluster.workflow import Workflow
import os
import shutil


class VarianceHomogWorkflow(Workflow):
    """
    tax4fun
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(VarianceHomogWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "update_info", "type": "string"},
            {"name": "table", "type": "infile", "format": "tool_lab.simple"},
            {"name": "group", "type": "infile", "format": "tool_lab.simple"},
            {"name": "method", "type": "string", "default": "bartlett"}, #bartlett\levene\ligner-killeen
            {"name": "main_id", "type": "string"}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())



    def run(self):

        self.tool = self.add_tool("tool_lab.variance_homog")
        options = {
            "table": self.option('table'),
            "group" : self.option('group'),
            "method" : self.option ('method')
        }
        self.tool.set_options(options)
        self.tool.on('end',self.set_db)

        self.tool.run()
        super(VarianceHomogWorkflow, self).run()

    def set_db(self):
        """
        """

        if os.path.exists(self.output_dir):
            shutil.rmtree(self.output_dir)
        shutil.copytree(self.tool.output_dir, self.output_dir)


        api = self.api.api('tool_lab.common')
        api.update_main('variance_homog', self.option('main_id'), {"main_id":self.option('main_id')})
        self.end()


    def end(self):

        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "方差齐性检验结果", 0]
        ])

        result_dir.add_regexp_rules([
            [".*/.*", "", "结果文件", 0],

        ])

        super(VarianceHomogWorkflow, self).end()

if __name__ == "__main__":
    from biocluster.wsheet import Sheet
    data = {
        'name': 'tool_lab.table2biom',
        'id': 'tsg_36964',
        'type': 'workflow',
        'options': {
            "main_id" : "5e4ce5e817b2bf4b326f4b20",
            "table" : "/mnt/ilustre/users/sanger-dev/sg-users/zouguanqing/small_tool/variance_homog/DATA.txt",
            "group" : "/mnt/ilustre/users/sanger-dev/sg-users/zouguanqing/small_tool/variance_homog/GROUP.txt"
            ##/mnt/ilustre/users/sanger-dev/i-sanger_workspace2/20210309/VarianceHomog_majorbio_293423_79500_198484
        }
    }

    wsheet = Sheet(data=data)

    wf = VarianceHomogWorkflow(wsheet)
    wf.run()

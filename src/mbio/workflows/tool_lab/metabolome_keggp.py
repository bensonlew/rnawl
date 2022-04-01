# -*- coding: utf-8 -*-


from biocluster.workflow import Workflow
import os
import shutil


class MetabolomeKeggpWorkflow(Workflow):
    """
    代谢kegg通路
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(MetabolomeKeggpWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "update_info", "type": "string"},
            {"name": "metab_name", "type": "string"},
            {"name": "main_id", "type": "string"}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())


    def run(self):
        self.anno = self.add_tool("tool_lab.anno_keggp")
        options = {
            "metab_name": self.option('metab_name')
        }

        self.anno.set_options(options)
        self.anno.on('end',self.set_db)

        self.anno.run()

        super(MetabolomeKeggpWorkflow, self).run()

    def set_db(self):
        """
        导表程序
        """
        anno_api = self.api.api('tool_lab.anno_keggp')

        anno_api.add_stat_detail(self.option('main_id'), self.anno.option('stat_out').path, self.anno.option('detail_out').path)

        if os.path.exists(self.output_dir):
            shutil.rmtree(self.output_dir)
        shutil.copytree(self.anno.output_dir, self.output_dir)
        anno_api.update_main(self.option('main_id'))

        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "KEGG通路注释结果", 0]

        ])

        result_dir.add_regexp_rules([
            [".*/detail.xls", "xls", "KEGG代谢物详情表", 0],
            [".*/stat.xls", "xls", "KEGG代谢物统计表", 0]
        ])

        super(MetabolomeKeggpWorkflow, self).end()

if __name__ == "__main__":
    from biocluster.wsheet import Sheet
    data = {
        'name': 'tool_lab.metabolome_keggp',
        'id': 'tsg_36964',
        'type': 'workflow',
        'options': {
            "metab_name" : "Nicotinamide adenine dinucleotide;;;NADH",
            "main_id" : "5e4ce5e817b2bf4b326f4b20"
        }
    }

    wsheet = Sheet(data=data)

    wf = MetabolomeKeggpWorkflow(wsheet)
    wf.run()

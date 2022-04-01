# -*- coding: utf-8 -*-


from biocluster.workflow import Workflow
import os
import shutil


class MetabolomeHmdbWorkflow(Workflow):
    """
    代谢kegg通路
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(MetabolomeHmdbWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "update_info", "type": "string"},
            {"name": "metab_name", "type": "string"},
            {"name": "main_id", "type": "string"}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())


    def run(self):
        self.anno = self.add_tool("tool_lab.anno_hmdb")
        options = {
            "metab_name": self.option('metab_name')
        }
        self.anno.set_options(options)
        self.anno.on('end',self.set_db)

        self.anno.run()

        super(MetabolomeHmdbWorkflow, self).run()

    def set_db(self):
        """
        导表程序
        """
        anno_api = self.api.api('tool_lab.anno_hmdb')

        anno_api.add_detail(self.option('main_id'),self.anno.option('out_table').path)
        anno_api.add_stat(self.option('main_id'),self.anno.output_dir+'/HmdbSuperclass.xls','superclass')
        anno_api.add_stat(self.option('main_id'),self.anno.output_dir+'/HmdbClass.xls','class')
        anno_api.add_stat(self.option('main_id'),self.anno.output_dir+'/HmdbSubclass.xls','subclass')

        if os.path.exists(self.output_dir):
            shutil.rmtree(self.output_dir)
        shutil.copytree(self.anno.output_dir, self.output_dir)
        anno_api.update_main(self.option('main_id'))


        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "HMDB注释结果", 0]

        ])

        result_dir.add_regexp_rules([
            [".*/detail.xls", "xls", "HMDB详情表", 0],
            [".*/HmdbSuperclass.xls", "xls", "HMDB Superclass统计表", 0],
            [".*/HmdbClass.xls", "xls", "HMDB Class统计表", 0],
            [".*/HmdbSubclass.xls", "xls", "HMDB Subclass统计表", 0]
        ])

        super(MetabolomeHmdbWorkflow, self).end()

if __name__ == "__main__":
    from biocluster.wsheet import Sheet
    data = {
        'name': 'tool_lab.metabolome_hmdb',
        'id': 'tsg_36964',
        'type': 'workflow',
        'options': {
            "metab_name" : "1-Methylhistidine;;;1,3-Diaminopropane",
            "main_id" : "5e4ce5e817b2bf4b326f4b20"
        }
    }

    wsheet = Sheet(data=data)

    wf = MetabolomeHmdbWorkflow(wsheet)
    wf.run()

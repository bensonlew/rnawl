# -*- coding: utf-8 -*-


from biocluster.workflow import Workflow
import os
import shutil


class KeggcWorkflow(Workflow):
    """
    化合物注释
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(KeggcWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "table", "type": "infile", "format": "tool_lab.simple"},
            {"name": "update_info", "type": "string"},
            #{"name": "name", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "search_type", "type": "string", "default":"entry"}  #entry （对应compound id）； cas_id； formula；name

        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())
        self.database_list = ["CBR","BP","EDC",'P','PC','L']

        self.database_desc ={
            "CBR" : "Compounds_Biological",
            "BP" : "Bioactive_Peptides",
            "EDC": "Endocrine_Disrupting",
            "P" : "Pesticides" ,
            "PC" : "Phytochemical",
            "L" : "Lipids"
        }

    def run(self):
        self.tool_map = {}
        for database in self.database_list:
            anno = self.add_tool("tool_lab.anno_keggc")
            options = {
                "metab_table": self.option('table'),
                "database_name" : database,
                "search_type" : self.option("search_type")
            }

            anno.set_options(options)
            self.tool_map[database] = anno

        self.on_rely(self.tool_map.values(), self.set_db)
        for tool in self.tool_map.values():
            tool.run()

        super(KeggcWorkflow, self).run()

    def set_db(self):
        """
        导表程序
        """
        anno_api = self.api.api('tool_lab.anno_keggc')

        true_databases = []
        for database in self.database_list:
            anno_tool = self.tool_map[database]
            ret = anno_api.add_stat_detail(self.option('main_id'), anno_tool.option('stat_out').path, database)
            if ret != 'no data':
                true_databases.append(database)
            else:
                continue
            database_desc = self.database_desc[database]
            if os.path.exists(self.output_dir+'/'+database_desc):
                shutil.rmtree(self.output_dir+'/'+database_desc)
            shutil.copytree(anno_tool.output_dir, self.output_dir+'/'+database_desc)

        anno_api.update_main(self.option('main_id'),true_databases)


        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "KEGG化合物分类注释结果", 0]

        ])

        result_dir.add_regexp_rules([
            [".*/level.xls", "xls", "化合物分类层级表", 0],
            [".*/stat.xls", "xls", "化合物分类统计表", 0]
        ])

        super(KeggcWorkflow, self).end()

if __name__ == "__main__":
    from biocluster.wsheet import Sheet
    data = {
        'name': 'keggc',
        'id': 'tsg_36964',
        'type': 'workflow',
        'options': {
            "table" : "/mnt/ilustre/users/sanger-dev/workspace/20200228/Metabolome_tsg_36991/Preprocess/output/pos/c.list",
            "main_id" : "5e4ce5e817b2bf4b326f4b20",

        }
    }

    wsheet = Sheet(data=data)

    wf = KeggcWorkflow(wsheet)
    wf.run()

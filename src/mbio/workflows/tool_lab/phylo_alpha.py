# -*- coding: utf-8 -*-


from biocluster.workflow import Workflow
import os
import shutil


class PhyloAlphaWorkflow(Workflow):
    """
    alpha系统发育多样性
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(PhyloAlphaWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "update_info", "type": "string"},
            {"name": "main_id", "type": "string"},
            {"name": "table", "type": "infile","format":"tool_lab.simple"},
            {"name": "opt_file", "type": "string","default":"is_tree"},
            {"name": "is_tree", "type": "infile","format":"tool_lab.simple"},
            {"name": "is_seq", "type": "infile","format":"tool_lab.simple"},
            {"name": "index", "type": "string","default":"pd"}

        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())


    def run(self):
        self.tool = self.add_tool("tool_lab.phylo_alpha")
        file_type = self.option('opt_file')
        if file_type == 'is_tree':
            in_file = self.option("is_tree")
        else:
            in_file = self.option("is_seq")

        options = {
            "table": self.option('table'),
            "opt_file": in_file,
            "file_type": file_type,
            "index": self.option('index')
        }
        self.tool.set_options(options)
        self.tool.on('end',self.set_db)

        self.tool.run()

        super(PhyloAlphaWorkflow, self).run()

    def set_db(self):
        """
        导表程序
        """
        if os.path.exists(self.output_dir):
            shutil.rmtree(self.output_dir)
        shutil.copytree(self.tool.output_dir, self.output_dir)

        self.end()

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "系统发育多样性结果", 0]
        ])

        result_dir.add_regexp_rules([
            [".*/.*xls", "xls", "详情表", 0]
        ])

        super(PhyloAlphaWorkflow, self).end()

if __name__ == "__main__":
    from biocluster.wsheet import Sheet
    data = {
        'name': 'tool_lab.phylo_alpha',
        'id': 'tsg_36964',
        'type': 'workflow',
        'options': {
            "main_id" : "5e4ce5e817b2bf4b326f4b20",
            "table" : "/mnt/ilustre/users/sanger-dev/sg-users/zouguanqing/small_tool/phylo_alpha/otu_table.xls",
            #"opt_file" : "/mnt/ilustre/users/sanger-dev/sg-users/zouguanqing/small_tool/phylo_alpha/otu_phylo.tre",
            "is_seq" : "/mnt/ilustre/users/sanger-dev/sg-users/zouguanqing/small_tool/phylo_alpha/otu_reps.fasta",
            "opt_file" : "is_seq",
            "index": "pd"
        }
    }

    wsheet = Sheet(data=data)

    wf = PhyloAlphaWorkflow(wsheet)
    wf.run()

# -*- coding: utf-8 -*-


from biocluster.workflow import Workflow
import os
import shutil


class Tax4funWorkflow(Workflow):
    """
    tax4fun
    """

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(Tax4funWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "update_info", "type": "string"},
            {"name": "table", "type": "infile", "format": "meta.otu.otu_table"},
            {"name": "main_id", "type": "string"}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())

    def deal_infile(self,infile):
        new_file = infile+'_new'
        with open(infile) as fr, open(new_file,'w') as fw:
            fw.write('\t'.join(fr.readline().split('\t')[:-1])+'\n')
            for line in fr:
                sp = line.strip().split('\t')
                new_sp0 = sp[-1]+'; '+sp[0]
                new_line = new_sp0 + '\t' + '\t'.join(sp[1:-1])+'\n'
                fw.write(new_line)
        return new_file


    def run(self):
        new_table = self.deal_infile(self.option('table').path)

        self.tool = self.add_tool("tool_lab.tax4fun")
        options = {
            "in_otu_table": new_table
        }
        self.tool.set_options(options)
        self.tool.on('end',self.set_db)

        self.tool.run()

        super(Tax4funWorkflow, self).run()

    def set_db(self):
        """
        """

        if os.path.exists(self.output_dir):
            shutil.rmtree(self.output_dir)
        shutil.copytree(self.tool.output_dir, self.output_dir)
        if os.path.exists(self.output_dir+'/predictions_Enzyme.xls'):
            os.remove(self.output_dir+'/predictions_Enzyme.xls')

        api = self.api.api('tool_lab.common')
        api.update_main('tax4fun', self.option('main_id'), {"main_id":self.option('main_id')})
        self.end()


    def end(self):

        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "tax4fun结果", 0]
        ])

        result_dir.add_regexp_rules([
            [".*/.*", "", "结果文件", 0],

        ])

        super(Tax4funWorkflow, self).end()

if __name__ == "__main__":
    from biocluster.wsheet import Sheet
    data = {
        'name': 'tool_lab.table2biom',
        'id': 'tsg_36964',
        'type': 'workflow',
        'options': {
            "main_id" : "5e4ce5e817b2bf4b326f4b20",
            "table" : "/mnt/ilustre/users/sanger-dev/sg-users/zouguanqing/small_tool/tax4fun/DATA.1.txt"
            ##/mnt/ilustre/users/sanger-dev/i-sanger_workspace2/20210309/Tax4fun_majorbio_293423_79500_198484
        }
    }

    wsheet = Sheet(data=data)

    wf = Tax4funWorkflow(wsheet)
    wf.run()

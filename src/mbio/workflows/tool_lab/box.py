# -*- coding: utf-8 -*-


from biocluster.workflow import Workflow
import gevent
import pandas as pd



class BoxWorkflow(Workflow):

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(BoxWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "table", "type": "infile", "format": "tool_lab.table"},
            {"name": "group", "type": "infile", "format": "tool_lab.group_table"},
            {"name": "group_checked", "type": "string"},
            {"name": "method", "type": "string", "default": "mean"},  #'sum', 'median','none' ,mean
            {"name": "log_change", "type": "string", "default": "none"}, #none,log2,log10
            {"name": "update_info", "type": "string"},
            {"name": "main_id", "type": "string"}
        ]
        self.add_option(options)
        self.revise_infiles()
        self.set_options(self._sheet.options())


    def run(self):
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_DATA_AFTER_END = False
        self.logger.info("start run workflow")
        gevent.spawn_later(5,self.set_db)
        super(BoxWorkflow, self).run()


    def set_db(self):
        """
        保存结果到mongo数据库中
        """

        api_name = self.api.api("tool_lab.box")
        self.logger.info("正在写入mongo数据库")
        main_id = self.option("main_id")
        if self.option('group').is_set:
            api_name.add_detail(main_id,self.option('table').path,self.option('group').path, self.option('method'),self.option("log_change"))
        else:
            table = pd.read_table(self.option('table').path,sep='\t',index_col=0)
            group = self.work_dir + '/group.xls'
            with open(group,'w') as fw:
                fw.write('#sample\tgroup\n')
                for s in table.columns:
                    fw.write("%s\t%s\n"%(s,s))
            api_name.add_detail(main_id,self.option('table').path,group,'none', self.option("log_change"))
        self.end()

    def end(self):
        super(BoxWorkflow, self).end()


if __name__ == '__main__':
    from biocluster.wsheet import Sheet
    data = {
        'name': 'test',
        'id': 'tsg_36964',
        'type': 'workflow',
        'options': {
            "table" : "/mnt/ilustre/users/sanger-dev/workspace/20200228/Metabolome_tsg_36991/Preprocess/output/pos/metab_abund.txt",
            "group" : "/mnt/ilustre/users/sanger-dev/workspace/20200228/Metabolome_tsg_36991/remote_input/group_table/group.txt",
            "main_id" : "5e9e6a6017b2bf2049a81be3"
        }
    }

    wsheet = Sheet(data=data)

    wf = BoxWorkflow(wsheet)
    wf.run()
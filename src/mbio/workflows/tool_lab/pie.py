# -*- coding: utf-8 -*-


from biocluster.workflow import Workflow
import gevent


class PieWorkflow(Workflow):

    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(PieWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "table", "type": "infile", "format": "tool_lab.simple"},
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
        super(PieWorkflow, self).run()


    def set_db(self):
        """
        保存结果到mongo数据库中
        """

        api_name = self.api.api("tool_lab.pie")
        self.logger.info("正在写入mongo数据库")
        main_id = self.option("main_id")
        api_name.add_detail(main_id,self.option('table').path)
        self.end()

    def end(self):
        super(PieWorkflow, self).end()


if __name__ == '__main__':
    from biocluster.wsheet import Sheet
    data = {
        'name': 'test',
        'id': 'tsg_36964',
        'type': 'workflow',
        'options': {
            "table" : "/mnt/ilustre/users/sanger-dev/workspace/20200228/Metabolome_tsg_36991/Preprocess/output/pos/metab_abund.txt",
            "main_id" : "5e9e6a6017b2bf2049a81be3"
        }
    }

    wsheet = Sheet(data=data)

    wf = PieWorkflow(wsheet)
    wf.run()
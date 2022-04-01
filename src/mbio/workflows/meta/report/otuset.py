# -*- coding: utf-8 -*-
# __author__ = 'xieshichang'
from biocluster.workflow import Workflow
from biocluster.core.exceptions import OptionError


class OtusetWorkflow(Workflow):
    """
    Otu集创建
    """
    def __init__(self, wsheet_object):
        self._sheet = wsheet_object
        super(OtusetWorkflow, self).__init__(wsheet_object)
        options = [
            {"name": "table_name", "type": "string"},
            {"name": "table_id", "type": "string"},
            {"name": "pvalue", "type": "float", "default": 0.0},
            {"name": "qvalue", "type": "float", "default": 0.0},
            {"name": "species_name", "type": "string", "default": ""},
            {"name": "label", "type": "string"},
            {"name": "lda", "type": "float", "default": 0.0},
            {"name": "top", "type": "int", "default": 0},
            {"name": "main_id", "type": "string"},
            {"name": "update_info", "type": "string"},
        ]
        self.add_option(options)
        self.set_options(self._sheet.options())
        self.IMPORT_REPORT_DATA = True
        self.IMPORT_REPORT_AFTER_END = False

    def check_options(self):
        """
        参数二次检查
        """
        if self.option("table_name") == 'sg_otu_venn':
            if not self.option('label'):
                self.set_error('在veen图分析页面必须配置参数 "label"')
        if self.option("table_name") in ["sg_species_difference_lefse"]:
            if (not self.option("pvalue")) and (not self.option("lda")):
                raise OptionError("{},{}不存在请检查".format(self.option("pvalue"), self.option("lda")))
        if self.option("table_name") in ["sg_randomforest"]:
            if (not self.option("top")):
                raise OptionError("{}不存在请检查".format(self.option("top")))
        if self.option("table_name") in ["sg_species_difference_check"]:
            if (not self.option("pvalue")) and (not self.option("qvalue")):
                raise OptionError("{},{},{}不存在请检查".format(self.option("pvalue"), self.option("qvalue"), self.option("species_name")))

    def run_otuset_create(self):
        '''
        在tool进行otu集构建
        '''
        opts = {}
        for k, v in self.get_option_object().items():
            self.logger.info("option name: {}, value: {}".format(k, v))
            if v.value and k != 'update_info':
                opts[k] = v.value
        self.logger.info("otuset_create opts {}".format(opts))
        self.otuset.set_options(opts)
        self.otuset.run()

    def end(self):
        for name, option in self.get_option_object().items():
            self.logger.info("option name: {}, value: {}".format(name, option))
        super(OtusetWorkflow, self).end()

    def run(self):
        """
        运行
        """
        self.logger.info("开始啦")
        self.otuset = self.add_tool('meta.otu.otuset_create')
        self.otuset.on('end', self.end)
        self.run_otuset_create()
        super(OtusetWorkflow, self).run()


if __name__ == "__main__":
    from biocluster.wpm.client import worker_client
    from mainapp.models.mongo.meta import Meta
    import datetime
    from otuset import OtusetWorkflow
    from biocluster.wsheet import Sheet
    import time
    import json
    data = {
        'type': 'workflow',
        'name': 'meta.report.otuset',
        'client': 'client03',
        'rerun': True,
    }

    test_dict = {
        'otuset_test_venn': {
            'table_name': 'sg_otu_venn',
            'table_id': '5f9632bf17b2bf257aabc492',
            'label': 'C & E',
        },
        'otuset_test_diff': {
            'table_name': 'sg_species_difference_check',
            'table_id': '5f521b7617b2bf6ec16fbaa5',
            'pvalue': 1,
            'qvalue': 1
        },
        'otuset_test_lefse': {
            'table_name': 'sg_species_difference_lefse',
            'table_id': '5f521b7717b2bf6ec16fbaaa',
            'lda': 1,
            'pvalue': 1,
        },
        'otuset_test_rf': {
            'table_name': 'sg_randomforest',
            'table_id': '5f97897d17b2bf1abe694530',
            'top': 20,
        }
    }

    meta = Meta()
    for i, opts in test_dict.items():
        table_info = meta.get_main_info(opts['table_id'], opts['table_name'])
        task_info = meta.get_task_info(table_info['task_id'])
        mongo_data = {
            'status': 'start',
            'desc': 'test',
            'name': opts['table_name'] + str(time.time()),
            'project_sn': task_info['project_sn'],
            'task_id': task_info['task_id'],
            'created_ts': datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        main_id = meta.insert_main_table('sg_otuset', mongo_data)
        update_info = {str(main_id): 'sg_otuset'}
        opts['main_id'] = str(main_id)
        opts['update_info'] = json.dumps(update_info)
        data['id'] = i
        data['options'] = opts
        worker = worker_client()
        result = worker.add_task(data)
        print result
        time.sleep(10)

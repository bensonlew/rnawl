# -*- coding: utf-8 -*-
# __author__ = 'guhaidong'

from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from mbio.packages.metabolome.metabsetc import Metabsetc
from mbio.packages.metabolome.common import Relation


class KeggcAgent(Agent):
    """
    代谢集kegg代谢物注释
    last_modify: 2018.6.19
    """

    def __init__(self, parent):
        super(KeggcAgent, self).__init__(parent)
        options = [
            {"name": "metabset", "type": "infile", "format": "metabolome.metabset"},  # 代谢集文件
            {"name": "anno_overview", "type": "infile", "format": "metabolome.overview_table"},  # 代谢集总览表
            {"name": "stat_out", "type": "outfile", "format": "sequence.profile_table"}
        ]
        self.add_option(options)
        self.step.add_steps("keggc")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.keggc.start()
        self.step.update()

    def stepfinish(self):
        self.step.keggc.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 1
        self._memory = '5G'

    def end(self):
        super(KeggcAgent, self).end()


class KeggcTool(Tool):
    def __init__(self, config):
        super(KeggcTool, self).__init__(config)

    def run(self):
        """
        运行
        :return:
        """
        super(KeggcTool, self).run()
        obj = Metabsetc()
        self.logger.info('metabset: '+self.option('metabset').path)
        self.logger.info('anno_overview: '+self.option('anno_overview').path)

        try:
            obj.run(self.option('metabset').path, self.option('anno_overview').path, self.work_dir + '/level.xls')
        except Exception,e:
            self.logger.error(e)
            self.set_error("keggc运行错误: %s", variables=(e), code="34700901")
        self.logger.info("keggc运行完毕")
        self.id_to_name()
        self.set_output()
        self.end()

    def set_output(self):
        self.logger.info("设置结果目录")
        self.option('stat_out').set_path(self.output_dir + '/level.xls')

    def id_to_name(self):
        # 增加一列代谢物名称
        old_file = os.path.join(self.work_dir, 'level.xls')
        self.level_out_final = os.path.join(self.output_dir, 'level.xls')
        table_path = self.option("anno_overview").prop["path"]
        self.logger.info(table_path)
        self.metab_trans = Relation()
        map_table, id_name_dict = self.metab_trans.get_dataframe_and_dict(table_path, metab_name="metab")
        self.metab_trans.add_metabolites_column(id_name_dict, old_file, self.level_out_final)


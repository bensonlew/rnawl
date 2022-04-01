# -*- coding: utf-8 -*-


import os
from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool
from mbio.packages.tool_lab.anno_keggc import AnnoKeggc



class AnnoKeggcAgent(Agent):
    """
    代谢组化合物注释
    """

    def __init__(self, parent):
        super(AnnoKeggcAgent, self).__init__(parent)
        options = [
            {"name": "metab_table", "type": "infile", "format": "tool_lab.simple"},  # 预处理的metab_desc结果
            {"name": "level_out", "type": "outfile", "format": "sequence.profile_table"},  # 注释层级结果
            {"name": "stat_out", "type": "outfile", "format": "sequence.profile_table"},  # 注释统计结果
            {"name": "database_name", "type": "string", "default" : "CBR"},
            {"name": "search_type", "type": "string", "default":""}
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
        检查参数
        :return:
        """
        if not self.option('metab_table'):
            raise OptionError('必须输入代谢物详情', code="34700501")
        return True

    def set_resource(self):
        """
        设置所需资源，需在子类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 1
        self._memory = "5G"

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", ""],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(AnnoKeggcAgent, self).end()


class AnnoKeggcTool(Tool):
    def __init__(self, config):
        super(AnnoKeggcTool, self).__init__(config)
        self.level_out = os.path.join(self.work_dir, 'level.xls')
        self.stat_out = os.path.join(self.work_dir, 'stat.xls')

    def run(self):
        """
        运行
        :return:
        """
        super(AnnoKeggcTool, self).run()
        obj = AnnoKeggc()
        map_dic ={
            "CBR" : "kegg_compound",
            "BP" : "kegg_compound_br08005",
            "EDC": "kegg_compound_br08006",
            "P" : "kegg_compound_br08007" ,
            "PC" : "kegg_compound_br08003",
            "L" : "kegg_compound_br08002"
        }
        self.logger.info(self.option("metab_table").prop["path"])
        try:
            database_name = map_dic[self.option('database_name')]
            obj.run(self.option('metab_table').path, self.level_out, self.stat_out,database_name, self.option("search_type"))
        except Exception as e:
            self.logger.error(e)
            self.set_error(e, code="34700501")
        self.set_output()
        self.end()



    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        for f in ['level.xls', 'stat.xls']:
            out = self.output_dir + '/' + f
            if os.path.exists(out):
                os.remove(out)
            os.link(self.work_dir+'/'+ f, out)
        self.option('level_out').set_path(self.output_dir+'/level.xls')
        self.option('stat_out').set_path(self.output_dir+'/stat.xls')

        self.logger.info("设置anno_keggc分析结果目录成功")

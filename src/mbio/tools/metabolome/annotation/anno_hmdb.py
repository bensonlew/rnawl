# -*- coding: utf-8 -*-
# __author__ = 'shaohua.yuan'

import os, shutil
from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool
import unittest
from mbio.packages.metabolome.common import Relation


class AnnoHmdbAgent(Agent):
    """
    HMDB分析
    version: v1.0
    author: guhaidong
    last_modify: 2018.06.15
    """

    def __init__(self, parent):
        super(AnnoHmdbAgent, self).__init__(parent)
        options = [
            {"name": "anno_overview", "type": "infile", "format": "sequence.profile_table"},
            #{"name": "metabset", "type": "infile", "format": "metabolome.metabset"},
            {"name": "metabset", "type": "infile", "format": "metabolome.metabset,metabolome.metab_abun"},
            {"name": "type", "type": "string", "default":"anno"},
            # name = anno, anno_overview为overview_anno表，工作流注释用，name = annohmdb或metabsetanno，anno_overview为hmdb注释结果表，筛选用
            {"name": "new_overview", "type": "outfile", "format": "metabolome.overview_table"},  # 新总览注释结果
            {"name": "level_out", "type": "outfile", "format": "sequence.profile_table"} , # hmdb注释结果
            {"name": "level_out_origin", "type": "outfile", "format": "sequence.profile_table"}  # 工作流使用所有代谢物的注释结果
        ]
        self.add_option(options)

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option('anno_overview'):
            raise OptionError('must input anno_overview')
        return True

    def set_resource(self):
        """
        设置所需资源，需在子类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 1
        self._memory = "5G"

    def end(self):
        super(AnnoHmdbAgent, self).end()


class AnnoHmdbTool(Tool):
    def __init__(self, config):
        super(AnnoHmdbTool, self).__init__(config)
        self.python_path =  "/miniconda2/bin/python"
        self.script = self.config.PACKAGE_DIR + '/metabolome/anno_hmdb.py'

    def run(self):
        """
        运行
        :return:
        """
        super(AnnoHmdbTool, self).run()
        self.logger.info("开始运行命令！")
        self.start_anno()
        self.id_to_name()
        self.set_output()

    def start_anno(self):
        self.logger.info(self.option("anno_overview"))
        anno_overview = self.option("anno_overview").prop["path"]
        self.logger.info(anno_overview)
        cmd = self.python_path + ' {} -i {} -o {} '.format(self.script, anno_overview, self.output_dir)
        if self.option("type") == "anno":
            cmd += " --origin "
        if self.option("type") == "annohmdb":
            cmd += " --header "
        if self.option("metabset"):
            metabset = self.option("metabset").prop["path"]
            cmd += " -s {} ".format(metabset)
        if self.config.DBVersion:
            cmd += ' -mongo ' + str(self.config.DBVersion)
        self.logger.info(cmd)
        command = self.add_command('hmdb_anno', cmd, ignore_error=True).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("hmdb_anno succeed")
        else:
            self.set_error("hmdb_anno failed")
            raise Exception("hmdb_anno failed")

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        self.metab_trans = Relation()
        anno_overview = self.option("anno_overview").prop["path"]
        if self.option("metabset") and self.option("type") != "anno" :
            if os.path.exists(self.output_dir + "/HmdbLevel.xls"):
                self.metab_trans.add_oid_fun(self.output_dir + "/HmdbLevel.xls", anno_overview, link_k='Metab_id',oldfile_link_id=1)
        else:
            if os.path.exists(self.output_dir + "/HmdbLevel.xls"):
                self.metab_trans.add_oid_fun(self.output_dir + "/HmdbLevel.xls", anno_overview, link_k='metab_id',oldfile_link_id=1)  #20191617
        if os.path.exists(self.output_dir + "/HmdbLevel.xls"):
            self.option('level_out').set_path(self.output_dir + "/HmdbLevel.xls")

        if self.option("type") == "anno":
            self.metab_trans.add_oid_fun(self.output_dir + "/HmdbLevel_Origin.xls", anno_overview, link_k='metab_id',oldfile_link_id=1)  #20191617
            self.option('new_overview').set_path(self.output_dir + "/anno.xls")
            self.option('level_out_origin').set_path(self.output_dir + "/HmdbLevel_Origin.xls")
        self.logger.info("设置hmdb_anno分析结果目录成功")
        self.end()

    def id_to_name(self):
        # 增加一列代谢物名称, old_path为移动位置，new_path为最终位置
        old_files = ["HmdbClass.xls","HmdbSubclass.xls","HmdbSuperclass.xls"]
        table_path = self.option("anno_overview").prop["path"]
        self.logger.info(table_path)
        self.metab_trans = Relation()
        if self.option("type") == "anno":
            map_table, id_name_dict = self.metab_trans.get_dataframe_and_dict(table_path, metab_name="metab")
        else:
            map_table, id_name_dict = self.metab_trans.get_dataframe_and_dict(table_path, id_name ="Metab_id")
        for eachfile in old_files:
            old_path = os.path.join(self.work_dir, eachfile)
            new_path = os.path.join(self.output_dir, eachfile)
            if os.path.exists(old_path):
                os.remove(old_path)
            if os.path.exists(new_path):
                shutil.move(new_path, old_path)
                self.metab_trans.add_metabolites_column(id_name_dict, old_path, new_path, id_name="Metab_ids")




class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            #"id": "CreatTable" + str(random.randint(1, 10000)),
            "id": "AnnoHMDB",
            "type": "tool",
            "name": "metabolome.annotation.anno_hmdb",
            "instant": True,
            "options": dict(
                type = "anno",
                anno_overview="/mnt/ilustre/users/sanger-dev/workspace/20190117/Metabolome_tsg_33258/Anno/AnnoOverview/output/anno.xls",
                #metabset="/mnt/ilustre/users/sanger-dev/workspace/20190109/AnnoHmdb_tsg_33213_99153_966940/remote_input/metabtable/metab_abund.txt",
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
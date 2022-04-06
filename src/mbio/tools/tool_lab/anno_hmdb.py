# -*- coding: utf-8 -*-


import os, shutil
from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool
import unittest
from mbio.packages.metabolome.common import Relation


class AnnoHmdbAgent(Agent):


    def __init__(self, parent):
        super(AnnoHmdbAgent, self).__init__(parent)
        options = [
            {"name": "metab_name", "type": "string"},
            {"name": "out_table", "type": "outfile", "format": "sequence.profile_table"}
        ]
        self.add_option(options)

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option('metab_name'):
            raise OptionError('must input metab_name')
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
        self.script = self.config.PACKAGE_DIR + '/tool_lab/anno_hmdb.py'

    def run(self):
        """
        运行
        :return:
        """
        super(AnnoHmdbTool, self).run()
        self.logger.info("开始运行命令！")
        self.start_anno()

        self.set_output()

    def start_anno(self):
        tmp = self.work_dir+'/tmp.xls'
        with open(tmp,'w') as fw:
            fw.write('metab_id\tname\n')
            id = 0
            for m in self.option('metab_name').split('\n'):
                if m.strip() == '':
                    continue
                id+=1
                fw.write('metab_'+str(id)+'\t'+m+'\n')

        cmd = self.python_path + ' {} -i {} -o {} '.format(self.script, tmp, self.output_dir)

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

        if os.path.exists(self.output_dir + "/Hmdb_detail.xls"):
            self.option('out_table').set_path(self.output_dir + "/Hmdb_detail.xls")


        self.logger.info("设置hmdb_anno分析结果目录成功")
        self.end()




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
            "name": "tool_lab.anno_hmdb",
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
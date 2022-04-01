# !usr/bin/python
# -*- coding: utf-8 -*-
# __author__ = "XueQinwen"

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import shutil
import unittest
import random

class RelationshipSankeyAgent(Agent):
    def __init__(self, parent):
        super(RelationshipSankeyAgent, self).__init__(parent)
        options=[
            {"name":"input_table", "type": "infile", "format": "small_rna.common"},
            ]
        self.add_option(options)

    def check_options(self):
        if not self.option("input_table"):
            raise OptionError("请输入input表格")

        
    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 2
        self._memory = "5G"
    
    def end(self):
        super(RelationshipSankeyAgent, self).end()

class RelationshipSankeyTool(Tool):
    def __init__(self, config):
        super(RelationshipSankeyTool,self).__init__(config)
        self.input_table = self.option("input_table").prop['path']
    
    def run(self):
        """
        运行
        """
        self.get_otu_info()
        self.write_table()
        self.set_output()
        self.end()

    def get_otu_info(self):
        otu_info = {}
        self.sankey_infos = []
        with open(self.input_table,"r") as otu:
            line=otu.readline()
            while 1:
                line = otu.readline()
                if not line:
                    break
                fd = line.rstrip().split('\t')
                sankey_info = {
                    "source":fd[1],
                    "target":fd[2],
                    "value":1,
                }
                self.sankey_infos.append(sankey_info)
                if otu_info.has_key(fd[0]):
                    if otu_info[fd[0]].has_key(fd[1]):
                        otu_info[fd[0]][fd[1]] += 1
                    else:
                        otu_info[fd[0]][fd[1]] = 1
                else:
                    otu_info[fd[0]]={fd[1]:1}
        for key in sorted(otu_info.keys()):
            for key2 in sorted(otu_info[key].keys()):
                sankey_info= {}
                sankey_info['source']=key
                sankey_info['target']=key2
                if otu_info[key][key2] == 0:
                    continue
                sankey_info['value']=otu_info[key][key2]
                self.sankey_infos.append(sankey_info)

    
    def write_table(self):
        self.logger.info(self.sankey_infos)
        self.table_path = os.path.join(self.work_dir, "plot.txt")
        with open(self.table_path,"w") as tb:
            for i in self.sankey_infos:
                tb.write("{}\t{}\t{}\n".format(i['source'],i['target'],i['value']))

    def set_output(self):
        '''
        将结果文件赋值到output文件夹下面
        :return:
        '''
        self.logger.info("设置结果目录")
        if len(self.output_dir) > 0:
            shutil.rmtree(self.output_dir)
        os.mkdir(self.output_dir)
        try:
            os.link(self.table_path,
                        os.path.join(self.output_dir,"plot.txt"))
        except Exception as e:
            self.set_error("设置结果目录失败{}".format(e))
        
class TestFunction(unittest.TestCase):
    """
    测试脚本用
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            'id': 'sankey_relationship' + str(random.randint(1, 100000)),
            "type": "tool",
            "name": "tool_lab.relationship_sankey",
            "options": {
                 "input_table":"/mnt/ilustre/users/sanger-dev/sg-users/xueqinwen/tool_lab_tools/relationship_sankey/otu.txt",
                }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
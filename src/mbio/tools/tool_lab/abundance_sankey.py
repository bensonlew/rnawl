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

class AbundanceSankeyAgent(Agent):
    def __init__(self, parent):
        super(AbundanceSankeyAgent, self).__init__(parent)
        options=[
            {"name":"otu_table", "type": "infile", "format": "meta.otu.otu_table"},
            {"name":"group_table","type": "infile", "format": 'toolapps.group_table,meta.otu.group_table'},
            ]
        self.add_option(options)

    def check_options(self):
        if not self.option("otu_table"):
            raise OptionError("请输入otu表格")
        if not self.option("group_table"):
            raise OptionError("请输入分组表格")
        
    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 2
        self._memory = "5G"
    
    def end(self):
        super(AbundanceSankeyAgent, self).end()

class AbundanceSankeyTool(Tool):
    def __init__(self, config):
        super(AbundanceSankeyTool,self).__init__(config)
        self.otu_table = self.option("otu_table").prop['path']
        self.group_table = self.option("group_table").prop['path']
    
    def run(self):
        """
        运行
        """
        self.get_group_info()
        self.get_otu_info()
        self.write_table()
        self.set_output()
        self.end()

    def get_otu_info(self):
        self.otu_info = {}
        with open(self.otu_table,"r") as otu:
            line=otu.readline()
            title = line.rstrip().split('\t')
            group = title[1:]
            while 1:
                line = otu.readline()
                if not line:
                    break
                fd = line.rstrip().split('\t')
                whole_name = fd[0]
                names = whole_name.split(';')
                last_name = whole_name.split(';')[-1]
                for name in names:
                    if self.otu_info.has_key(name):
                        continue
                    else:
                        self.otu_info[name]={}
                info = fd[1:]
                sum_info = 0
                for i in range(len(info)):
                    if not self.otu_info[last_name].has_key(self.group_info[group[i]]):
                        self.otu_info[last_name][self.group_info[group[i]]] = float(info[i])
                    else:
                        self.otu_info[last_name][self.group_info[group[i]]] += float(info[i])
                    sum_info+=float(info[i])
                for i in range(len(names)):
                    if i == len(names)-1:
                        break
                    if not self.otu_info[names[i]].has_key(names[i+1]):
                        self.otu_info[names[i]][names[i+1]]=sum_info
                    else:
                        self.otu_info[names[i]][names[i+1]]+=sum_info
        self.sankey_infos = []
        for key in sorted(self.otu_info.keys()):
            for key2 in sorted(self.otu_info[key].keys()):
                sankey_info= {}
                sankey_info['source']=key
                sankey_info['target']=key2
                if self.otu_info[key][key2] == 0:
                    continue
                sankey_info['value']=self.otu_info[key][key2]
                self.sankey_infos.append(sankey_info)

    
    def get_group_info(self):
        self.group_info = {}
        with open(self.group_table,'r') as group_table:
            line = group_table.readline()
            while 1:
                line = group_table.readline()
                if not line:
                    break
                fd = line.rstrip().split('\t')
                sample=fd[0]
                group = fd[1]
                # if self.group_info.has_key(group):
                #     self.group_info[group].apend(sample)
                # else:
                #     self.group_info[group] = [sample]
                if self.group_info.has_key(sample):
                    self.logger.info("分组表有重复的样本")
                else:
                    self.group_info[sample]=group

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
            'id': 'sankey_abundance' + str(random.randint(1, 100000)),
            "type": "tool",
            "name": "tool_lab.abundance_sankey",
            "options": {
                 "otu_table":"/mnt/ilustre/users/sanger-dev/sg-users/xueqinwen/tool_lab_tools/abudance_sankey/otu.txt",
                 "group_table":"/mnt/ilustre/users/sanger-dev/sg-users/xueqinwen/tool_lab_tools/abudance_sankey/group.txt",
                }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
# -*- coding: utf-8 -*-
# __author__ = 'ysh'
# last_modify:20190409

from biocluster.module import Module
import os
import shutil
from biocluster.core.exceptions import OptionError
import unittest


class SignalpModule(Module):
    def __init__(self, work_id):
        super(SignalpModule, self).__init__(work_id)
        options = [
            {"name": "query", "type": "string"}, # 多个时逗号分隔
            {"name": "type", "type": "string", "default": "bac"},  # 菌株类型gram-，gram+，bac, euk
            {"name": "sample", "type": "string"},  # 样品名称, 多个时逗号分隔
            {"name": "signalp", "type": "outfile", "format": "sequence.profile_table"},
            {"name": "signalp1", "type": "outfile", "format": "sequence.profile_table"},
        ]
        self.add_option(options)
        self.tool_list = []

    def check_options(self):
        if not self.option("query"):
            raise OptionError("必须设置输入序列文件")
        if not self.option("sample"):
            raise OptionError("必须设置样品名称")
        if self.option("type") not in ["gram-", "gram+", "euk", "bac"]:
            raise OptionError("菌株选择类型错误，只能为gram-，gram+，euk, bac")
        return True

    def run_singalp(self):
        self.spe_list = self.option('sample').split(",")
        self.query = self.option("query").split(",")
        self.logger.info(self.query)
        self.sig_dir = {}
        n = 0
        if self.option("type") == "bac":
            self.type_list = ["gram-","gram+"]
        else:
            self.type_list = [self.option("type")]
        for i in self.spe_list:
            for j in self.type_list:
                self.sig = self.add_tool("align.signalp_anno")
                options = {
                    'query': self.query[n],
                    'type': j,
                    "out_format": "short",
                    "sample": i,
                    "d_score": True
                }
                self.sig.set_options(options)
                self.sig_dir[i + j] = self.sig.output_dir
                self.tool_list.append(self.sig)
            n = n + 1
        if len(self.tool_list) > 1:
            self.on_rely(self.tool_list, self.set_output)
            self.logger.info(self.tool_list)
        else:
            self.tool_list[0].on('end', self.set_output)
        for tool in self.tool_list:
            tool.run()

    def set_output(self):
        for each in self.spe_list:
            for n in self.type_list:
                sample_dir = os.path.join(self.output_dir, each)
                file_dir = self.sig_dir[each + n]
                self.logger.info(file_dir)
                self.linkdir(file_dir)
        if os.path.exists(self.output_dir + "/" + self.option("sample") + "_Gram+_SignalP.txt"):
            self.option("signalp", self.output_dir + "/" + self.option("sample") + "_Gram+_SignalP.txt")
        if os.path.exists(self.output_dir + "/" + self.option("sample") + "_Gram-_SignalP.txt"):
            self.option("signalp1", self.output_dir + "/" + self.option("sample") + "_Gram-_SignalP.txt")
        self.end()

    def run(self):
        super(SignalpModule, self).run()
        self.run_singalp()

    def end(self):
        super(SignalpModule, self).end()

    def linkdir(self, dirpath):
        """
        link一个文件夹下的所有文件到本module的output目录
        """
        allfiles = os.listdir(dirpath)
        newdir = self.output_dir
        if not os.path.exists(newdir):
            os.mkdir(newdir)
        oldfiles = [os.path.join(dirpath, i) for i in allfiles]
        newfiles = [os.path.join(newdir, i) for i in allfiles]
        for newfile in newfiles:
            if os.path.exists(newfile):
                os.remove(newfile)
        for i in range(len(allfiles)):
            os.link(oldfiles[i], newfiles[i])

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
            "id": "Singalp_module",
            "type": "module",
            "name": "bacgenome.signalp",
            "instant": True,
            "options": dict(
                query="s3://bacgenome/files/m_188/188_5c6cde8a58bab/tsg_33557/workflow_results//B11/assembly_predict/predict/CDS_predict/B11_CDS.fnn",
                samples="B11",

            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
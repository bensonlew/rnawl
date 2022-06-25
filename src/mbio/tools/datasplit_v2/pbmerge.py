# -*_ coding: utf-8 -*-
# __author__ = 'XueQinwen'
# last_modified: 20210907
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import re
import os
import unittest

class PbmergeAgent(Agent):
    """
    SMRT Link v10.1
    从subreads里提取ccs数据
    """
    def __init__(self, parent):
        super(PbmergeAgent,self).__init__(parent)
        options = [
            {"name":"ccs_bams","type":"string"},
            # {"name":"barcode_path","type":"string"},
            # {"name":"index_path","type":"string"}
        ]
        self.add_option(options)
        self.queue = "chaifen"  # 投递到指定的队列chaifen

    def set_resource(self):
        '''
        设置所需资源
        '''
        self._cpu = 20 
        self._memory = "100G"
    
    def end(self):
        super(PbmergeAgent,self).end()

class PbmergeTool(Tool):
    def __init__(self, config):
        super(PbmergeTool,self).__init__(config)
        self.pbmerge = "bioinfo/ref_rna_v3/HTSeq/miniconda3/bin/pbmerge"
        # self.lima = "program/SmrtLink/smrtlink/smrtcmds/bin/lima"
        self.samtools = "miniconda2/bin/samtools"

    def run(self):
        super(PbmergeTool, self).run()
        self.run_pbmerge()
        self.end()

    def run_pbmerge(self):
        self.ccs_new = os.path.join(self.work_dir,"ccs.new.bam")
        cmd = '{} -o {} '.format(self.pbmerge, self.ccs_new)
        for i in os.listdir(self.option('ccs_bams')):
            if i.endswith(".bam"):
                cmd += " {}/{}".format(self.option('ccs_bams'), i)
        self.logger.info(cmd)
        self.logger.info("开始合并ccs")
        command = self.add_command("pbmerge",cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("合并ccs完成")
        else:
            self.set_error("合并数据出现了问题")
    
    # def run_lima(self):
    #     self.split_dir = os.path.join(self.work_dir,"split")
    #     os.mkdir(self.split_dir)
    #     output_prefix = os.path.join(self.split_dir,"split.bam")
    #     cmd = "{} --different --split-bam-named --min-passes 1 --peek-guesss --num-threads 20 ".format(self.lima)
    #     cmd += "{} {} {}".format(self.ccs_new,self.option("barcode_path"),output_prefix)
    #     command = self.add_command("lima",cmd).run()
    #     self.wait(command)
    #     if command.return_code == 0:
    #         self.logger.info("拆分完成")
    #     else:
    #         self.set_error("拆分出现问题")

    # def get_index(self):
    #     self.index = {}
    #     with open(self.option("index_path"),'r') as ip:
    #         while 1:
    #             line = ip.readline()
    #             if not line:
    #                 break 
    #             fd = line.rstrip().split('\t')
    #             self.index[fd[0]]={}
    #             self.index[fd[0]][1]= fd[1]
    #             self.index[fd[0]][2]= fd[2]
    
    def set_output(self):
        # self.o_split_dir = os.path.join(self.output_dir,'split')
        # if not os.path.exists(self.o_split_dir):
        #     os.mkdir(self.o_split_dir)
        # for i in self.index.keys():
        #     try:
        #         os.link(os.path.join(self.split_dir,'split.{}--{}.bam'.format(self.index[i][1],self.index[i][2])),
        #             os.path.join(self.o_split_dir,'split.{}.bam'.format(i)))
        #         os.link(os.path.join(self.split_dir,'split.{}--{}.bam.pbi'.format(self.index[i][1],self.index[i][2])),
        #             os.path.join(self.o_split_dir,'split.{}.bam.pbi'.format(i)))
        #         os.link(os.path.join(self.split_dir,'split.{}--{}.subreadset.xml'.format(self.index[i][1],self.index[i][2])),
        #             os.path.join(self.o_split_dir,'split.{}.subreadset.xml'.format(i)))
        #     except:
                # raise OptionError("设置bam输出目录失败")
        try:
            os.link(os.path.join(self.ccs_new),
                os.path.join(self.output_dir,"ccs.new.bam"))
        except:
            raise OptionError("设置ccs.new.bam到结果目录失败")




class TestFunctionf(unittest.TestCase):
    '''
    测试脚本
    '''

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "pbmerge_" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "datasplit_v2.pbmerge",
            "options": {
                "ccs_bams":"/mnt/ilustre/users/sanger-dev/workspace/20210915/PacbioSplit_20181016PE300_20210915_083704/ccs_bam"
            }
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()

    
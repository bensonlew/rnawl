# -*- coding: utf-8 -*-

from biocluster.agent import Agent
from biocluster.tool import Tool
import os,re
import subprocess
from biocluster.core.exceptions import OptionError
import types
import pandas as pd




class FramebotAgent(Agent):
    def __init__(self, parent):
        super(FramebotAgent, self).__init__(parent)
        options = [
            {"name": "fasta", "type": "infile","format":"sequence.fasta"},
            {"name": "database", "type": "string"},  # 数据库选择
            {"name": "acid_length", "type": "int", "default": 80},  # 氨基酸长度阈值
            {"name": "seq_identity", "type": "float", "default": 0.4},
            {"name": "ref_acid", "type": "infile", "format": "sequence.fasta"},  # 参考氨基酸fasta文件
            {"name": "nucl_fasta", "type": "infile", "format": "sequence.fasta"},
            {"name": "prot_fasta", "type": "infile", "format": "sequence.fasta"},
            {'name': 'framebot_nucl_fasta', 'type': 'outfile', 'format': 'sequence.fasta'},  # 输出核酸序列
            {'name': 'framebot_prot_fasta', 'type': 'outfile', 'format': 'sequence.fasta'},  # 输出核酸序列
        ]
        self.add_option(options)
        self.step.add_steps('framebot')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.framebot.start()
        self.step.update()

    def step_end(self):
        self.step.framebot.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检查
        """
        if not self.option('fasta').is_set:
            raise OptionError('必须提供输入fasta序列')
        if self.option("database") == "custom_mode":
            if not self.option('ref_acid').is_set:
                raise OptionError('必须提供输入核酸序列')
        else:
            if self.option("database") not in ['fgr/amoA_archaea_202012', 'fgr/amoA_bacteria_202012',
                                           'fgr/amoA_AOB_like_202012', 'fgr/amoA_comammox_202012', 'fgr/nosZ_202012',
                                           'fgr/nosZ_atypical_1_202012', 'fgr/nosZ_atypical_2_202012', 'fgr/nirK_202012',
                                           'fgr/nirS_202012', 'fgr/mcrA_202012', 'fgr/nifH_202012', 'fgr/pmoA_202012',
                                           'fgr/mmoX_202012']:
                raise OptionError("数据库%s不被支持", variables=(self.option("database")), code="12700104")


    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 16
        self._memory = '150G'

    def end(self):
        super(FramebotAgent, self).end()


class FramebotTool(Tool):
    def __init__(self, config):
        super(FramebotTool, self).__init__(config)
        self._version = '1.0'
        self.framebor_path = self.config.SOFTWARE_DIR + '/bioinfo/meta/RDPTools-2.0.2/RDPTools/'
        self.java_path =  "/program/sun_jdk1.8.0/bin/java"
        self.python_path = "miniconda2/bin/python"
        self.gcc = self.config.SOFTWARE_DIR + '/gcc/5.1.0/bin'
        self.gcc_lib = self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64'
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)

    def run_framebot(self):

        if self.option("database") == "custom_mode":

            cmd = "{} -Xmx10240M -Xms10240M -jar {}FrameBot.jar framebot -o {} -l {} -i {} -N {} {}".format(self.java_path, self.framebor_path,
                                                                             "framebot", self.option("acid_length")
                                                                             , self.option("seq_identity"),
                                                                             self.option("ref_acid").prop["path"],
                                                                             self.option("fasta").prop["path"])
        else:
            acid_path = self.config.SOFTWARE_DIR + "/database/Framebot/fgr/seed/" + self.option("database") + ".seeds"
            cmd = "{} -Xmx10240M -Xms10240M -jar {}FrameBot.jar framebot -o {} -l {} -i {} -N {} {}".format(self.java_path, self.framebor_path,
                                                                             "framebot", self.option("acid_length")
                                                                             , self.option("seq_identity"),
                                                                             acid_path,
                                                                             self.option("fasta").prop["path"])
        self.logger.info(cmd)
        command = self.add_command('run_framebot', cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_framebot成功")
        else:
            self.set_error("run_framebot失败")

    def run(self):
        super(FramebotTool, self).run()
        self.run_framebot()
        self.set_output()

    def set_output(self):
        if os.path.getsize(self.work_dir + '/framebot_corr_nucl.fasta') == 0:
            self.set_error("FrameBot矫正后序列数为0")
        with open(self.work_dir + '/framebot_corr_nucl.fasta') as f, open(self.work_dir + '/framebot_corr_nucl_new.fasta',"w") as t:
            data = f.readlines()
            for i in data:
                if i.startswith(">"):
                    t.write(i)
                else:
                    new = i.strip().upper()
                    t.write(new+"\n")
        nucl_fasta = os.path.join(self.output_dir, "framebot_nucl.fasta")
        prot_fasta = os.path.join(self.output_dir, "framebot_prot.fasta")
        for i in nucl_fasta,prot_fasta:
            if os.path.exists(i):
               os.remove(i)
        os.link(self.work_dir + '/framebot_corr_nucl_new.fasta', nucl_fasta)
        os.link(self.work_dir + '/framebot_corr_prot.fasta', prot_fasta)
        self.option('framebot_nucl_fasta',nucl_fasta)
        self.option('framebot_prot_fasta', prot_fasta)
        self.end()
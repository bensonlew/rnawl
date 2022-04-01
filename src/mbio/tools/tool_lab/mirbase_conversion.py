# -*- coding: utf-8 -*-
# __author__ = 'zhangyitong'

import unittest
import os
import glob
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError


class MirbaseConversionAgent(Agent):
    def __init__(self, parent):
        super(MirbaseConversionAgent, self).__init__(parent)
        options = [
            {"name": "infile", "type": "infile", 'format': 'medical_transcriptome.common'},
            {"name": "instr", "type": "string", 'default': ''},
            {"name": "target_ver", "type": "string", "default": "v22"},
        ]
        self.add_option(options)
        self.step.add_steps("mirbase_conversion")
        self.on('start', self.step_start)
        self.on('end', self.step_finish)

    def step_start(self):
        self.step.mirbase_conversion.start()
        self.step.update()

    def step_finish(self):
        self.step.mirbase_conversion.finish()
        self.step.update()

    def check_options(self):
        if not self.option("infile").is_set and not self.option("instr"):
            raise OptionError("必须设置输入含有miRNA name/Accessions 的文件或字符串")
        if not self.option("target_ver"):
            raise OptionError("必须设置输入目标miRBase版本")
        return True

    def set_resource(self):
        self._cpu = 1
        self._memory = "10G"

    def end(self):
        # result_dir = self.add_upload_dir(self.output_dir)
        # result_dir.add_relpath_rules([
        #     [".", "", "结果输出目录"]
        # ])
        # result_dir.add_regexp_rules([
        #     [r"disgenet_enrichment.xls$", "xls", "DisGeNET富集分析结果"]
        # ])
        super(MirbaseConversionAgent, self).end()


class MirbaseConversionTool(Tool):
    def __init__(self, config):
        super(MirbaseConversionTool, self).__init__(config)
        self._version = "v1.0"
        software_dir = self.config.SOFTWARE_DIR
        self.gcc = software_dir + '/gcc/5.1.0/bin'
        self.gcc_lib = software_dir + '/gcc/5.1.0/lib64'
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)
        self.mirbase_conversion_path = self.config.PACKAGE_DIR + "/tool_lab/mirbase_conversion.r"
        # self.r_path = self.config.SOFTWARE_DIR + "/bioinfo/ref_rna_v3/miniconda3/bin/Rscript"
        self.r_path = self.config.SOFTWARE_DIR + "/program/R-3.5.1/bin/Rscript"

    def run(self):
        super(MirbaseConversionTool, self).run()
        self.run_mirbase_conversion()
        self.set_output()
        self.end()

    def run_mirbase_conversion(self):
        cmd = '{} {}'.format(self.r_path, self.mirbase_conversion_path)
        if self.option('infile').is_set:
            cmd += ' -f {}'.format(self.option('infile').prop["path"])
        elif self.option('instr'):
            cmd += ' -s {}'.format(self.option('instr'))
        cmd += ' -t {}'.format(self.option('target_ver'))
        cmd_name = 'convert_mirbase_id'
        command = self.add_command(cmd_name, cmd, shell=True)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("{} 运行成功".format(cmd_name))
        elif command.return_code is None:
            self.logger.warn("运行{}出错，返回值为None，尝试重新运行".format(cmd_name))
            command.rerun()
            self.wait()
            if command.return_code is 0:
                self.logger.info("{} 运行成功".format(cmd_name))
            else:
                self.set_error("运行%s>>>%s出错", variables=(cmd_name, cmd))
        else:
            self.set_error("运行%s>>>%s出错", variables=(cmd_name, cmd))

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        files = glob.glob(self.work_dir + '/*.xls')
        for each in files:
            name = os.path.basename(each)
            link = os.path.join(self.output_dir, name)
            if os.path.exists(link):
                os.remove(link)
            os.link(each, link)
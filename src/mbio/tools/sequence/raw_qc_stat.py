# -*- coding: utf-8 -*-
# __author__ :zhouxuan

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class RawQcStatAgent(Agent):
    """
    PE质控fastq序列质控前后信息统计
    version 1.0
    author: zhouxuan
    last_modify: 20170526
    """

    def __init__(self, parent):
        super(RawQcStatAgent, self).__init__(parent)
        options = [
            {"name": "base_info_dir", "type": "infile", "format": "sequence.baif_dir"},  # 碱基质量统计表的文件夹
            {"name": "sickle_dir", "type": "infile", 'format': "sequence.fastq_dir"},  # 经过质控的fastq文件的文件夹
        ]
        self.add_option(options)
        self.step.add_steps("sickle_stat")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.sickle_stat.start()
        self.step.update()

    def stepfinish(self):
        self.step.sickle_stat.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option('base_info_dir'):
            raise OptionError("必须输入base_info_dir文件", code="34002301")
        if not self.option('sickle_dir'):
            raise OptionError("必须输入sickle_dir文件", code="34002302")
        return True

    def set_resource(self):  # 后续需要测试确认
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        """
        self._cpu = 3
        self._memory = '5G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"]
        ])
        result_dir.add_regexp_rules([
            ["", "", ""],
        ])
        super(RawQcStatAgent, self).end()


class RawQcStatTool(Tool):
    def __init__(self, config):
        super(RawQcStatTool, self).__init__(config)
        self._version = "v1.0"
        #self.script_path = "bioinfo/seq/scripts/readStat.pl"
        self.perl_path = '/program/perl-5.24.0/bin/perl'
        self.script_path = self.config.PACKAGE_DIR + "/sequence/scripts/readStat.pl"

    def run(self):
        super(RawQcStatTool, self).run()
        self.read_stat()
        self.set_output()
        self.end()

    def read_stat(self):
        stat_list = os.path.join(self.work_dir, 'stat_list')
        sickle_list = os.path.join(self.work_dir, 'sickle_list')
        self.get_list(stat_list, self.option('base_info_dir').prop['path'])
        self.get_list(sickle_list, self.option('sickle_dir').prop['path'])
        cmd = "{} {} {} {} {}".format(self.perl_path, self.script_path, stat_list, sickle_list, self.output_dir + "/reads")
        self.logger.info("start read_stat")
        command = self.add_command("read_stat", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("read_stat done")
        else:
            command.rerun()
            self.wait(command)
            if command.return_code == 0:
                self.logger.info("read_stat done")
            else:
                self.set_error("read_stat error", code="34002301")
                raise Exception("read_stat error")

    def get_list(self, list_name, dir_path):
        stat_name = os.listdir(dir_path)
        with open(list_name, 'a') as w:
            for name in stat_name:
                sickle_path = os.path.join(dir_path, name)
                w.write(sickle_path + "\n")

    def set_output(self):
        if len(os.listdir(self.output_dir)) == 2:
            self.logger.info("read_stat比对的结果文件正确生成")
        else:
            self.set_error("hmm比对结果出错", code="34002302")

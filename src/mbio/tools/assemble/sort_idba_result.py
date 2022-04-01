# -*- coding: utf-8 -*-
# __author__ = 'guhaidong'

import os
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool


class SortIdbaResultAgent(Agent):
    """
    进行newbler拼接
    version: v1.0
    author: guhaidong
    last_modify: 2017.09.12
    """

    def __init__(self, parent):
        super(SortIdbaResultAgent, self).__init__(parent)
        options = [
            {"name": "idba_contig", "type": "infile", "format": "sequence.fasta_dir"},  # 输入idba拼接结果路径
            {"name": "newbler", "type": "string"},  # 输入newbler拼接结果输出路径
            {"name": "min_contig", "type": "int", "default": 300},  # 输入最短contig长度，默认300
            # {"name": "predict", "type": "outfile", "format": "sequence.fasta_dir"}, # 输出fasta结果文件夹(供基因预测用)
        ]
        self.add_option(options)
        self.step.add_steps("SortIdbaResult")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.SortIdbaResult.start()
        self.step.update()

    def stepfinish(self):
        self.step.SortIdbaResult.finish()
        self.step.update()

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option('idba_contig'):
            raise OptionError('必须输入idba拼接结果路径', code="31301801")
        self.logger.info(self.option('newbler'))
        if not self.option('newbler'):
            raise OptionError('必须输入newbler拼接结果状态文件路径', code="31301802")
        return True

    def set_resource(self):
        """
        设置所需资源，需在子类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 1
        self._memory = "10G"  # 改回 by GHD @ 20180428
        # tmp_mem = 10 * (self._rerun_time + 1)  # 每次重运行的内存增加10G by GHD @ 20180320
        # self._memory = '%sG' % tmp_mem
        # self.logger.info('sort_idba_result use memory : ' + self._memory)

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(SortIdbaResultAgent, self).end()


class SortIdbaResultTool(Tool):
    def __init__(self, config):
        super(SortIdbaResultTool, self).__init__(config)
        self.perl_path = '/program/perl/perls/perl-5.24.0/bin/perl '
        # self.sort_result_path = self.config.SOFTWARE_DIR + '/bioinfo/metaGenomic/scripts/sort_idba_result.pl '
        self.sort_result_path = self.config.PACKAGE_DIR + '/metagenomic/scripts/sort_idba_result.pl '

    def sort_run(self):
        """
        对混拼结果进行整理
        :return:
        """
        cmd = self.perl_path + self.sort_result_path + ' %s   %s   %s   %s   %s' \
                                  % (self.option('idba_contig').prop['path'] + '/more.list',
                                     self.option('idba_contig').prop['path'] + '/less.list',
                                     self.option('min_contig'),
                                     self.option('newbler') + '/454ReadStatus.txt',
                                     self.work_dir)
        self.logger.info("开始整合混拼结果")
        command = self.add_command("sort_idba_result", cmd, ignore_error=True)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("整合混拼结果完成")
        elif command.return_code == 137:
            self.add_state('memory_limit', 'memory is low!')   # add memory limit error by guhaidong @ 20180320
        else:
            self.set_error("整合混拼结果出错!", code="31301801")

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        out_predict_dir = self.output_dir + '/predict'
        newbler_contig = self.output_dir + '/predict/newbler.contig.fa'
        newbler_mix = self.output_dir + '/Newbler_Mix.contig.fa'
        if os.path.exists(out_predict_dir):
            shutil.rmtree(out_predict_dir)
        os.mkdir(out_predict_dir)
        if os.path.exists(newbler_mix):
            os.remove(newbler_mix)
        newbler_output_dir = self.option('newbler')
        os.link(newbler_output_dir + '/newbler.contig.fa', newbler_contig)
        for file in os.listdir(self.work_dir):
            if "Newbler_Mix" in file:
                os.link(self.work_dir + '/' + file, self.output_dir + '/' + file)
            elif "contig.fa" in file:
                os.link(self.work_dir + '/' + file, out_predict_dir + '/' + file)
        # self.option('predict').set_path(out_predict_dir)
        self.logger.info("设置结果目录成功")

    def run(self):
        """
        运行
        :return:
        """
        super(SortIdbaResultTool, self).run()
        self.sort_run()
        self.set_output()
        self.end()

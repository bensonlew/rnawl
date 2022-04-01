# -*- coding: utf-8 -*-
# __author__ = 'guhaidong'

import os
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool

class MegahitAgent(Agent):
    """
    进行megahit拼接
    version: v1.0
    author: guhaidong
    last_modify: 2017.09.08
    """
    def __init__(self, parent):
        super(MegahitAgent, self).__init__(parent)
        options = [
            {"name": "fastq1", "type": "infile", "format": "sequence.fastq"},  # 输入文件,l
            {"name": "fastq2", "type": "infile", "format": "sequence.fastq"},  # 输入文件,r
            {"name": "fastqs", "type": "infile", "format": "sequence.fastq"},  # 输入文件,s 可不传
            {"name": "sample_name", "type": "string"},  #输入样品名
            {"name": "cpu", "type": "int", "default": 5},  # 拼接线程数，默认5
            {"name": "mem", "type": "int", "default": 10},  # 拼接使用内存，默认10
            {"name": "mem_mode", "type": "string", "default": "mem"},
            # 拼接内存模式,mem表示根据'mem'参数，minimum表示最小内存，moderate表示普通内存
            {"name": "min_contig", "type": "int", "default": 300},   # 最短contig值
            {"name": "mink", "type": "int", "default": 47},  # 最小kmer值
            {"name": "maxk", "type": "int", "default": 97},  # 最大kmer值
            {"name": "step", "type": "int", "default": 10},  # kmer步长
            {"name": "contig", "type": "outfile", "format": "sequence.fasta"},  # 输出文件,sample.contig.fa
        ]
        self.add_option(options)
        self.step.add_steps("Megahit")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)
        self._memory_increase_step = 150  # 每次重运行增加内存150G by guhaidong @ 20180427--20190530改为增加150Gqingchen.zhang

    def stepstart(self):
        self.step.Megahit.start()
        self.step.update()

    def stepfinish(self):
        self.step.Megahit.finish()
        self.step.update()

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option('fastq1'):
            raise OptionError('必须输入fastq1', code="31301101")
        if not self.option('fastq2'):
            raise OptionError('必须输入fastq2', code="31301102")
        if not self.option('sample_name'):
            raise OptionError('必须输入样本名', code="31301103")
        if self.option('mem_mode') not in ['mem', 'minimum', 'moderate']:
            raise OptionError('内存模式错误，选择[mem|minimum|moderate]之一', code="31301104")
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = self.option('cpu')
        # tmp_mem = int(self.option('mem')) + 50 * self._rerun_time  # 每次因拼接失败而重运行的内存增加50G by GHD @ 20180320
        # self._memory = '%sG' % tmp_mem
        # self.logger.info('megahit use memory : ' + self._memory)
        self._memory = "{}G".format(self.option('mem'))

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(MegahitAgent, self).end()


class MegahitTool(Tool):
    def __init__(self, config):
        super(MegahitTool, self).__init__(config)
        self.gcc = self.config.SOFTWARE_DIR + '/gcc/5.1.0/bin:' + self.config.SOFTWARE_DIR + "/program/Python/bin"
        self.gcc_lib = self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64:' + self.config.SOFTWARE_DIR + "/program/Python/lib"
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)
        self.megahit_path = '/bioinfo/metaGenomic/megahit/'

    def megahit_run(self):
        """
        进行megahit拼接
        :return:
        """
        if os.path.exists(self.output_dir + '/' + self.option('sample_name') + '.contig.fa'):
            return
        if os.path.exists(self.work_dir + '/run'):
            shutil.rmtree(self.work_dir + '/run')
        cmd = self.megahit_path + 'megahit -1 %s -2 %s '\
           % (self.option('fastq1').prop['path'], self.option('fastq2').prop['path'],)
        if self.option('fastqs').is_set:
            cmd += '-r %s ' % (self.option('fastqs').prop['path'])
        cmd += '-o %s --k-min %s --k-max %s --k-step %s --min-contig-len %s '\
            % (self.work_dir + '/run',
            self.option('mink'),
            self.option('maxk'),
            self.option('step'),
            self.option('min_contig'))
        if self.option('mem_mode') == 'mem':
            member_byte_format = self.option('mem') * 1000000000
            cmd += '-m %s' % (member_byte_format)
        elif self.option('mem_mode') == 'minimum':
            mem_mode_code = 0
            cmd += '--mem-flag %s' % (mem_mode_code)
        elif self.option('mem_mode') == 'moderate':
            mem_mode_code = 1
            cmd += '--mem-flag %s' % (mem_mode_code)
        self.logger.info("运行megahit拼接")
        command = self.add_command("megahit", cmd, ignore_error=True)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行megahit完成")
        elif command.return_code in [1, -9, -7, 250,137,245, 255, 247, 249]:  # 加入return_code检测，megahit在sanger超出内存的返回值为250
            # 加入255 by ghd @ 20180712.z,加入247 @20190221
            # 加入249 by qingchen.zhang @ 20190614 软件返回码为-7
            # 加入249 by qingchen.zhang @ 20190709 返回码为-7，软件返回码为-6
            # 加入245 by qingchen.zhang @ 20191204 软件返回码为-11
            self.logger.info("return code: %s" % command.return_code)
            self.add_state('memory_limit', 'memory is low!')   # add memory limit error by guhaidong @ 20180320
        elif command.return_code is None:
            self.logger.info("return code: %s" % command.return_code)
            self.add_state('memory_limit', 'memory is low!')
        else:
            self.logger.info("return code: %s" % command.return_code)
            self.set_error("megahit运行出错!", code="31301101")

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        out_fa = self.output_dir + '/' + self.option('sample_name') + '.contig.fa'
        if os.path.exists(out_fa):
            os.remove(out_fa)
        os.link(self.work_dir + '/run/final.contigs.fa', out_fa)
        self.option('contig').set_path(out_fa)
        self.logger.info("设置megahit分析结果目录成功")

    def run(self):
        """
        运行
        :return:
        """
        super(MegahitTool, self).run()
        self.megahit_run()
        self.set_output()
        self.end()
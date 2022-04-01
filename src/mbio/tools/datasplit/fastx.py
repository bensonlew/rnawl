# -*- coding: utf-8 -*-
# __author__ = 'xuting'

"""fastxtoolkit  用于统计碱基质量信息"""
import os
import errno
from biocluster.config import Config
from biocluster.tool import Tool
from biocluster.agent import Agent
from biocluster.core.exceptions import OptionError


class FastxAgent(Agent):
    """
    fastxtoolkit
    """
    def __init__(self, parent=None):
        super(FastxAgent, self).__init__(parent)
        self._run_mode = "ssh1"
        options = [
            {'name': 'sample_info', 'type': "infile", 'format': 'datasplit.miseq_split'},  # 样本拆分信息表
            {'name': 'unzip_path', 'type': "string"}  # bcl2fastq的下机输出目录
        ]
        self.add_option(options)

    def check_option(self):
        """
        参数检测
        """
        if not self.option('sample_info').is_set:
            raise OptionError("参数sample_info不能为空")
        if not self.option('unzip_path'):
            raise OptionError("参数unzip_path不能为空")
        return True

    def set_resource(self):
        """
        设置所需要的资源
        """
        self._cpu = 3
        self._memory = ''


class FastxTool(Tool):
    """
    """
    def __init__(self, config):
        super(FastxTool, self).__init__(config)
        self._version = 1.0
        self.fastx_dir = "bioinfo/seq/fastx_toolkit_0.0.14/"
        self.python_dir = "program/Python/bin/python"
        self.q20q30_stat = os.path.join(Config().SOFTWARE_DIR, "datasplit/bin/q20q30_stat.py")
        self.gnuplot = os.path.join(Config().SOFTWARE_DIR, "bioinfo/plot/gnuplot/bin")
        self.lib_path = os.path.join(Config().SOFTWARE_DIR, "library/lib")
        self.option('sample_info').get_info()
        file_name = os.listdir(self.option("unzip_path"))
        self.fastqs = list()
        for n in file_name:
            self.fastqs.append(os.path.join(self.option("unzip_path"), n))
        self.fastx = list()

    def make_ess_dir(self):
        """
        为软件fastxtoolkit的运行创建必要的运行目录
        """
        fastx_dir = os.path.join(self.work_dir, "fastx")
        dir_list = [fastx_dir]
        for name in dir_list:
            try:
                os.makedirs(name)
            except OSError as exc:
                if exc.errno == errno.EEXIST and os.path.isdir(name):
                    pass
                else:
                    raise OSError("创建目录失败")

    def fastxtoolkit(self):
        """
        统计碱基质量
        """
        cmd_list = list()
        i = 0
        log_path = os.path.join(self.work_dir, "log.sh")
        with open(log_path, 'wb') as w:
            my_str = "fastxtoolkit"
            w.write(my_str.center(79, "#"))
            w.write("\n")
        for fastq in self.fastqs:
            i += 1
            file_name = os.path.join(self.work_dir, "fastx", os.path.basename(fastq) + ".fastxstat")
            self.fastx.append(file_name)
            cmd = self.fastx_dir + "fastx_quality_stats" + " -i " + fastq + " -o " + file_name
            self.logger.debug(cmd)
            with open(log_path, 'wb') as w:
                w.write(cmd)
                w.write(cmd)
            command = self.add_command("fastx_quality_stats" + str(i), cmd)
            cmd_list.append(command)
        for mycmd in cmd_list:
            self.logger.info("开始运行fastx_quality_stats")
            mycmd.run()
        self.wait()
        for mycmd in cmd_list:
            if mycmd.return_code == 0:
                self.logger.info("运行fastx_quality_stats完成")
            else:
                self.set_error("运行fastx_quality_stats出错")

    def fastx_nucl_dist(self):
        """
        根据碱基质量统计文件绘制质量分布图和box图
        """
        cmd_list = list()
        i = 0
        self.set_environ(PATH=self.gnuplot, LD_LIBRARY_PATH=self.lib_path, LIBRARY_PATH=self.lib_path)
        for fastxstat in self.fastx:
            i += 1
            nucl_name = os.path.join(self.work_dir, "fastx", os.path.basename(fastxstat) + ".nucl.png")
            box_name = os.path.join(self.work_dir, "fastx", os.path.basename(fastxstat) + ".box.png")
            cmd = (self.fastx_dir + "script/fastx_nucleotide_distribution_graph.sh -i " + fastxstat
                   + " -o " + nucl_name)
            command = self.add_command("fastx_nucleotide_distribution_graph.sh" + str(i), cmd)
            cmd_list.append(command)
            cmd = (self.fastx_dir + "script/fastq_quality_boxplot_graph.sh -i " + fastxstat
                   + " -o " + box_name)
            command = self.add_command("fastq_quality_boxplot_graph.sh" + str(i), cmd)
            cmd_list.append(command)
        for mycmd in cmd_list:
            self.logger.info("开始运行" + mycmd.name)
            mycmd.run()
        self.wait()
        for mycmd in cmd_list:
            if mycmd.return_code == 0:
                self.logger.info(mycmd.name + "完成")
            else:
                self.set_error(mycmd.name + "发生错误")

    def q20_q30(self):
        cmd_list = list()
        i = 0
        for fastq in self.fastqs:
            i += 1
            file_name = os.path.join(self.work_dir, "fastx", os.path.basename(fastq) + ".q20q30")
            cmd = (self.python_dir + " " + self.q20q30_stat + " -i " + fastq + " -o " + file_name)
            command = self.add_command("q20q30_stat" + str(i), cmd)
            cmd_list.append(command)
        for mycmd in cmd_list:
            self.logger.info("开始运行" + mycmd.name)
            mycmd.run()
        self.wait()
        for mycmd in cmd_list:
            if mycmd.return_code == 0:
                self.logger.info(mycmd.name + " 统计完成")
            else:
                self.set_error(mycmd.name + " 统计出错")

    def run(self):
        super(FastxTool, self).run()
        self.make_ess_dir()
        self.fastxtoolkit()
        self.fastx_nucl_dist()
        self.q20_q30()
        self.end()

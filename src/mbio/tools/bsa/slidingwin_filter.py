# -*- coding: utf-8 -*-
# __author__ = 'wangzhaoyue'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class SlidingwinFilterAgent(Agent):
    """
    bsa分析阈值过滤
    version v1.0.1
    author: wangzhaoyue
    last_modify: 2018.02.24
    """
    def __init__(self, parent):
        super(SlidingwinFilterAgent, self).__init__(parent)
        options = [
            {"name": "slidingwin_file", "type": "infile", "format": "bsa.vcf"},  # sliding-win.result文件
            {"name": "slid_file", "type": "infile", "format": "bsa.vcf"},  # sliding-win.slid.result文件,用于计算阈值对应关系
            {"name": "region_type", "type": "string", "default": "quantile"},  # 阈值类型，分位值quantile/index值index
            {"name": "region_value", "type": "string"}  # 阈值
        ]
        self.add_option(options)
        self.step.add_steps("slidingwin_filter")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.slidingwin_filter.start()
        self.step.update()

    def stepfinish(self):
        self.step.slidingwin_filter.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option('slidingwin_file'):
            raise OptionError('必须输入sliding-win.result', code="31501101")
        if not self.option('slid_file'):
            raise OptionError('必须输入sliding-win.slid.result', code="31501102")
        if self.option('region_type') not in ["quantile", "index"]:
            raise OptionError('阈值类型只有quantile和index两种', code="31501103")
        if self.option('region_value'):
            if float(self.option('region_value')) < 0 or float(self.option('region_value')) > 1:
                raise OptionError('阈值范围在0-1之间', code="31501104")
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 1
        self._memory = "5G"

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(SlidingwinFilterAgent, self).end()


class SlidingwinFilterTool(Tool):
    def __init__(self, config):
        super(SlidingwinFilterTool, self).__init__(config)
        self._version = "v1.0.1"
        self.R_path = 'program/R-3.3.1/bin/Rscript '
        self.slidingthread_path = self.config.PACKAGE_DIR + '/bsa/slidingthread.R'
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/library/lib64')
        self.region_value = 0
        self.column = 0
        self.quantile = 0
        self.index = 0

    def run(self):
        """
        运行
        :return:
        """
        super(SlidingwinFilterTool, self).run()
        self.logger.info(self.option('region_type'))
        self.logger.info(self.option('region_value'))
        self.filter_NA()
        self.run_slidingthread()
        self.get_quantile_index()
        self.select_slidingwin_result()
        self.set_output()
        self.end()

    def run_slidingthread(self):
        """
        输入滑窗结果，算出每个阈值对应的index值
        :return:
        """
        with open(self.option("slidingwin_file").prop['path'])as fr:
            firstline = fr.readline()
            tmp = firstline.strip().split("\t")
            if len(tmp) == 9:   # 当只有一个index值时取第4列
                self.column = 4
            else:   # 取delta_index值,取第6列
                self.column = 6
        cmd = self.R_path + self.slidingthread_path + ' --infile {} --outfile {} --col {}'.format(
            self.work_dir + '/sliding-win.slid.result.filter', self.work_dir + '/slidingwin.thesold', self.column)
        self.logger.info('slidingthread.R,计算阈值和index值的对应关系')
        command = self.add_command("slidingthread_cmd", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("slidingthread运行完成")
        else:
            self.set_error("slidingthread运行出错!", code="31501101")

    def get_quantile_index(self):
        """
        得到当前阈值及index值
        :return:
        """
        if self.option('region_type') == 'quantile':
            self.region_value = float(self.option('region_value')) if self.option('region_value') else 0.999
        else:
            self.region_value = float(self.option('region_value')) if self.option('region_value') else 0.8  # 根据参数确定阈值
        with open(self.work_dir + '/slidingwin.thesold.quantile')as fr, open(self.work_dir + '/quantile.index', 'w+')as fw:
            lines = fr.readlines()
            if self.option('region_type') == 'quantile':
                self.quantile = float(str(self.region_value)[:6])   # 对输入的分位数截取，保留四位小数
                for line in lines[1:]:
                    tmp = line.strip().split('\t')
                    if float(tmp[0]) == self.quantile:   # 根据分位数，找到对应的index值，截取，保留四位小数
                        self.index = float(tmp[1])
            else:
                self.index = float(self.region_value)  # 对输入的index截取，保留四位小数
                quantile_list = []
                for line in lines[1:]:
                    tmp = line.strip().split('\t')
                    if float(tmp[1][:6]) == float(str(self.index)[:6]):  # 根据index，找到对应的分位数值，放在quantile_list中
                        quantile_list.append(float(tmp[0]))
                if len(quantile_list) != 0:
                    quantile_list.sort()
                    self.quantile = quantile_list[-1]
                else:
                    self.quantile = 0
            fw.write("quantile\tindex\n")
            fw.write(str(self.quantile) + '\t' + str(self.index) + '\n')

    def select_slidingwin_result(self):
        """
        根据类型及对应的阈值，筛选表格中对应列大于等于此阈值的行存入新文件,
        """
        with open(self.option('slidingwin_file').prop['path'])as fr, open(self.work_dir + '/sliding-win.threshold.select', 'w+')as fw:
            select_line = []
            lines = fr.readlines()
            fw.write(lines[0])
            for line in lines[1:]:
                tmp = line.strip().split("\t")
                if float(tmp[-6]) >= self.index:
                    select_line.append(line)
                    fw.write(line)
            if len(select_line) == 0:
                self.set_error("筛选结果为空，请重新确定阈值再运行！", code="31501102")
                self.set_error("筛选结果为空，请重新确定阈值再运行！", code="31501106")

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        for root, dirs, files in os.walk(self.output_dir):
            for names in files:
                os.remove(os.path.join(root, names))
        self.logger.info("设置结果目录")
        try:
            os.link(self.work_dir + "/sliding-win.threshold.select", self.output_dir + "/sliding-win.threshold.select")
            os.link(self.work_dir + "/quantile.index", self.output_dir + "/quantile.index")
            self.logger.info("设置slidingwin_filter分析结果目录成功")

        except Exception as e:
            self.logger.info("设置slidingwin_filter分析结果目录失败{}".format(e))
            self.set_error("设置slidingwin_filter分析结果目录失败%s",variables=(e), code="31501103")
        self.logger.info('设置文件夹路径成功')

    def filter_NA(self):
        """
        过滤掉sliding-win.slid.result文件中第四列为NA的
        :return:
        """
        with open(self.option("slid_file").prop['path'])as fr, open(self.work_dir + '/sliding-win.slid.result.filter', 'w+')as fw:
            lines = fr.readlines()
            for line in lines:
                tmp = line.strip().split('\t')
                if tmp[-2] != 'NA':
                    fw.write("\t".join(tmp) + "\n")

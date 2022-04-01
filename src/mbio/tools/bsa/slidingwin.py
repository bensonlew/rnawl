# -*- coding: utf-8 -*-
# __author__ = 'wangzhaoyue'

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os


class SlidingwinAgent(Agent):
    """
    bsa分析亲本子代深度过滤
    version v1.0.1
    author: wangzhaoyue
    last_modify: 2018.02.22
    """
    def __init__(self, parent):
        super(SlidingwinAgent, self).__init__(parent)
        options = [
            {"name": "index_file", "type": "infile", "format": "bsa.vcf"},  # index-calc.result.index文件
            {"name": "col", "type": "string", "default": "1,2,3"},  # the col of chr pos index
            {"name": "method", "type": "string", "default": "bp"},  # 模式，bp/num
            {"name": "win", "type": "string", "default": "2000000"},  # the window size
            {"name": "step", "type": "string", "default": "10000"},  # the step size
            {"name": "abs", "type": "string", "default": "0"}  # 没有亲本的时候该值为1，有亲本的时候该值为0
        ]
        self.add_option(options)
        self.step.add_steps("slidingwin")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.slidingwin.start()
        self.step.update()

    def stepfinish(self):
        self.step.slidingwin.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option('index_file'):
            raise OptionError('必须输入index-calc.result.index', code="31501001")
        if self.option('method') not in ["bp", "num"]:
            raise OptionError('模式只有bp和num两种', code="31501002")
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 1
        self._memory = "10G"

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(SlidingwinAgent, self).end()


class SlidingwinTool(Tool):
    def __init__(self, config):
        super(SlidingwinTool, self).__init__(config)
        self._version = "v1.0.1"
        self.R_path = 'program/R-3.3.1/bin/Rscript '
        self.slidingwin_path = self.config.PACKAGE_DIR + '/bsa/slidingwin.R'
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/library/lib64')

    def run(self):
        """
        运行
        :return:
        """
        super(SlidingwinTool, self).run()
        self.run_slidingwin()
        self.set_output()
        self.end()

    def run_slidingwin(self):
        """
        运行slidingwin.R,计算滑窗
        """
        cmd = self.R_path + self.slidingwin_path + ' --infile {} --outfile {} --col {}  --win {} --step {} --method {} ' \
                                                   '--abs {}'\
            .format(self.option('index_file').prop['path'], self.work_dir + '/sliding-win', self.option('col'),
                    self.option('win'), self.option('step'), self.option('method'), self.option("abs"))
        self.logger.info('运行slidingwin.R,计算滑窗')
        command = self.add_command("slidingwin_cmd", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("slidingwin运行完成")
        else:
            self.set_error("slidingwin运行出错!", code="31501001")

    def set_output(self):
        """
        将结果文件复制到output文件夹下面 --这里加上output中文件已经存在要删除的判断 modified by HONGDONG
        :return:
        """
        for root, dirs, files in os.walk(self.output_dir):
            for names in files:
                os.remove(os.path.join(root, names))
        self.logger.info("设置结果目录")
        with open(self.work_dir + "/sliding-win.result", "r") as f:
            lines = f.readlines()
            if len(lines) < 2:
                self.set_error("当前滑窗方式中突变位点个数太少，导致sliding-win.result文件为空，请重设滑窗方式和滑窗策略", code="31501002")
        try:
            os.link(self.work_dir + "/sliding-win.result", self.output_dir + "/sliding-win.result")
            os.link(self.work_dir + "/sliding-win.slid.result", self.output_dir + "/sliding-win.slid.result")
            self.logger.info("设置slidingwin分析结果目录成功")

        except Exception as e:
            self.logger.info("设置slidingwin分析结果目录失败{}".format(e))
            self.set_error("设置slidingwin分析结果目录失败%s",variables=(e), code="31501003")

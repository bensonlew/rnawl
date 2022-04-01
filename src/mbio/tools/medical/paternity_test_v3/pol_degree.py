# -*- coding: utf-8 -*-
# __author__ = 'hongyu.chen'

from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.config import Config
import datetime


class PolDegreeAgent(Agent):
    """
    (亲子)多重、杂捕污染度计算
    用法：
        必须设置ref_dp_xls、sample_id参数
        如果调用tool时不重写dp的默认值，计算杂捕的污染度
        如果设置dp的值，计算多重的污染度
    """
    def __init__(self, parent=None):
        super(PolDegreeAgent, self).__init__(parent)
        options = [
            {'name': 'ref_dp_xls', 'type': 'infile', 'format': 'paternity_test_V2.xls'},
            {'name': 'dp', 'type': 'float', 'default': 50},
            {'name': 'sample_id', 'type': 'string'}
        ]
        self.add_option(options)
        self.step.add_steps('PolDegree')
        self.on('start', self.step_start)
        self.on('end', self.step_end)
        self.queue = "gada"

    def step_start(self):
        self.step.PolDegree.start()
        self.step.update()

    def step_end(self):
        self.step.PolDegree.finish()
        self.step.update()

    def check_options(self):
        if not self.option("ref_dp_xls").is_set:
            raise OptionError("必须设置ref.dp.xls文件")

        if not self.option('sample_id'):
            raise OptionError('必须给定sample_id')

        dp = self.option("dp")
        if dp < 0:
            raise OptionError("dp的值必须大于零")

    def set_resource(self):
        self._cpu = 1
        self._memory = "10G"

    def end(self):
        super(PolDegreeAgent, self).end()


class PolDegreeTool(Tool):
    def __init__(self, config):
        super(PolDegreeTool, self).__init__(config)
        self.dict = {'0.3': '0.7', '0.35': '0.65', '0.4': '0.6'}
        self.pol_degree = list()
        dp = self.option("dp")

        if dp <= 100:
            self.max = 0.95
            self.min = 0.05
        elif dp <= 300:
            self.max = 0.97
            self.min = 0.03
        else:
            self.max = 0.98
            self.min = 0.02

    def caculate_pol_degree(self):
        for key, value in self.dict.items():
            num = num_pol = 0
            with open(self.option("ref_dp_xls").prop["path"], 'r') as f:
                for line in f:
                    line_split = line.rstrip("\n").split(' ')
                    if (float(line_split[5]) > self.min and float(line_split[5])< float(key)) \
                        or (float(line_split[5]) > float(value) and float(line_split[5]) < self.max):
                        num += 1
                        num_pol += 1
                    else:
                        num += 1
                pol_percent = float(num_pol)/num if num != 0 else 0.0
                result = (key, pol_percent)
                self.pol_degree.append(result)

    def judge_result(self):
        checker = 0
        pollution = ''
        for item in self.pol_degree:
            if item[0] == '0.4' and item[1] > 0.1:
                checker += 1
            if item[0] == '0.35' and item[1] > 0.05:
                checker += 1
            if item[0] == '0.3' and item[1] > 0.03:
                checker += 1
        if checker == 0:
            pollution = '/'
        if checker == 1:
            pollution = '+'
        if checker == 2:
            pollution = '++'
        if checker == 3:
            pollution = '+++'
        sample_id = self.option('sample_id')

        update_data = {
            'sample_id': sample_id,
            'pollution': pollution,
            'date': datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        data_list_update = list()
        data_list_update.append(update_data)
        query_tuple = ('sample_id', 'sample_id')
        self.api.api('medical.paternity_test_v3.paternity_test_v3').add_to_mongo(data_list_update, 'sg_sample_qc', 'update', query_tuple)

    def run_pol_degree(self):
        self.caculate_pol_degree()
        self.judge_result()
        self.end()

    def run(self):
        super(PolDegreeTool, self).run()
        self.run_pol_degree()

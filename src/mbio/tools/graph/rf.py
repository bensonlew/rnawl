# coding=utf-8
"""
a random forest tool
version 0.0.0
author: yingnn
date: 2017.10.16

"""
import os
import math
from unipath import Path
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import rf_const as CONSTS


class RfAgent(Agent):
    '''random forest tool agent

    '''
    def __init__(self, parent):
        super(RfAgent, self).__init__(parent)
        options = [{"name": CONSTS.input,
                    "type": "infile",
                    "format": "toolapps.table_rf"}, ]
        self.add_option(options)

        self.step.add_steps(CONSTS.step_rf)

        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.step_rf.start()
        self.step.update()

    def step_end(self):
        self.step.step_rf.finish()
        self.step.update()

    def check_options(self):
        if not self.option(CONSTS.input):
            raise OptionError("input file error")

    def set_resource(self):
        self._cpu = CONSTS.n_cpu
        total = os.path.getsize(self.option(CONSTS.input).prop["path"])
        total = 20
        memory = "{}G".format(total)
        self._memory = memory

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_info = [[".", "", u"分析结果"],
                       ["./train_data_describe.xls", 'xls',
                        u'训练集信息，信息包括样本数，样本种类，特征（变量）数，'
                        u'样本（行），特征（列）的缺失值（NA）占比'],
                       ["./train_bad.xls", 'xls',
                        u'缺失值过多（> 50%）而弃用的训练集样本'],
                       ["./test_data_describe.xls", 'xls',
                        u'测试集（如果有测试集）信息，信息包括样本数，特征数，'
                        u'缺失值占比'],
                       ["./test_bad.xls", 'xls',
                        u'缺失值过多（> 50%）的测试集样本，仍用于类别的预测，'
                        u'缺失值被填充，分类结果可靠性差'],
                       ["./data_used.xls", 'xls',
                        u'模型用到的训练集样本和个数'],
                       ["./test_probability_predict_rf.xls", 'xls',
                        u'模型预测到的测试样本的类别和可能性'],
                       ["./feature_importances_rf.xls", 'xls',
                        u'特征的重要性'],
                       ["./rfe_scores.xls", 'xls',
                        u'使用重要的特征，'
                        u'以精度（accuracy，多于两个类别）或 AUC （两个类别）为指标，'
                        u'十字交叉验证模型准确率；'
                        u'横轴为使用到的排名靠前的特征的个数，纵轴为指标值'], ]
        result_dir.add_relpath_rules(result_info)
        super(RfAgent, self).end()


class RfTool(Tool):
    """random forest tool

    """
    def __init__(self, config):
        super(RfTool, self).__init__(config)

    def run_main(self):
        input = self.option(CONSTS.input).prop['new_table']
        cmd = os.path.join(self.config.SOFTWARE_DIR, CONSTS.sh)
        cmd = ' '.join([CONSTS.sh, input])
        cmd = self.add_command(CONSTS.step_rf, cmd)
        cmd.run()
        self.logger.info('%s start running' % cmd)
        self.wait(cmd)
        if cmd.return_code == 0:
            self.logger.info('%s finished' % cmd.name)
        else:
            self.set_error('something wrong with %s running' % cmd.name)
            raise Exception("error occurs")

    def set_output(self):
        files = ['train_data_describe.txt',
                 'data_used.txt',
                 'train_bad.txt',
                 'test_bad.txt',
                 'feature_importances_rf.txt',
                 'test_probability_predict_rf.txt',
                 'test_data_describe.txt',
                 'rfe_scores.txt',
                 'model.txt']
        for file in files:
            if os.path.exists(os.path.join(self.work_dir, file)):
                p = Path(os.path.join(self.work_dir, file))
                p1 = Path(self.output_dir)
                p1 = p1.child(p.stem + '.xls')
                p.hardlink(p1)
                # p.hardlink(os.path.join(self.output_dir, file))


    def run(self):
        super(RfTool, self).run()
        self.logger.info('main func run')
        self.run_main()
        self.logger.info('set output files')
        self.set_output()
        self.logger.info('end func run')
        self.end()


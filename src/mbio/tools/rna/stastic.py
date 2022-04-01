## !/mnt/ilustre/users/sanger-dev/app/program/Python/bin/python
# -*- coding: utf-8 -*-
# __author__ = "moli.zhou"

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from biocluster.config import Config
import os
import shutil

class StasticAgent(Agent):
    """
		进行个数的统计
    """
    def __init__(self, parent):
        super(StasticAgent, self).__init__(parent)
        options = [#输入的参数
            {"name": "data","type": "infile", "format": "ref_rna.protein_regulation.txt"}, #待统计数据
            {"name": "row", "type": "int"} #对数据中的第几列进行统计
        ]
        self.add_option(options)
        self.step.add_steps("Stastic")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.Stastic.start()
        self.step.update()

    def stepfinish(self):
        self.step.Stastic.finish()
        self.step.update()


    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option('data').is_set:
            raise OptionError("请确认输入")
        if not self.option('row'):
            raise OptionError('请输入要统计的是第几列')
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 10
        self._memory = '100G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"]
        ])
        result_dir.add_regexp_rules([
            ["stastic.txt", "txt", "统计信息"],
        ])
        super(StasticAgent, self).end()


class StasticTool(Tool):
    """
    蛋白质互作组预测tool
    """
    def __init__(self, config):
        super(StasticTool, self).__init__(config)
        self._version = '1.0.1'
        self.python_path = 'program/Python/bin/'
        self.script_path = Config().SOFTWARE_DIR + '/bioinfo/rna/scripts/'
        # self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/program/perl/perls/perl-5.24.0/lib/site_perl/5.24.0/Bio')

    # python phmmer_process.py 1e-180  PlantTFDB-all_TF_pep.fas test.fas planttfdb_family_vs_tfid.txt
    def run_stas(self):

        cmd = '{}python {}stastic_script.py {} {}'.\
            format(self.python_path,self.script_path,self.option('data').prop['path'], self.option('row'))

        self.logger.info(cmd)
        self.logger.info("开始运行统计")
        cmd = self.add_command("cmd", cmd).run()
        self.wait(cmd)
        if cmd.return_code == 0:
            self.logger.info("运行个数统计成功")
        else:
            self.logger.info("运行个数统计出错")

    def set_output(self):
        """
        将结果文件link到output文件夹下面
        :return:
        """
        for root, dirs, files in os.walk(self.output_dir):
            for names in files:
                os.remove(os.path.join(root, names))
        self.logger.info("设置结果目录")

        f = 'stastic_result.txt'
        os.link(self.work_dir + '/' + f, self.output_dir + '/' + f)
        self.logger.info('设置文件夹路径成功')

    def run(self):
        super(StasticTool, self).run()
        self.run_stas()
        self.set_output()
        self.end()

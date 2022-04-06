## !/mnt/ilustre/users/sanger-dev/app/miniconda2/bin/python
# -*- coding: utf-8 -*-
# __author__ = "moli.zhou"
#last_modify:20161108

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from biocluster.config import Config
import os
import sys
import shutil

class TfPredictAgent(Agent):
    """
    利用hmmer软件，进行转录因子预测
    version v1.0
    author: moli.zhou
    last_modify: 2016.11.4
    """
    def __init__(self, parent):
        super(TfPredictAgent, self).__init__(parent)
        options = [#输入的参数
            {"name": "query_amino", "type": "infile", "format": "sequence.fasta"},  # 上游输入的氨基酸文件（含与差异基因的对应）
            {"name": "diff_gene_id", "type": "string"},
            {"name": "database", "type": "string", "default": "AnimalTFDB"}, #还有PlantTFDB和AnimalTFDB
        ]
        self.add_option(options)
        self.step.add_steps("TfPredict")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.TfPredict.start()
        self.step.update()

    def stepfinish(self):
        self.step.TfPredict.finish()
        self.step.update()


    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        database_list = ["PlantTFDB", "AnimalTFDB", "iTAK"]
        # if not self.option('query_amino').is_set:
        #     raise OptionError("必须输入氨基酸序列")
        if not self.option('diff_gene_id'):
            raise OptionError("请输入差异基因对应id关系")
        if self.option('database') not in database_list:  # species的判定有问题
            raise OptionError("database选择不正确")
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
            ["TfPredict.txt", "txt", "转录因子预测信息"],
        ])
        super(TfPredictAgent, self).end()


class TfPredictTool(Tool):
    """
    蛋白质互作组预测tool
    """
    def __init__(self, config):
        super(TfPredictTool, self).__init__(config)
        self._version = '1.0.1'
        self.python_path = 'program/Python/bin/'
        self.perl_path = 'program/perl/perls/perl-5.24.0/bin/'
        self.script_path = Config().SOFTWARE_DIR + '/bioinfo/rna/scripts/'
        self.ref_path = Config().SOFTWARE_DIR + '/database/refGenome/TF/'
        self.itak_path = Config().SOFTWARE_DIR + '/bioinfo/rna/iTAK-1.6b/'
        self.set_environ(LD_LIBRARY_PATH=self.config.SOFTWARE_DIR + '/program/perl/perls/perl-5.24.0/lib/site_perl/5.24.0/Bio')

    # python phmmer_process.py 1e-180  PlantTFDB-all_TF_pep.fas test.fas planttfdb_family_vs_tfid.txt
    def run_tf(self):
        if self.option("database") == 'PlantTFDB':
            ref = self.ref_path + "plant/planttfdb.hmm"
            family = self.ref_path + "plant/family_vs_DBD_plant.txt"
            tf_cmd = "{}python {}TF_process_plant.py {} {} {} {}"\
                .format(self.python_path,self.script_path, ref, self.option("query_amino").prop['path'],family,
                        self.option("diff_gene_id"))

        elif self.option("database") == 'iTAK':
            tf_cmd = '{}python {}TF_process_iTAK.py {} {} {} {}'.\
                format(self.python_path, self.script_path, Config().SOFTWARE_DIR+'/'+self.perl_path,self.itak_path,
                       self.option("query_amino").prop['path'],self.option("diff_gene_id"))

        elif self.option('database') == 'AnimalTFDB':
            ref = self.ref_path + "animal/animaltfdb.hmm"
            family = self.ref_path + "animal/family_vs_DBD_animal_2.0.txt"
            tf_cmd = '{}python {}TF_process_animal.py {} {} {} {}'\
                .format(self.python_path,self.script_path,ref,self.option("query_amino").prop['path'],family,
                        self.option("diff_gene_id"))

        self.logger.info(tf_cmd)
        self.logger.info("开始运行TFPredict")
        cmd = self.add_command("tf_cmd", tf_cmd).run()
        self.wait(cmd)
        if cmd.return_code == 0:
            self.logger.info("运行TFPredict成功")
        else:
            self.logger.info("运行TFPredict出错")

    def set_output(self):
        """
        将结果文件link到output文件夹下面
        :return:
        """
        for root, dirs, files in os.walk(self.output_dir):
            for names in files:
                os.remove(os.path.join(root, names))
        self.logger.info("设置结果目录")

        if self.option("database") == 'PlantTFDB':
            f = 'TF_result.txt'
            if os.path.exists(f):   
                os.link(self.work_dir + '/' + f, self.output_dir + '/' + f)
                self.logger.info('设置文件夹路径成功')
            else:
                self.logger.info('没有找到对应的转录因子家族')
                
        elif self.option("database") == 'AnimalTFDB':
            f = 'TF_result.txt'
            if os.path.exists(f):   
                os.link(self.work_dir + '/' + f, self.output_dir + '/' + f)
                self.logger.info('设置文件夹路径成功')
            else:
                self.logger.info('没有找到对应的转录因子家族')
        elif self.option("database") == 'iTAK':
            f = 'TF_result.txt'
            if os.path.exists(f):
                os.link(self.work_dir + '/' + f, self.output_dir + '/' + f)
                self.logger.info('设置文件夹路径成功')
            else:
                self.logger.info('没有找到对应的转录因子家族')

    def run(self):
        super(TfPredictTool, self).run()
        self.run_tf()
        self.set_output()
        self.end()


# -*- coding: utf-8 -*-
# __author__ = 'zhangyitong'


import unittest
import os
import glob
import re
import pandas as pd
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError


class SelectTableAgent(Agent):
    def __init__(self, parent):
        super(SelectTableAgent, self).__init__(parent)
        options = [
            {"name": "origin_table", "type": "infile", "format": "sequence.profile_table"},
            {"name": "group", "type": "infile", "format": "meta.otu.group_table"},
            {"name": "group_method", "type": "int", "default": 0},  # 1：求和，2：求均值，3：求中位数
            {"name": "scale", "type": "bool", "default": True},  ## 是否标准化,代谢用
            {"name": "merge", "type": "string", 'default': 'nomerge'},
            {"name": "select_table", "type": "outfile", "format": "sequence.profile_table"},
            {"name": "select_origin_abu", "type": "outfile", "format": "sequence.profile_table"},  ## 有scale时使用
            {"name": "map_table", "type": "outfile", "format": "ref_rna_v2.common"},  # file for id mapping

        ]
        self.add_option(options)
        self.step.add_steps("table_select")
        self.on('start', self.step_start)
        self.on('end', self.step_finish)

    def step_start(self):
        self.step.table_select.start()
        self.step.update()

    def step_finish(self):
        self.step.table_select.finish()
        self.step.update()

    def check_options(self):
        if not self.option("origin_table").is_set:
            raise OptionError("必须设置输入数据表")
        if not self.option("group").is_set:
            raise OptionError("必须设置输入样本分组文件")
        return True

    def set_resource(self):
        self._cpu = 2
        self._memory = "20G"

    def end(self):
        super(SelectTableAgent, self).end()


class SelectTableTool(Tool):
    def __init__(self, config):
        super(SelectTableTool, self).__init__(config)
        self._version = "v1.0"
        self.set_environ(PATH=os.path.join(self.config.SOFTWARE_DIR, 'program/Python/bin'))
        self.program = {
            'python': 'miniconda2/bin/python',
        }
        self.script = {
            'profile': os.path.join(self.config.PACKAGE_DIR, 'tool_lab/circular_heatmap/profile_select.py'),
        }
        self.file = {
            'processed': os.path.join(self.work_dir, 'processed_table.txt'),
            'mapping': os.path.join(self.work_dir, 'mapped_id.txt'),
            'outfile': os.path.join(self.output_dir, 'select_table.xls')
        }
        self.n = 1

    def run(self):
        super(SelectTableTool, self).run()
        self.id_processing()
        self.run_group()
        self.set_output()
        self.end()

    def convert(self, map_dict, raw):
        new_id = re.sub(r'\W', "", raw)
        if new_id not in map_dict.keys():
            map_dict[new_id] = raw
        else:
            new_id = new_id + '_' + str(self.n)
            map_dict[new_id] = raw
            self.n += 1
        return new_id

    def id_processing(self):
        map_dict = dict()
        raw_df = pd.read_table(self.option('origin_table').prop['path'], header=0)
        col1 = raw_df.columns[0]
        raw_df[col1] = raw_df[col1].apply(lambda x: self.convert(map_dict, x))
        raw_df.to_csv(self.file['processed'], sep='\t', index=False, header=True)
        with open(self.file['mapping'], 'w') as m:
            m.write('metab_id\tMetabolite\n')
            for k, v in map_dict.items():
                v = v.replace(';', '_').replace(',', '_')
                v = v.replace('.', '_@').replace('\'', '’').replace('(', '（').replace(')', '）').replace(':', '：')
                m.write(k + '\t' + v + '\n')

    def run_group(self):
        self.logger.info("start table_select")
        cmd = '{} {} -i {} -o {} '.format(self.program['python'], self.script['profile'], self.file['processed'],
                                          self.file['outfile'])
        cmd += "-g {} -gm {} ".format(self.option('group').prop['path'], str(self.option("group_method")))
        cmd += "-merge {}".format(self.option("merge"))
        if self.option('scale'):
            cmd += ' --scale -odir {}'.format(self.output_dir)
        self.logger.info(cmd)
        command = self.add_command('table_select', cmd, ignore_error=True).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("table_select succeed")
        elif command.return_code in [-9, 1]:  # modified return_code by guhaidong @ 20180628
            self.add_state("memory_limit", "memory is low!")
        else:
            self.set_error("table_select failed", code="34002501")

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        if os.path.exists(self.file['outfile']):
            self.logger.info("设置select_table成功")
            self.option("select_table", self.file['outfile'])
        if os.path.exists(self.file['mapping']):
            self.logger.info("设置map_table成功")
            self.option("map_table", self.file['mapping'])
        if self.option('scale'):
            before_scale_table = self.output_dir + "/select_before_scale.xls"
            if os.path.exists(before_scale_table):
                self.option("select_origin_abu").set_path(before_scale_table)


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run script to do test.
    """

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "target_depth_bwa_{}_{}".format(random.randint(1000, 9999), random.randint(1000, 9999)),
            "type": "tool",
            "name": "tool_lab.target_depth.bwa",
            "instant": False,
            "options": dict(
                sample_name='HY_1',
                fastq_l='/mnt/ilustre/users/sanger-dev/sg-users/zhangyitong/308_TagGeneInRawSeq/script/bwa/HY_1.1.fq',
                fastq_r='/mnt/ilustre/users/sanger-dev/sg-users/zhangyitong/308_TagGeneInRawSeq/script/bwa/HY_1.2.fq',
                seq_file='/mnt/ilustre/users/sanger-dev/sg-users/zhangyitong/308_TagGeneInRawSeq/script/bwa/NDM_5.fasta',
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTests([TestFunction("test")])
    unittest.TextTestRunner(verbosity=2).run(suite)
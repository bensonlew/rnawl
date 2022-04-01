# -*- coding: utf-8 -*-
# __author__ = 'zoujiaxun'

import ConfigParser
import glob
import os
import re
import shutil
import subprocess
import unittest
import configparser
from biocluster.agent import Agent
from biocluster.config import Config
from biocluster.core.exceptions import OptionError
from biocluster.tool import Tool
from skbio.parse.sequences import parse_fasta
from Bio import SeqIO


class KnownPirnaNewAgent(Agent):
    """
    已知miRNA鉴定
    """

    def __init__(self, parent):
        super(KnownPirnaNewAgent, self).__init__(parent)
        options = [
            {'name': 'species', 'type': 'string'},
            {'name': 'fa', 'type': 'infile', 'format': 'small_rna.fasta'},
            {'name': 'length', 'type': 'string'},
            {'name': 'config', 'type': 'infile', 'format': 'small_rna.common'},

            {'name': 'pirna_table', 'type': 'outfile', 'format': 'small_rna.common'},
            {'name': 'pirna_stat', 'type': 'outfile', 'format': 'small_rna.common'},
        ]
        self.add_option(options)
        self.step.add_steps("known_pirna")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.known_pirna.start()
        self.step.update()

    def stepfinish(self):
        self.step.known_pirna.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        # if self.option("species").lower() == "null" and not self.option("orgnism_list").is_set:
        #     raise OptionError("必须指定具体物种")
        # if self.option("orgnism_list").is_set and self.option("whole"):
        #     self.logger.info("当物种列表和全物种同时提供的时候，以全物种进行分析")
        # if not self.option("clean_fa").is_set:
        #     raise OptionError("必须提供质控后的FASTA文件")
        # if not self.option("config").is_set:
        #     raise OptionError("必须提供配置文件")
        # with open(self.option("clean_fa").prop["path"], "r") as f:
        #     line = f.readline()
        #     if not re.compile(r'^>\S+_x\d+').match(line):
        #         raise OptionError("质控后的序列文件格式不对")

        self.set_resource()
        return True

    def set_resource(self):
        """
        设置所需资源，需在类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 1
        self._memory = '50G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "已知miRNA鉴定结果目录"]
        ])
        result_dir.add_regexp_rules([
            [r"pdfs", "", "已知miRNA二级结构文件目录"],
            ["known_mirna_count.xls", "xls", "已知miRNA定量count表"],
            ["known_mirna_norm.xls", "xls", "已知miRNA定量norm表"],
            ["filtered.fa", "", "过滤FASTA文件"],
        ])
        super(KnownPirnaNewAgent, self).end()


class KnownPirnaNewTool(Tool):
    def __init__(self, config):
        super(KnownPirnaNewTool, self).__init__(config)
        self.cmd_path = "bioinfo/align/bowtie-1.2.3-linux-x86_64/bowtie"
        self.database = "/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/smallrna/piRNA/"
        self.file = {
            'pirna_table': os.path.join(self.output_dir, 'known_pirna_{}.bwt'.format(self.option('length'))),
            'pirna_stat': os.path.join(self.output_dir, 'known_pirna_{}.txt'.format(self.option('length')))
        }

    def run(self):
        """
        运行
        :return:
        """
        super(KnownPirnaNewTool, self).run()
        self.bowtie()
        self.stat()
        self.set_output()
        self.end()

    def bowtie(self):
        self.logger.info("开始进行已知piRNA鉴定")
        database_path = os.path.join(self.database, self.option('species'), self.option('species')+'_len', '{}_len_{}.fa'.format(self.option('species'), self.option('length')))
        cmd = '{} -f -v 1 --nofw {} {} {}'.format(self.cmd_path, database_path, self.option('fa').path, self.file['pirna_table'])
        cmd_name = 'bowtie'
        command = self.add_command(cmd_name, cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("{} 运行成功".format(cmd_name))
        elif command.return_code is None:
            self.logger.warn("运行{}出错，返回值为None，尝试重新运行".format(cmd_name))
            command.rerun()
            self.wait()
            if command.return_code is 0:
                self.logger.info("{} 运行成功".format(cmd_name))
            else:
                self.set_error("运行%s>>>%s出错", variables=(cmd_name, cmd), code="33704306")
        else:
            self.set_error("运行%s>>>%s出错", variables=(cmd_name, cmd), code="33704307")

    def stat(self):
        conf = configparser.ConfigParser()
        conf.read(self.option('config').path)
        sample = dict()
        for i in conf['NAME']:
            sample[i.upper()] = conf.get('NAME', i)
        pirna_table = self.file['pirna_table']

        with open(pirna_table, 'r') as pir, open(self.file['pirna_stat'], 'w') as ps:
            ps.write('pirna_id' + '\t' + 'count' + '\t' + 'sample' + '\n')
            for line in pir.readlines():
                sample_id = line.strip().split('\t')[0]
                sample2sample = sample[re.findall(r"^(.*?)_", sample_id)[0]]
                pir_id = line.strip().split('\t')[2]
                count = re.findall(r"_x(.*?)$", sample_id)[0]
                ps.write(pir_id + '\t' + count + '\t' + sample2sample + '\n')

    def set_output(self):
        """
        将结果文件link到output文件夹下面
        :return:
        """
        self.option('pirna_table').set_path(self.file['pirna_table'])
        self.option('pirna_stat').set_path(self.file['pirna_stat'])


class TestFunction(unittest.TestCase):

    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "KnownPirna_" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "small_rna.srna.known_pirna",
            "instant": False,
            "options": dict(
                species="ssc",
                clean_fa="/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/smallrna/test/uniq_1000000.fa",
                config="/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/smallrna/test/qc_file.config",
                # #mirdeep2_version='0.1.3',
                # whole=False,
                # database="pmiren"
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()

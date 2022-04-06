# -*- coding: utf-8 -*-
# __author__ = 'hao.gao'

from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool
import os

class RrnaAlignAgent(Agent):
    '''
    16s比对数据库（GTDB或自定义数据库）
    '''
    def __init__(self, parent):
        super(RrnaAlignAgent, self).__init__(parent)
        options = [
            {"name": "type", "type": "string", "default":"GTDB"},  # GTDB or custom
            {"name": "sample", "type": "string"},
            {'name': 'input', 'type': 'infile', "format": "sequence.fasta"},  ## 比对基因组的16s序列
            {'name': 'custom_fa', 'type': 'infile', "format": "sequence.fasta"},  ##custom时的数据库序列
            {'name': 'blast', 'type': 'outfile', "format": "sequence.profile_table"},  ## 16s比对的结果
        ]
        self.add_option(options)

    def check_options(self):
        """
        参数二次检查
        :return:
        """
        if self.option('type') not in ["GTDB", "custom"]:
            raise OptionError("必须输入正确的type：GTDB or custom !")
        else:
            if self.option('type') not in ["custom"] and self.option('custom_fa').is_set:
                raise OptionError("必须输入正确的type：{}错误的!".format(self.option('type')))
        if not self.option('input').is_set:
            raise OptionError("必须输入input!")
        return True

    def set_resource(self):
        """
        设置资源
        :return:
        """
        self._cpu = '2'
        self._memory = '10G'

    def end(self):
        """
        运行结束
        :return:
        """
        super(RrnaAlignAgent, self).end()

class RrnaAlignTool(Tool):
    """
    16s比对数据库（GTDB或自定义数据库）
    """
    def __init__(self, config):
        super(RrnaAlignTool, self).__init__(config)
        self.python = "/miniconda2/bin/python"
        self.blast_path = "/bioinfo/align/ncbi-blast-2.3.0+/bin/"
        self.gtdb_database = self.config.SOFTWARE_DIR + "/database/GTDB/GTDB_rrna/bac120_ssu_reps_r95"


    def run(self):
        """
        运行
        :return:
        """
        super(RrnaAlignTool, self).run()
        if self.option("type") in ['GTDB']:
            self.run_blast()
            self.set_output()
        elif self.option("type") in ['custom']:
            self.run_makedb()
            self.run_blast()
            self.set_output()
        self.end()

    def run_makedb(self):
        """
        自定义时构建数据库索引
        :return:
        """
        cmd = '{}makeblastdb -dbtype nucl -in {} -out {}'.format(self.blast_path, self.option("custom_fa").prop['path'],self.work_dir+"/ref")
        self.logger.info(cmd)
        command = self.add_command("run_makedb", cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info('run_makedb运行成功！')
        else:
            self.set_error('run_makedb运行失败！')

    def run_blast(self):
        """
        GTDB or custom时的blast的比对
        :return:
        """
        sample =self.option("sample")
        database = ''
        if self.option("type") in ['GTDB']:
            database =self.gtdb_database
        elif self.option("type") in ['custom']:
            database = self.work_dir+"/ref"
        cmd = '{}blastn -query {} -db {} -out {} -evalue 1e-5 -outfmt 6 -num_alignments 30 -num_threads 2'.format(self.blast_path, self.option("input").prop['path'],
                                                                 database, self.output_dir + "/"+ sample +".16s.blast.xls")
        self.logger.info(cmd)
        command = self.add_command("run_blast", cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info('run_blast运行成功！')
        else:
            self.set_error('run_blast运行失败！')

    def set_output(self):
        sample =self.option("sample")
        if os.path.getsize(self.output_dir + "/"+ sample +".16s.blast.xls") >0:
            self.option("blast", self.output_dir + "/"+ sample +".16s.blast.xls")
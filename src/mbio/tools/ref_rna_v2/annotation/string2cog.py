# -*- coding: utf-8 -*-
# __author__ = 'wangbixuan'
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
import xml.etree.ElementTree as ET
import subprocess


class String2cogAgent(Agent):
    """
    author: wangbixuan
    last_modified: 20160719
    """

    def __init__(self, parent):
        super(String2cogAgent, self).__init__(parent)
        options = [
            {"name": "blastout", "type": "infile", "format": "ref_rna_v2.blast_xml"},
            {"name": "cog_list", "type": "outfile", "format": "ref_rna_v2.cog_list"},
            {"name": "cog_table", "type": "outfile", "format": "ref_rna_v2.cog_table"},
            # {"name": "E_value", "type": "float", "default": "1E-6"},
            # {"name": "min_Coverage", "type": "string", "default": "50%"}
        ]
        self.add_option(options)
        self.step.add_steps('string2cog')
        self.on('start', self.step_start)
        self.on('end', self.step_end)

    def step_start(self):
        self.step.string2cog.start()
        self.step.update()

    def step_end(self):
        self.step.string2cog.finish()
        self.step.update()

    def check_options(self):
        if self.option("blastout").is_set:
            document = ET.parse(self.option("blastout").prop['path'])
            root = document.getroot()
            db = root.find('BlastOutput_db')
            if db.text.endswith('string'):
                pass
            else:
                raise OptionError("BLAST比对数据库不支持", code = "33702801")
        else:
            raise OptionError("必须提供BLAST结果文件", code = "33702802")

    def set_resource(self):
        self._cpu = 20
        self._memory = '5G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ['cog_list.xls', 'xls', 'COG注释表'],
            ['cog_summary.xls', 'xls', 'COG注释总结'],
            ['cog_table.xls', 'xls', 'COG注释详细表']
        ])
        super(String2cogAgent, self).end()


class String2cogTool(Tool):

    def __init__(self, config):
        super(String2cogTool, self).__init__(config)
        self._version = '1.0'  # to be changed
        # self.cmd_path = '{}/miniconda2/bin/python {}/bioinfo/annotation/scripts/string2cog.py'.format(
        #    self.config.SOFTWARE_DIR, self.config.SOFTWARE_DIR)
        self.cmd_path = '{}/miniconda2/bin/python {}/bioinfo/annotation/scripts/string2cog_mongo.py'.format(
            self.config.SOFTWARE_DIR, self.config.SOFTWARE_DIR)  # 修改（sqlite3数据库换成mongo数据库；2016.11.15；zengjing）

    def run(self):
        super(String2cogTool, self).run()
        self.run_string2cog()

    def run_string2cog(self):
        cmd = self.cmd_path
        cmd += ' %s %s' % (self.option("blastout").prop['path'], self.work_dir)
        self.logger.info('运行string2cog.py')
        self.logger.info(cmd)
        try:
            subprocess.check_output(cmd, shell=True)
            self.logger.info('运行string2cog.py完成')
            outfiles = ['cog_list.xls', 'cog_summary.xls', 'cog_table.xls']
            for item in outfiles:
                linkfile = self.output_dir + '/' + item
                if os.path.exists(linkfile):
                    os.remove(linkfile)
                os.link(self.work_dir + '/' + item, linkfile)
            self.option('cog_list', self.output_dir + '/cog_list.xls')
            self.option('cog_table', self.output_dir + '/cog_table.xls')
            self.end()
        except subprocess.CalledProcessError:
            self.set_error('运行string2cog.py出错', code = "33702803")

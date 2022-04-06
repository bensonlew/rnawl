# -*- coding: utf-8 -*-
# __author__ = 'zengjing'
# modified 2017.05.04
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
from biocluster.core.exceptions import OptionError
import xml.etree.ElementTree as ET
import subprocess


class String2cogv9Agent(Agent):
    """
    输入比对到string库的table文件或xml文件进行cog注释及统计
    """

    def __init__(self, parent):
        super(String2cogv9Agent, self).__init__(parent)
        options = [
            {"name": "blastout", "type": "infile", "format": "prok_rna.blast_xml"},
            {"name": "string_table", "type": "infile", "format": "prok_rna.blast_table"},
            {"name": "cog_list", "type": "outfile", "format": "prok_rna.cog_list"},
            {"name": "cog_table", "type": "outfile", "format": "prok_rna.cog_table"},
            {"name": "cog_summary", "type": "outfile", "format": "prok_rna.cog_summary"}
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
        if not self.option("blastout").is_set and not self.option("string_table").is_set:
            raise OptionError("必须提供比对到string库的xm文件或table文件")
        # if self.option("blastout").is_set:
        #     document = ET.parse(self.option("blastout").prop['path'])
        #     root = document.getroot()
        #     db = root.find('BlastOutput_db')
        #     if db.text.endswith('string'):
        #         pass
        #     else:
        #         raise OptionError("BLAST比对数据库不支持")

    def set_resource(self):
        self._cpu = 10
        self._memory = '20G'

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
        super(String2cogv9Agent, self).end()


class String2cogv9Tool(Tool):

    def __init__(self, config):
        super(String2cogv9Tool, self).__init__(config)
        self._version = '1.0'  # to be changed
        self.python = self.config.SOFTWARE_DIR + '/miniconda2/bin/python'
        # self.cog_xml = self.config.SOFTWARE_DIR + '/bioinfo/annotation/scripts/string2cog_v9.py'
        self.cog_xml = self.config.PACKAGE_DIR + "/denovo_rna_v2/String2Cog.pl"
        self.cog_table = self.config.PACKAGE_DIR + "/denovo_rna_v2/cog_annot.py"
        self.perl = '/program/perl-5.24.0/bin/perl'
        self.db = self.config.SOFTWARE_DIR + "/database/COG/cog.db"

    def run(self):
        super(String2cogv9Tool, self).run()
        self.run_string2cog()

    def run_string2cog_old(self):
        if self.option("blastout").is_set:
            cmd = '{} {} {} {}'.format(self.python, self.cog_xml, self.option('blastout').prop['path'], self.work_dir)
        else:
            cmd = '{} {} {} {}'.format(self.python, self.cog_table, self.option('string_table').prop['path'], self.work_dir)
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
            self.option('cog_summary', self.output_dir + '/cog_summary.xls')
            self.option('cog_table', self.output_dir + '/cog_table.xls')
            self.end()
        except subprocess.CalledProcessError:
            self.set_error('运行string2cog.py出错')

    def run_string2cog(self):
        if self.option("blastout").is_set:
            cmd = '{} {} -i {} --format blastxml -db {} -e 1e-3 -o {}/tmp_out'.format(self.perl, self.cog_xml, self.option('blastout').prop['path'], self.db, self.work_dir)
        else:
            cmd = '{} {} -i {} --format blasttable -db {} -e 1e-3 -o {}/tmp_out'.format(self.perl, self.cog_xml, self.option('blastout').prop['path'], self.db, self.work_dir)
        self.logger.info('运行string2cog.pl')
        self.logger.info(cmd)
        cmd_obj = self.add_command("string_cog", cmd).run()
        self.wait(cmd_obj)
        if cmd_obj.return_code == 0:
            self.logger.info('运行string2cog.pl完成')
            for file in os.listdir(self.work_dir + "/tmp_out"):
                file_path = os.path.join(self.work_dir + "/tmp_out", file)
                if file == "cog_list.xls":
                    os.link(file_path, self.output_dir + "/cog_list.xls")
                if file == "cog_summary.xls":
                    os.link(file_path, self.output_dir + "/cog_summary.xls")
                if file == "cog_table.xls":
                    os.link(file_path, self.output_dir + "/cog_table.xls")
            self.option('cog_list', self.output_dir + '/cog_list.xls')
            self.option('cog_summary', self.output_dir + '/cog_summary.xls')
            self.option('cog_table', self.output_dir + '/cog_table.xls')
            self.end()
        else:
            self.set_error("string2cog出错")
            raise Exception("string2cog出错")

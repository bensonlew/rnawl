# -*- coding: utf-8 -*-
# __author__ = 'zouxuan'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re


class SignalpAnnoAgent(Agent):
    """
    分泌蛋白预测 v1.0
    author: zouxuan
    last_modify: 20180203
    """

    def __init__(self, parent):
        super(SignalpAnnoAgent, self).__init__(parent)
        options = [
            {"name": "query", "type": "infile", "format": "sequence.fasta"},
            {"name": "type", "type": "string", "default": "gram-"},  # 菌株类型gram-，gram+，euk
            {"name": "out_format", "type": "string", "default": "summary"},
            {"name": "result", "type": "outfile", "format": "sequence.profile_table"},
            {"name": "sample", "type": "string"},  # 样品名称
            {"name": "d_score", "type": "bool", "default": False}  # 是否添加D值
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("query").is_set:
            raise OptionError("必须设置输入序列文件", code="31102001")
        if not self.option("sample"):
            raise OptionError("必须设置样品名称", code="31102002")
        if self.option("type") not in ["gram-", "gram+", "euk"]:
            raise OptionError("菌株选择类型错误，只能为gram-，gram+，euk+", code="31102003")
        return True

    def set_resource(self):
        self._cpu = 1
        self._memory = '10G'
        ##guanqing.zou 20180905
        if self.option("type") == 'euk':
            self._cpu = 2
            self._memory = '8G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        super(SignalpAnnoAgent, self).end()


class SignalpAnnoTool(Tool):
    def __init__(self, config):
        super(SignalpAnnoTool, self).__init__(config)
        self._version = "1.0"
        self.signalp_path = self.config.SOFTWARE_DIR + '/bioinfo//Genomic/Sofware/signalp-4.1/signalp'
        self.signalp_sh = 'bioinfo/Genomic/Sofware/signalp-4.1/signalp.sh'
        self.python_script_path = self.config.PACKAGE_DIR + '/bacgenome/add_gene_info.py '
        self.python_path = '/program/Python/bin/python'

    def run(self):
        """
        运行
        :return:
        """
        super(SignalpAnnoTool, self).run()
        self.signalp_predict()
        self.end()

    def signalp_predict(self):
        linkfile = self.work_dir + '/' + os.path.basename(self.option('query').prop['path'])
        if os.path.exists(linkfile):
            os.remove(linkfile)
        os.link(self.option('query').prop['path'], linkfile)
        cmd = '{} {} {} {} {} {} {}'.format(self.signalp_sh, self.signalp_path,
                                         self.option('type'),
                                         self.option('out_format'), self.work_dir, linkfile,
                                         self.work_dir + '/' + self.option('type') + '_SignalP.txt')
        command = self.add_command('predict', cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("signalp预测成功")
        else:

            self.set_error("signalp预测失败", code="31102001")
            self.set_error("signalp预测失败", code="31102002")

        signalp_outdir = self.output_dir + '/' + self.option("sample") + '_' + self.option(
            'type').capitalize() + '_SignalP.txt'
        if self.option('type') == 'euk':  #zouguanqing 20180611  真菌输出文件名称不包含type
            signalp_outdir = self.output_dir + '/' + self.option("sample") + '_SignalP.xls'

        self.get_result(self.work_dir + '/' + self.option('type') + '_SignalP.txt', signalp_outdir)

        with open(signalp_outdir) as fr:
            lines = fr.readlines()
            lines_num = len(lines)

        if lines_num < 2:
            os.remove(signalp_outdir)
        else:
            cmd2 = '{} {} -i {} -f {} -s {} '.format(self.python_path, self.python_script_path, signalp_outdir,
                                                     self.option('query').prop['path'],
                                                     self.option('sample'))
            command2 = self.add_command('get_file', cmd2).run()
            self.wait(command2)
            if command2.return_code == 0:
                self.logger.info("文件生成成功")
            else:
                self.set_error("文件生成失败", code="31102003")
                self.set_error("文件生成失败", code="31102004")

    def get_result(self, in_file, out_file):
        regex = re.compile('\s+')
        with open(in_file, 'rb') as f, open(out_file, 'w') as o:
            if self.option("d_score"):
                o.write("Gene ID\tCmax\tC pos\tYmax\tY pos\tSmax\tS pos\tSmean\tNetworks-used\tD-score\n")
            else:
                o.write("Gene ID\tCmax\tC pos\tYmax\tY pos\tSmax\tS pos\tSmean\tNetworks-used\n")
            lines = f.readlines()
            for line in lines[2:]:
                line = regex.split(line.strip())
                if line[9] == 'Y':
                    if self.option("d_score"):
                        w_line = '\t'.join(
                            [line[0], line[1], line[2], line[3], line[4], line[5], line[6], line[7], line[11],line[8]])
                    else:
                        w_line = '\t'.join(
                            [line[0], line[1], line[2], line[3], line[4], line[5], line[6], line[7], line[11]])
                    o.write(w_line + '\n')

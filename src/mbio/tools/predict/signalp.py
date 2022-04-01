# -*- coding: utf-8 -*-
# __author__ = 'zouxuan'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re
from mbio.packages.metagenomic.common import link_file,link_dir
import glob


class SignalpAgent(Agent):
    """
    分泌蛋白预测 v1.0
    author: guhaidong
    last_modify: 20181114
    """

    def __init__(self, parent):
        super(SignalpAgent, self).__init__(parent)
        options = [
            {"name": "query", "type": "infile", "format": "sequence.fasta"},
            {"name": "query_dir", "type": "infile", "format": "sequence.fasta_dir"},
            {"name": "type", "type": "string", "default": "gram-"},  # 菌株类型gram-，gram+，euk
            {"name": "out_format", "type": "string", "default": "summary"},
            {"name": "result", "type": "outfile", "format": "sequence.profile_table"},
            {"name": "sample", "type": "string", "default": "signalp"},  # 样品名称
            {"name": "mem", "type": "int", "default": 3}  # 预测使用的内存
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option("query").is_set and not self.option("query_dir").is_set:
            raise OptionError("must type into query option", code="33301101") # 必须设置输入序列文件

        if self.option("type") not in ["gram-", "gram+", "euk"]:
            raise OptionError("isolation type must be gram-,gram+ or euk", code="33301102") # 菌株选择类型错误，只能为gram-，gram+，euk

        return True

    def set_resource(self):
        if self.option("query").is_set:
            self._cpu = 2
        else:
            self._cpu = 10
        self._memory = str(self.option('mem')) + "G"

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        super(SignalpAgent, self).end()


class SignalpTool(Tool):
    def __init__(self, config):
        super(SignalpTool, self).__init__(config)
        self._version = "1.0"
        self.signalp_path = self.config.SOFTWARE_DIR + '/bioinfo//Genomic/Sofware/signalp-4.1/signalp'
        self.signalp_sh = 'bioinfo//Genomic/Sofware/signalp-4.1/signalp.sh'
        # self.python_script_path = self.config.PACKAGE_DIR + '/bacgenome/add_gene_info.py '
        self.python_path = '/program/Python/bin/python'

    def run(self):
        """
        运行
        :return:
        """
        super(SignalpTool, self).run()
        self.signalp_predict()
        self.end()

    def signalp_predict(self):
        if self.option('query').is_set:
            linkfile = self.work_dir + '/' + os.path.basename(self.option('query').prop['path'])
            link_file(self.option('query').prop['path'], linkfile)
            cmd = '{} {} {} {} {} {} {}'.format(self.signalp_sh, self.signalp_path,
                                         self.option('type'),
                                         self.option('out_format'), self.work_dir, linkfile,
                                         self.work_dir + '/' + self.option('type') + '_SignalP.txt')
            command = self.add_command('predict', cmd).run()
            self.wait(command)
            if command.return_code == 0:
                self.logger.info("signalp预测成功")
            else:
                self.set_error("signalp预测失败", code="33301101")
                self.set_error("signalp预测失败", code="33301102")
        else:
            linkfile = self.work_dir + '/input_fa'
            link_dir(self.option('query_dir').prop['path'], linkfile)
            count = 0
            command_list = []
            for file in os.listdir(linkfile):
                count += 1
                file_path = os.path.join(linkfile, file)
                cmd = '{} {} {} {} {} {} {}'.format(self.signalp_sh, self.signalp_path,
                                         self.option('type'),
                                         self.option('out_format'), self.work_dir, file_path,
                                         self.work_dir + '/' + str(count) + self.option('type') + '_SignalP.txt')
                command = self.add_command('predict' + str(count), cmd).run()
                command_list.append(command)
                if len(command_list) == 10 or count == len(os.listdir(linkfile)):
                    self.wait(*command_list)
                    for one in command_list:
                        if one.return_code == 0:
                            self.logger.info("signalp %s success" % one.name)
                        else:
                            self.set_error("signalp %s error" , variables=( one.name), code="33301103")
                    command_list = []
        signalp_outdir =  self.output_dir + '/' + self.option("sample") + '_' + self.option('type').capitalize() + '_SignalP.txt'
        self.get_result(signalp_outdir)
        self.option("result", signalp_outdir)
        # self.get_result(self.work_dir + '/' + self.option('type') + '_SignalP.txt',signalp_outdir)
        # cmd2 = '{} {} -i {} -f {} -s {} '.format(self.python_path, self.python_script_path,signalp_outdir,
        #                                          self.option('query').prop['path'],
        #                                          self.option('sample'))
        # command2 = self.add_command('get_file', cmd2).run()
    # self.wait(command2)
        # if command2.return_code == 0:
        #     self.logger.info("文件生成成功")
        # else:
        #     self.set_error("文件生成失败", code="31102003")
        #     self.set_error("文件生成失败", code="31102004")

    def get_result(self, out_file):
        with open(out_file, "w") as o:
            o.write("Gene ID\tCmax\tC pos\tYmax\tY pos\tSmax\tS pos\tSmean\tNetworks-used\n")
        regex = re.compile('\s+')
        for base_name in glob.glob(self.work_dir + '/*SignalP.txt'):
            in_file = os.path.join(self.work_dir, base_name)
            with open(in_file, 'rb') as f, open(out_file, 'a') as o:
                lines = f.readlines()
                for line in lines[2:]:
                    line = regex.split(line.strip())
                    query = line[0].rsplit("_1", 1)[0]  # 蛋白id轉為基因id
                    if line[9] == 'Y':
                        w_line = '\t'.join(
                            [query, line[1], line[2], line[3], line[4], line[5], line[6], line[7], line[11]])
                            # [line[0], line[1], line[2], line[3], line[4], line[5], line[6], line[7], line[11]])
                        o.write(w_line + '\n')

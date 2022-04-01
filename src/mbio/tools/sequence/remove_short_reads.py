# -*- coding: utf-8 -*-
# __author__ = 'zhouxuan'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os

class RemoveShortReadsAgent(Agent):
    """
    剔除fastq中较短的序列
    author: zhouxuan
    last_modify: 2017.06016
    """

    def __init__(self, parent):
        super(RemoveShortReadsAgent, self).__init__(parent)
        options = [
            {"name": "fastq_dir", "type": "infile", "format": "sequence.fastq_dir"},  # fastq文件夹
            {"name": "reasult_dir", "type": "outfile", 'format': "sequence.fastq_dir"}  # 输出文件夹
            ]
        self.add_option(options)

    def check_options(self):
        if not self.option("fastq_dir").is_set:
            raise OptionError("必须正确设置输入文件", code="34002401")
        return True

    def set_resource(self):
        self._cpu = 10
        self._memory = '15G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
            ])
        super(RemoveShortReadsAgent, self).end()


class RemoveShortReadsTool(Tool):
    def __init__(self, config):
        super(RemoveShortReadsTool, self).__init__(config)
        self._version = "1.0"
        self.set_environ(PATH=self.config.SOFTWARE_DIR + '/program/perl/perls/perl-5.24.0/bin')
        self.set_environ(PERLBREW_ROOT=self.config.SOFTWARE_DIR + '/program/perl')
        #self.perl_script = self.config.SOFTWARE_DIR + '/bioinfo/seq/scripts/remove_short_reads.pl'
        #self.perl_script_2 = self.config.SOFTWARE_DIR + "/bioinfo/seq/scripts/remove_short_reads.pair.pl"
        self.perl_script = self.config.PACKAGE_DIR + '/sequence/scripts/remove_short_reads.pl'
        self.software = 'program/parafly-r2013-01-21/bin/bin/ParaFly'
        self.perl_script_2 = self.config.PACKAGE_DIR + "/sequence/scripts/remove_short_reads.pair.pl"

    def run(self):
        """
        运行
        :return:
        """
        super(RemoveShortReadsTool, self).run()
        self.run_remove()
        self.set_output()
        self.end()

    def run_remove(self):
        list_path = os.path.join(self.option('fastq_dir').prop['path'], 'list.txt')
        result_list = os.path.join(self.output_dir, 'list.txt')
        sample_name = []
        l_dict = {}
        cmd_list = []
        with open(list_path, 'r') as l, open(result_list, 'w') as w:
            for line in l:
                line = line.strip("\n").split('\t')
                if len(line) == 2:
                    file_name = line[1] + ".sickle.s.fq"
                    cmd = '{} {} 50 {}'.format(self.perl_script,
                                               os.path.join(self.option('fastq_dir').prop['path'], line[0]),
                                               os.path.join(self.output_dir, file_name))
                    cmd_list.append(cmd)
                    w.write(line[1] + '.sickle.s.fq\t' + line[1] + '\n')
                else:
                    if line[1] not in sample_name:
                        sample_name.append(line[1])
                    l_dict[line[1]+'_'+line[2]] = os.path.join(self.option('fastq_dir').prop['path'], line[0])
                    if line[2] == 'l':
                        w.write(line[1] + '.sickle.1.fq\t' + line[1] + '\tl\n')
                    elif line[2] == 'r':
                        w.write(line[1] + '.sickle.2.fq\t' + line[1] + '\tr\n')
                    elif line[2] == 's':  # metagenome'PSE format add by zhujuan 20170925
                        file_name = line[1] + ".sickle.s.fq"
                        w.write(line[1] + '.sickle.s.fq\t' + line[1] + '\ts\n')
                        cmd = '{} {} 50 {}'.format(self.perl_script,
                                                   os.path.join(self.option('fastq_dir').prop['path'], line[0]),
                                                   os.path.join(self.output_dir, file_name))
                        cmd_list.append(cmd)
        for i in sample_name:
            cmd = '{} {} {} 50 {}'.format(self.perl_script_2, l_dict[i + "_l"], l_dict[i + '_r'],
                                          os.path.join(self.output_dir, i + ".sickle"))
            cmd_list.append(cmd)
        n = len(cmd_list)/15
        if len(cmd_list)%15 != 0:
            n += 1
        for i in range(0, n):
            cmd_file = os.path.join(self.work_dir, 'list_{}.txt'.format(i+1))
            wrong_cmd = os.path.join(self.work_dir, 'failed_cmd_{}.txt'.format(i+1))
            with open(cmd_file, 'w') as c:
                for n in range(0, 15):
                    if len(cmd_list) == 0:
                        break
                    cmd = cmd_list.pop()
                    c.write(cmd + "\n")
            p_cmd = '{} -c {} -CPU 10 -failed_cmds {}'.format(self.software, cmd_file, wrong_cmd)
            command = self.add_command('all_cmd_{}'.format(i+1), p_cmd).run()
            self.wait(command)
            if command.return_code == 0:
                self.logger.info("运行{}完成".format(command.name))
            else:
                self.set_error("运行%s出错", variables=(command.name), code="34002401")
                raise Exception("运行{}出错".format(command.name))

    def set_output(self):
        self.logger.info('开始设置输出结果文件')
        try:
            self.option('reasult_dir', self.output_dir)
            self.logger.info("设置输出结果文件正常")
        except Exception as e:
            raise Exception("设置输出结果文件异常——{}".format(e))

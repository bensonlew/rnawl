# -*- coding: utf-8 -*-
# __author__ = 'wangzhaoyue'

import os
import re
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool
from mbio.packages.ref_rna.trans_step import step_count

class MetageneAgent(Agent):
    """
    使用metagene软件进行基因预测
    version: Metagene
    author: wangzhaoyue
    last_modify: 2017.06.19
    """
    def __init__(self, parent):
        super(MetageneAgent, self).__init__(parent)
        options = [
            {"name": "cut_more_scaftig", "type": "infile", "format": "sequence.fasta"},
            # 输出文件，去掉小于最短contig长度的序列
            {"name": "sample_name", "type": "string"},  # 样本的名称
            {"name": "min_gene", "type": "string", "default": "100"},  # 输入最短基因长度，如100
            {"name": "fna", "type": "outfile", "format": "sequence.fasta"},  # 输出文件，样本的核酸序列
            {"name": "cut_more_fna", "type": "outfile", "format": "sequence.fasta"},  # 输出文件，样本去除最小值后的核酸序列
            #{"name": "faa", "type": "outfile", "format": "sequence.fasta"},  # 输出文件，样本的蛋白序列
        ]
        self.add_option(options)
        self.step.add_steps("Metagene")
        self._memory_increase_step = 20  # 每次重运行增加内存20G by guhaidong @ 20180428
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.Metagene.start()
        self.step.update()

    def stepfinish(self):
        self.step.Metagene.finish()
        self.step.update()

    def check_options(self):
        """
        检查参数
        :return:
        """
        if not self.option('cut_more_scaftig'):
            raise OptionError('必须输入样本去掉小于最短contig长度的序列文件', code="32201101")
        if not self.option('sample_name'):
            raise OptionError('必须输入样本的名称', code="32201102")
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        last modified: 20180320
        :return:
        """
        self._cpu = 1
        self._memory = "3G"  # 3G改到5G by GHD @ 20180428
        # tmp_mem = 5 * (self._rerun_time + 1)  #  用1G进行测试
        # if self._rerun_time == 3:
        #     tmp_mem = 50  # 预测设定最高值为50G
        # self._memory = '%sG' % (tmp_mem)

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(MetageneAgent, self).end()


class MetageneTool(Tool):
    def __init__(self, config):
        super(MetageneTool, self).__init__(config)
        self._version = "metagene"
        self.sh_path = 'program/sh ' + self.config.PACKAGE_DIR + '/gene_structure/scripts/'
        self.metagene_path = self.config.SOFTWARE_DIR + '/bioinfo/metaGenomic/metagene/'
        self.perl_path = '/miniconda2/bin/perl '
        # self.metagene_seqs_path = self.config.SOFTWARE_DIR + '/bioinfo/metaGenomic/scripts/metagene_seqs.pl '
        self.metagene_seqs_path = self.config.PACKAGE_DIR + '/gene_structure/scripts/metagene_seqs.pl '
        # self.cut_more_path = self.config.SOFTWARE_DIR + '/bioinfo/metaGenomic/scripts/cut_more.pl '
        self.cut_more_path = self.config.PACKAGE_DIR + '/sequence/scripts/cut_more.pl '
        self.emboss_path = "bioinfo/seq/EMBOSS-6.6.0/emboss/"
        self.rerun_num = 0

    def run(self):
        """
        运行
        :return:
        """
        super(MetageneTool, self).run()
        self.run_metagene()
        self.end()

    def treat_error(self, rerun_cmd, message):  # 尝试重复运行tool，目前此方法不需要 @ 20180320
        self.logger.info("开始执行错误cmd")
        if self.rerun_num < 2:
            # rerun_cmd.run()
            self.add_state('memory_limit', 'just wanna error!')
            self.rerun_num += 1
            self.logger.info("开始进行第%s次重复运行" % self.rerun_num)
            #self.wait(rerun_cmd)
            #if rerun_cmd.return_code == 0:
            #    self.logger.info("运行%s完成" % message)
            #else:
            #    self.treat_error(rerun_cmd, message)
        else:
            self.set_error("运行%s出错!", variables=(message), code="32201101")

    def run_metagene(self):
        """
        metagene [input] -m > [output]
        :return:
        """
        cmd = self.sh_path + 'metagene.sh %s %s %s' % (self.metagene_path, self.option('cut_more_scaftig').prop['path'],
                                                       self.work_dir + '/' + self.option('sample_name') + '.metagene.csv')
        command = self.add_command("metagene", cmd, ignore_error=True, shell=False)  # shell=False by ghd @ 20190301
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行metagene的cmd完成")
            self.run_metageneseqs()
        # elif command.return_code == 137:  # there is other return_code proudcted through memeory_limit by xieshichang 20200506
        elif command.return_code in [137, 134]:  # by xieshichang 20200506
            self.add_state('memory_limit', 'memory is low!')   # add memory limit error by guhaidong @ 20180320
        else:
            self.set_error("运行metagene的cmd出错!", code="32201102")

    def run_metageneseqs(self):
        """
        MetageneSeqs  -m [csv] -f [scaftig] -o [fna]
        :return:
        """
        cmd = self.perl_path + self.metagene_seqs_path + '-m %s -f %s -o %s' % (self.work_dir + '/' + self.option('sample_name') + '.metagene.csv', self.option('cut_more_scaftig').prop['path'], self.work_dir + '/' + self.option('sample_name') + '.metagene.fna')
        command = self.add_command("metageneseqs运行", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行MetageneSeqs的cmd完成")
            self.run_cut_more()
        else:
            self.set_error("运行MetageneSeqs的cmd运行出错!", code="32201103")

    def run_cut_more(self):
        """
        perl cut_more.pl [run_MetageneSeqs的输出文件] [最短contig长度] [输出文件的名称前缀]
        :return:
        """
        cmd = self.perl_path + self.cut_more_path + self.work_dir + '/' + self.option('sample_name') + '.metagene.fna' + ' ' + self.option('min_gene') + ' ' + self.option('sample_name') + '.metagene.fna'
        command = self.add_command("cut_more", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行cut_more完成")
            self.set_output()
        else:
            self.set_error("运行cut_more运行出错!", 32201104, code="32201104")
    '''
    def run_transeq(self):
        """
        transeq -sequence [fna.more] -table 11 -trim -outseq [faa]
        :return:
        """
        cmd = self.emboss_path + 'transeq -sequence %s -table 11 -trim -outseq %s' % (
            self.work_dir + '/' + self.option('sample_name') + '.metagene.fna.more' + self.option('min_gene'),
            self.work_dir + '/' + self.option('sample_name') + '.metagene.faa')
        command = self.add_command("transeq运行", cmd)
        command.run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("运行transeq的cmd完成")
        else:
            self.set_error("运行transeq的cmd运行出错!")

    def step_count(self):
        """
        调用函数，对挑选出来的序列,整理到文件夹汇总，并进行步长统计
        :return:
        """
        self.logger.info("开始预测部分进行步长统计")
        step_list = [200, 700]
        all_files = os.listdir(self.work_dir)
        for files in all_files:
            m = re.search(r'(\S+)\.metagene\.fna\.more', files)
            if m:
                fna_file = self.work_dir + '/' + files
                for step in step_list:
                    output1 = self.work_dir + '/' + m.group(1) + '_fna_step_' + str(step) + '.txt'
                    output2 = self.output_dir + '/' + m.group(1) + '_fna_step_' + str(step) + '.final.txt'
                    step_count(fna_file, output1, 20, step, output2)
            else:
                pass
        self.logger.info("预测部分步长统计结束")
    '''
    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        outfasta1 = self.output_dir + '/' + self.option('sample_name') + '.metagene.fna'
        outfasta2 = self.output_dir + '/' + self.option('sample_name') + '.metagene.more' + self.option('min_gene') + '.fa'
        if os.path.exists(outfasta1):
            os.remove(outfasta1)
        if os.path.exists(outfasta2):
            os.remove(outfasta2)
        #os.link(self.work_dir + '/' + self.option('sample_name') + '.metagene.fna',
        #        outfasta1)
        os.link(self.work_dir + '/' + self.option('sample_name') + '.metagene.fna.more' + self.option('min_gene'),
                outfasta2)
        #os.link(self.work_dir + '/' + self.option('sample_name') + '.metagene.faa',
        #        self.output_dir + '/' + self.option('sample_name') + '.metagene.faa')
        #self.option('fna').set_path(self.output_dir + '/' + self.option('sample_name') + '.metagene.fna')
        self.option('cut_more_fna').set_path(
             outfasta2)
        #    self.output_dir + '/' + self.option('sample_name') + '.metagene.more' + self.option('min_gene') + '.fa')
        #self.option('faa').set_path(self.output_dir + '/' + self.option('sample_name') + '.metagene.faa')
        self.logger.info("设置Metagene分析结果目录成功")

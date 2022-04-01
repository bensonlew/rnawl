# -*- coding: utf-8 -*-
# __author__ = 'guhaidong'
import os
from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool
import re


class IdbaAgent(Agent):
    """
    运用IDBA软件自带的fq2fa将序列进行合并
    version v1.0.1
    author: guhaidong
    last_modify: 2017.09.12
    """

    def __init__(self, parent):
        super(IdbaAgent, self).__init__(parent)
        options = [
            {"name": "fastq1", "type": "infile", "format": "sequence.fastq"},  # 输入文件,l
            {"name": "fastq2", "type": "infile", "format": "sequence.fastq"},  # 输入文件,r
            {"name": "fastqs", "type": "infile", "format": "sequence.fastq"},  # 输入文件,s 可不传
            {"name": "sample_name", "type": "string"},  # 拼接样品名称
            {"name": "min_contig", "type": "int", "default": 300},  # 最短contig值
            {"name": "mem", "type": "int", "default": 250},  # 拼接使用内存，不填写则根据样品自动判断
            {"name": "split_num", "type": "int", "default": 1},  # 拼接前拆分reads份数，1表示不拆分
            {"name": "mink", "type": "int", "default": 47},  # 最小kmer值
            {"name": "maxk", "type": "int", "default": 97},  # 最大kmer值
            {"name": "step", "type": "int", "default": 10},  # kmer步长
            {"name": "contig", "type": "outfile", "format": "sequence.fasta"},  # 输出文件,sample.contig.fa
        ]
        self.add_option(options)
        self.step.add_steps("IdbaFq2fa")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)
        self._memory_increase_step = 50  # 每次重运行增加内存50G by guhaidong @ 20180427

    def stepstart(self):
        self.step.IdbaFq2fa.start()
        self.step.update()

    def stepfinish(self):
        self.step.IdbaFq2fa.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option('fastq1'):
            raise OptionError('必须输入*l.fastq文件', code="31301001")
        if not self.option('fastq2'):
            raise OptionError('必须输入*r.fastq文件', code="31301002")
        if not self.option('contig'):
            raise OptionError('必须传入输出文件名', code="31301003")
        if self.option('split_num') < 1:
            raise OptionError('reads拆分数目不能小于1', code="31301004")
        if self.option('mem') > 250:
            raise OptionError('内存超出250G', code="31301005")
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 16
        if self.option('mem') < 2:
            self._memory = "2G"
        else:
            # tmp_mem = int(self.option('mem')) + 50 * self._rerun_time  # 每次因拼接失败而重运行的内存增加50G by GHD @ 20180320
            # self._memory = '%sG' % tmp_mem
            self._memory = '%sG' % (self.option('mem'))  # 改回 by guhaidong @ 20180427
        if self.option('split_num') > 1:
            self._memory = '250G'
        self.logger.info('idba use memory : ' + self._memory)

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["", "", ""]
        ])
        super(IdbaAgent, self).end()


class IdbaTool(Tool):
    def __init__(self, config):
        super(IdbaTool, self).__init__(config)
        self.gcc = self.config.SOFTWARE_DIR + '/gcc/5.1.0/bin'
        self.gcc_lib = self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64'
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)
        self.perl_path = '/program/perl/perls/perl-5.24.0/bin/perl '
        # self.split_fasta_path = self.config.SOFTWARE_DIR + '/bioinfo/metaGenomic/scripts/split_fa.pl '
        self.split_fasta_path = self.config.PACKAGE_DIR + '/sequence/scripts/split_fa.pl '
        # self.combine_contig_path = self.config.SOFTWARE_DIR + '/bioinfo/metaGenomic/scripts/combine_contig.pl '
        self.combine_contig_path = self.config.PACKAGE_DIR + '/metagenomic/scripts/combine_contig.pl '
        self._version = "v1.1.1"
        self.idba_path = '/bioinfo/metaGenomic/MaxBin-2.2.5/auxiliary/idba-1.1.1/bin/'
        # self.sh_path = '/bioinfo/metaGenomic/scripts/'
        # self.sh_path = "../../../../../.." + self.config.PACKAGE_DIR + '/sequence/scripts/'
        self.sh_path = self.config.PACKAGE_DIR + "/sequence/scripts/"
        if not self.option('sample_name'):
            self.sample_name = os.path.basename(self.option('fastq1').prop['path']).split('.sickle.l.fastq')[0]
        else:
            self.sample_name = self.option('sample_name')

    def run(self):
        """
        运行
        :return:
        """
        super(IdbaTool, self).run()
        self.run_fq2fa_merge()
        if self.option('fastqs').is_set:
            self.run_fq2fa_paired()
            self.run_fq2fa_cat()
        else:
            # sample_name = os.path.basename(self.option('fastq1').prop['path']).split('.sickle.l.fastq')[0]
            os.link(self.work_dir + '/' + self.sample_name + '.pair.tmp.fa',
                    self.work_dir + '/' + self.sample_name + '.fa')
        if self.option('split_num') > 1:
            self.run_split()
            self.run_idba(self.option('split_num'))
            self.run_combine_contig()
        else:
            self.run_idba()
        self.set_output()
        self.end()

    def run_fq2fa_merge(self):
        """
        运行IDBA的Fq2fa，将序列1,2进行合并
        """
        # sample_name = os.path.basename(self.option('fastq1').prop['path']).split('.sickle.l.fastq')[0]
        cmd1 = self.idba_path + 'fq2fa --merge %s %s %s' % (self.option('fastq1').prop['path'],
                                                            self.option('fastq2').prop['path'],
                                                            self.work_dir + '/' + self.sample_name + '.pair.tmp.fa')
        self.logger.info('运行Idba的Fq2fa，将序列1,2进行合并')
        command = self.add_command("fq2fa_cmd1", cmd1).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("1,2合并完成")
        else:
            self.set_error("1,2合并出错!", code="31301001")

    def run_fq2fa_paired(self):
        """
        如有S端，运行IDBA的Fq2fa，将序列s进行转换
        :return:
        """
        # sample_name = os.path.basename(self.option('fastq1').prop['path']).split('.sickle.l.fastq')[0]
        cmd2 = self.idba_path + 'fq2fa --paired %s %s' % (self.option('fastqs').prop['path'],
                                                          self.work_dir + '/' + self.sample_name + '.single.tmp.fa')
        self.logger.info('运行Idba的Fq2fa，将序列s进行转换')
        command = self.add_command("fq2fa_cmd2", cmd2).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("s转换完成")
        else:
            self.set_error("s转换出错!", code="31301002")

    def run_fq2fa_cat(self):
        """
        如有S端，需将S端转换的结果与1,2合并结果再合并
        :return:
        """
        self.logger.info('将所有序列合并')
        # sample_name = os.path.basename(self.option('fastq1').prop['path']).split('.sickle.l.fastq')[0]
        result_fa = self.work_dir + '/' + self.sample_name + '.fa'
        if os.path.exists(result_fa):
            os.remove(result_fa)
        fa_list = [self.work_dir + '/' + self.sample_name + '.pair.tmp.fa',
                   self.work_dir + '/' + self.sample_name + '.single.tmp.fa']
        cmd3 = self.sh_path + 'cat_seq.sh'
        for files in fa_list:
            # cmd3 = self.sh_path + 'cat_seq.sh %s %s ' % (files, result_fa)
            cmd3 += ' ' + files
        cmd3 += ' ' + result_fa
        self.logger.info("start cat seqs to {}".format(result_fa))
        command = self.add_command("fq2fa_cmd3", cmd3, shell=False).run()  # shell=True by ghd @ 20190301
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("run_fq2fa_cat done")
        else:
            self.set_error("run_fq2fa_cat error", code="31301003")

    def run_split(self):
        """
        对pe reads 和 se reads分别进行拆分
        :param num:拆分数量
        :return:
        """
        cmd = self.perl_path + self.split_fasta_path + ' n ' + str(self.option('split_num')) + ' ' + \
              self.work_dir + '/' + self.sample_name + '.fa ' + self.sample_name + ' ' + self.work_dir + '/split'
        self.logger.info("拆分fasta")
        command = self.add_command("split", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("拆分fasta完成")
        else:
            self.set_error("拆分fasta失败！", code="31301004")

    def run_idba(self, number=0):  # 增加对样品分割后分别拼接的选择
        """
        运行IDBA软件，对序列进行组装
        :return:
        """
        # sample_name = os.path.basename(self.option('fastq1').prop['path']).split('.sickle.l.fastq')[0]
        if os.path.exists(self.work_dir + '/' + self.sample_name + '/contig.fa'):
            self.logger.info("拼接结果已存在，跳过拼接")
            return
        if not number:
            cmd = self.idba_path + 'idba_ud -r %s -o %s  --pre_correction --num_threads 16 --mink %s --maxk %s --step %s --min_contig %s' % (
                self.work_dir + '/' + self.sample_name + '.fa',
                self.work_dir + '/' + self.sample_name,
                self.option('mink'),
                self.option('maxk'),
                self.option('step'),
                self.option('min_contig'))
            self.logger.info('运行Idba,进行组装')
            command = self.add_command("idba_cmd", cmd, ignore_error=True).run()
            self.wait(command)
            if command.return_code == 0:
                self.logger.info("idba组装完成")
            elif command.return_code in [-9,-6,137]:  # 加入return_code检测，idba在sanger超出内存的返回值为-9,有时为-6
                self.logger.info("return code: %s" % command.return_code)
                self.add_state('memory_limit', 'memory is low!')   # add memory limit error by guhaidong @ 20180320
            else:
                self.logger.info("return code: %s" % command.return_code)
                self.set_error("idba组装出错!", code="31301005")
        else:
            contig_list = open(self.work_dir + "/contig.list", "w")
            for i in range(1, number + 1):
                contig_list.write(self.work_dir + '/' + self.sample_name + '_split' + str(i) + '/contig.fa\n')
                if os.path.exists(self.work_dir + '/' + self.sample_name + '_split' + str(i) + '/contig.fa'):
                    self.logger.info(self.sample_name + '_split' + str(i) + "拼接结果已有，跳过")
                    continue
                cmd = self.idba_path + 'idba_ud -r %s -o %s --pre_correction --num_threads 16 --mink %s --maxk %s --step %s --min_contig %s' % (
                    self.work_dir + '/split/' + self.sample_name + '_' + str(i) + '.fa',  # work_dir/split/sample_n.fa
                    self.work_dir + '/' + self.sample_name + '_split' + str(i),
                    self.option('mink'),
                    self.option('maxk'),
                    self.option('step'),
                    self.option('min_contig')
                )
                command = self.add_command("idba_cmd_{}".format(i), cmd, ignore_error=True).run()
                self.wait(command)
                if command.return_code == 0:
                    self.logger.info("idba组装完成")
                elif command.return_code in [-9,-6,137]:  # 加入return_code检测，idba在sanger超出内存的返回值为-9,有时为-6
                    self.logger.info("return code: %s" % command.return_code)
                    self.add_state('memory_limit', 'memory is low!')   # add memory limit error by guhaidong @ 20180320
                else:
                    self.logger.info("return code: %s" % command.return_code)
                    self.set_error("idba组装出错！", code="31301006")
            contig_list.close()

    def run_combine_contig(self):
        """
        将所有的contig文件合并在一起，注意修改序列序号
        :return:
        """
        # 输入contig路径表 输出至定义的文件，没有路径创造路径
        # cmd = combine_contig.pl contig.list self.work_dir + '/' + self.sample_name + '/contig.fa'
        if not os.path.exists(self.work_dir + '/' + self.sample_name):
            os.makedirs(self.work_dir + '/' + self.sample_name)
        '''
        for i in range(1, self.option('split_num') + 1):
            if os.path.exists(self.work_dir + '/combine/' + self.sample_name + '_split' + str(i) + '.contig.fa'):
                os.remove(self.work_dir + '/combine/' + self.sample_name + '_split' + str(i) + '.contig.fa')
            os.link(self.work_dir + '/' + self.sample_name + '_split' + str(i) + '/contig.fa',
                    self.work_dir + '/combine/' + self.sample_name + '_split' + str(i) + '.contig.fa')
        '''
        cmd = self.perl_path + self.combine_contig_path + " -list %s -o %s" % (self.work_dir + "/contig.list",
                                                                               self.work_dir + "/" + self.sample_name + "/contig.fa")
        self.logger.info('正在合并contig结果')
        command = self.add_command("combine_contig", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("合并contig结果完成")
        else:
            self.set_error("合并contig解果出错！", code="31301007")

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        #     sample_name = os.path.basename(self.option('fastq1').prop['path']).split('.sickle.l.fastq')[0]
        output_file = self.output_dir + '/' + self.sample_name + '.contig.fa'
        if os.path.exists(output_file):
            os.remove(output_file)
        os.link(self.work_dir + '/' + self.sample_name + '/contig.fa', output_file)
        self.option('contig').set_path(output_file)
        self.logger.info('设置组装拼接分析结果目录成功')
        #     if self.option('fastqs').is_set:
        #     shutil.copy2(self.work_dir + "/" + sample_name + "/transcripts.gtf", self.output_dir + "/transcripts.gtf")
        #     self.logger.info("设置组装拼接分析结果目录成功")

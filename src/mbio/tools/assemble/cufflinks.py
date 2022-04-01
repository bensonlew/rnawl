# -*- coding: utf-8 -*-
# __author__ = 'wangzhaoyue'
import os
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.config import Config
from mbio.packages.ref_rna.filter_gtf import FilterGtf
import re


class CufflinksAgent(Agent):
    """
    有参转录组拼接
    version v1.0.1
    author: wangzhaoyue
    last_modify: 2016.09.09
    """
    def __init__(self, parent):
        super(CufflinksAgent, self).__init__(parent)
        options = [
            {"name": "sample_bam", "type": "infile", "format": "align.bwa.bam"},  # 所有样本比对之后的bam文件
            {"name": "ref_fa", "type": "infile", "format": "sequence.fasta"},  # 参考基因文件
            {"name": "ref_gtf", "type": "infile", "format": "gene_structure.gtf"},  # 参考基因的注释文件
            {"name": "cpu", "type": "int", "default": 10},  #cufflinks软件所分配的cpu数量
            {"name": "F", "type": "int", "default": 0.1},  # min-isoform-fraction
            {"name": "fr_stranded", "type": "string", "default": "fr-unstranded"},  # 是否链特异性
            {"name": "strand_direct", "type": "string", "default": "none"},  # 链特异性时选择正负链
            {"name": "sample_gtf", "type": "outfile", "format": "gene_structure.gtf"},  # 输出的转录本文件
        ]
        self.add_option(options)
        self.step.add_steps("cufflinks")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.cufflinks.start()
        self.step.update()

    def stepfinish(self):
        self.step.cufflinks.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option('sample_bam'):
            raise OptionError('必须输入样本文件为bam格式')
        if not self.option('ref_fa'):
            raise OptionError('必须输入参考序列ref.fa')
        if not self.option('ref_gtf'):
            raise OptionError('必须输入参考序列ref.gtf')
        if self.option("fr_stranded") != "fr-unstranded" and not self.option("strand_direct").is_set:
            raise OptionError("当链特异性时必须选择正负链")
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 11
        self._memory = "15G"

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["_out.gtf", "gtf", "样本拼接之后的gtf文件"]
        ])
        super(CufflinksAgent, self).end()


class CufflinksTool(Tool):
    def __init__(self, config):
        super(CufflinksTool, self).__init__(config)
        self._version = "v1.0.1"
        self.cufflinks_path = '/bioinfo/rna/cufflinks-2.2.1/'
        self.bioawk_path =  Config().SOFTWARE_DIR + '/bioinfo/seq/bioawk/'
        self.bedtools_path = Config().SOFTWARE_DIR + '/bioinfo/seq/bedtools-2.25.0/bin/'
        self.script_path = '/bioinfo/rna/scripts/'

    def run(self):
        """
        运行
        :return:
        """
        super(CufflinksTool, self).run()
        self.run_cufflinks()
        self.run_gtf_to_fa()
        self.set_output()
        self.end()

    def run_cufflinks(self):
        """
        运行cufflinks软件，进行拼接组装
        """
        cmd = ""
        #  修改cufflinks参数，避免坐标超出范围时运行不出来 刘彬旭
        if self.option('fr_stranded') == "fr-unstranded":
            sample_name = os.path.basename(self.option('sample_bam').prop['path']).split('.bam')[0]
            cmd = self.cufflinks_path + ('cufflinks -p %s -g %s --library-type %s -o  ' % (
                 self.option('cpu'), self.option('ref_gtf').prop['path'], self.option('fr_stranded')) +
                 sample_name) + ' %s' % (self.option('sample_bam').prop['path'])
        else:
            if self.option('strand_direct') == "firststrand":
                sample_name = os.path.basename(self.option('sample_bam').prop['path']).split('.bam')[0]
                cmd = self.cufflinks_path + ('cufflinks -p %s -g %s --library-type %s -o  ' % (
                    self.option('cpu'), self.option('ref_gtf').prop['path'],
                    self.option('strand_direct')) + sample_name) + ' %s' % (self.option('sample_bam').prop['path'])
            else:
                if self.option('strand_direct')=='secondstrand':
                    sample_name = os.path.basename(self.option('sample_bam').prop['path']).split('.bam')[0]
                    cmd = self.cufflinks_path +( 'cufflinks -p %s -g %s --library-type %s -o  ' % (
                        self.option('cpu'), self.option('ref_gtf').prop['path'],
                        self.option('strand_direct')) + sample_name) + ' %s' % (self.option('sample_bam').prop['path'])

        self.logger.info('运行cufflinks软件，进行组装拼接')
        command = self.add_command("cufflinks_cmd", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("cufflinks运行完成")
        else:
            self.set_error("cufflinks运行出错!")

        # 增加过滤超出基因组范围坐标步骤 刘彬旭
        cmd = self.script_path + "fasta_range.sh %s %s %s" % (
            self.bioawk_path, self.option('ref_fa').prop['path'], self.work_dir + "/" + sample_name + '.filter.bed')
        self.logger.info('运行biowat，生成基因组区间')
        command = self.add_command("ref_fa_to_bed_cmd", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("bioawk运行完成")
        else:
            self.set_error("bioawk运行出错!")

        cmd = self.script_path + "filter_gtf_by_range.sh %s %s %s %s" % (
            self.bedtools_path , self.work_dir + "/" + sample_name + '.filter.bed', self.work_dir + "/" + sample_name + '/transcripts.gtf',self.work_dir + "/" + sample_name + '/transcripts.filter.gtf')
        self.logger.info('运行bedtools，过滤基因组区间')
        command = self.add_command("bedtools_filter_cmd", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("bedtools运行完成")
            if os.path.getsize(self.work_dir + "/" + sample_name + '/transcripts.filter.gtf') == 0 :
                self.logger.info('bedtools过滤结果为空值, 使用package过滤')
                try:
                    FilterGtf().filter_gtf_by_bed(self.work_dir + "/" + sample_name + '/transcripts.gtf',
                                                  self.work_dir + "/" + sample_name + '/transcripts.filter.gtf',
                                                  self.work_dir + "/" + sample_name + '.filter.bed')
                except:
                    self.set_error("过滤gtf运行出错!")

        else:
            self.logger.info("bedtools运行出错")

    def run_gtf_to_fa(self):
        """
        运行gtf_to_fasta，转录本gtf文件转fa文件
        """
        sample_name = os.path.basename(self.option('sample_bam').prop['path']).split('.bam')[0]
        cmd = self.cufflinks_path + "gffread %s -g %s -w %s_out.fa" % (
            self.work_dir + "/" + sample_name + "/transcripts.filter.gtf", self.option('ref_fa').prop['path'],
        self.work_dir + "/" + sample_name)
        self.logger.info('运行gtf_to_fasta，形成fasta文件')
        command = self.add_command("gtf_to_fa_cmd", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("gtf_to_fasta运行完成")
        else:
            self.set_error("gtf_to_fasta运行出错!")

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        try:
            sample_name = os.path.basename(self.option('sample_bam').prop['path']).split('.bam')[0]
            if os.path.exists(self.output_dir + "/transcripts.gtf"):
                os.remove(self.output_dir + "/transcripts.gtf")
            if os.path.exists(self.output_dir + "/" + sample_name + "_out.fa"):
                os.remove(self.output_dir + "/" + sample_name + "_out.fa")
            os.link(self.work_dir + "/" + sample_name + "/transcripts.filter.gtf", self.output_dir + "/transcripts.gtf")
            shutil.move(self.output_dir + "/transcripts.gtf", self.output_dir + "/" + sample_name + "_out.gtf")
            os.link(self.work_dir + "/" + sample_name + "_out.fa", self.output_dir + "/" + sample_name + "_out.fa")
            self.option('sample_gtf').set_path(self.work_dir + "/" + sample_name + "/" + "transcripts.filter.gtf")
            self.logger.info("设置组装拼接分析结果目录成功")

        except Exception as e:
            self.logger.info("设置组装拼接分析结果目录失败{}".format(e))
            self.set_error("设置组装拼接分析结果目录失败{}".format(e))

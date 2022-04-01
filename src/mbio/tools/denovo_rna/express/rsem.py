# -*- coding: utf-8 -*-
# __author__ = 'qiuping'
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re


class RsemAgent(Agent):
    """
    调用align_and_estimate_abundance.pl脚本，运行rsem，进行表达量计算分析
    version v1.0
    author: qiuping
    last_modify: 2016.06.20
    """
    def __init__(self, parent):
        super(RsemAgent, self).__init__(parent)
        options = [
            {"name": "fq_type", "type": "string"},  # PE OR SE
            {"name": "rsem_fa", "type": "infile", "format": "sequence.fasta"},  # trinit.fasta文件
            {"name": "fq_l", "type": "infile", "format": "sequence.fastq"},  # PE测序，包含所有样本的左端fq文件的文件
            {"name": "fq_r", "type": "infile", "format": "sequence.fastq"},  # PE测序，包含所有样本的左端fq文件的文件
            {"name": "fq_s", "type": "infile", "format": "sequence.fastq"},  # SE测序，包含所有样本的fq文件的文件
            {"name": "fa_build", "type": "outfile", "format": "sequence.fasta"},  # trinit.fasta文件
            {"name": "only_bowtie_build", "type": "bool", "default": False},  # 为true时该tool只建索引
            {"name": "bowtie_build_rsem", "type": "bool", "default": False}  # 为true时,建完索引后运行rsem，否则默认有bowtie2的索引，只跑rsem
        ]
        self.add_option(options)
        self.step.add_steps("rsem")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.rsem.start()
        self.step.update()

    def stepfinish(self):
        self.step.rsem.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option('fq_type'):
            raise OptionError('必须设置测序类型：PE OR SE')
        if self.option('fq_type') not in ['PE', 'SE']:
            raise OptionError('测序类型不在所给范围内')
        if not self.option('only_bowtie_build'):
            if self.option("fq_type") == "PE" and not self.option("fq_r").is_set and not self.option("fq_l").is_set:
                raise OptionError("PE测序时需设置左端序列和右端序列输入文件")
            if self.option("fq_type") == "SE" and not self.option("fq_s").is_set:
                raise OptionError("SE测序时需设置序列输入文件")
        if not self.option("rsem_fa").is_set:
            raise OptionError("需设置rsem_fa输入文件")
        if not isinstance(self.option('only_bowtie_build'), bool):
            raise OptionError('only_bowtie_build只能为bool')
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 10
        self._memory = '20G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"]
        ])
        result_dir.add_regexp_rules([
            [r"results$", "xls", "rsem结果"]
        ])
        super(RsemAgent, self).end()


class RsemTool(Tool):
    """
    Lefse tool
    """
    def __init__(self, config):
        super(RsemTool, self).__init__(config)
        self._version = '1.0.1'
        self.rsem = "/bioinfo/rna/scripts/align_and_estimate_abundance.pl"
        self.rsem_path = self.config.SOFTWARE_DIR + '/bioinfo/rna/RSEM-1.2.31/bin'
        self.bowtie_path = self.config.SOFTWARE_DIR + '/bioinfo/align/bowtie2-2.2.9/'
        # self.bowtie2 = '/bioinfo/align/bowtie2-2.2.9/bowtie2-build'
        self.gcc = self.config.SOFTWARE_DIR + '/gcc/5.1.0/bin'
        self.gcc_lib = self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64'
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)
        self.set_environ(PATH=self.rsem_path)
        self.set_environ(PATH=self.bowtie_path)
        # self.rsem_fasta = self.work_dir + '/' + os.path.basename(self.option('rsem_fa').prop['path'])

    def bowtie_build(self):
        rsem_fasta = self.work_dir + '/' + os.path.basename(self.option('rsem_fa').prop['path'])
        if os.path.exists(rsem_fasta):
            os.remove(rsem_fasta)
            os.link(self.option('rsem_fa').prop['path'], rsem_fasta)
        else:
            os.link(self.option('rsem_fa').prop['path'], rsem_fasta)
        cmd = self.rsem + ' --transcripts %s --seqType fq --single test.fq --est_method  RSEM --output_dir %s --trinity_mode --aln_method bowtie2 --prep_reference' % (rsem_fasta, self.work_dir)
        self.logger.info('开始运行bowtie2建索引')
        bowtie_cmd = self.add_command('bowtie_build', cmd).run()
        self.wait()
        if bowtie_cmd.return_code == 0:
            self.logger.info("%s运行完成" % bowtie_cmd)
            self.option('fa_build', rsem_fasta)
        else:
            self.set_error("%s运行出错!" % bowtie_cmd)
            raise("%s运行出错!" % bowtie_cmd)

    def run_rsem(self, rsem_fasta):
        if self.option('fq_type') == 'SE':
            sample = os.path.basename(self.option('fq_s').prop['path']).split('_sickle_s.fastq')[0]
            rsem_cmd = self.rsem + ' --transcripts %s --seqType fq --single %s --est_method  RSEM --output_dir %s --thread_count 6 --trinity_mode --aln_method bowtie2 --output_prefix %s' % (rsem_fasta, self.option('fq_s').prop['path'], self.work_dir, sample)
        else:
            sample = os.path.basename(self.option('fq_l').prop['path']).split('_sickle_l.fastq')[0]
            rsem_cmd = self.rsem + ' --transcripts %s --seqType fq --right %s --left %s --est_method  RSEM --output_dir %s --thread_count 6 --trinity_mode --aln_method bowtie2 --output_prefix %s' % (rsem_fasta, self.option('fq_r').prop['path'], self.option('fq_l').prop['path'], self.work_dir, sample)

        self.logger.info("开始运行_rsem_cmd")
        cmd = self.add_command("rsem_cmd", rsem_cmd).run()
        self.wait()
        self.logger.info("....返回码为：%s" % cmd.return_code)
        if cmd.return_code == 0:
            self.logger.info("%s运行完成" % cmd)
        else:
            self.set_error("%s运行出错!" % cmd)
            raise("%s运行出错!" % cmd)

    def set_output(self):
        """
        将结果文件link到output文件夹下面
        :return:
        """
        for root, dirs, files in os.walk(self.output_dir):
            for names in files:
                os.remove(os.path.join(root, names))
        self.logger.info("设置结果目录")
        results = os.listdir(self.work_dir)
        try:
            for f in results:
                if re.search(r'results$', f):
                    os.link(self.work_dir + '/' + f, self.output_dir + '/' + f)
            self.logger.info("设置rsem分析结果目录成功")
        except Exception as e:
            self.logger.info("设置rsem分析结果目录失败{}".format(e))

    def run(self):
        super(RsemTool, self).run()
        if self.option('only_bowtie_build'):
            self.bowtie_build()
        else:
            if self.option('bowtie_build_rsem'):
                self.bowtie_build()
                self.run_rsem(self.option('fa_build').prop['path'])
            else:
                self.run_rsem(self.option('rsem_fa').prop['path'])
        self.set_output()
        self.end()

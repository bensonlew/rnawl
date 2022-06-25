# -*- coding: utf-8 -*-
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
import types
from biocluster.core.exceptions import OptionError
from mbio.packages.beta_diversity.plsda_r import plsda
import subprocess


class CodonBiasAgent(Agent):
    """

    """
    def __init__(self, parent):
        super(CodonBiasAgent, self).__init__(parent)
        options = [
            {"name": "fasta", "type": "infile", "format": "sequence.fasta"}
        ]
        self.add_option(options)


    def check_options(self):
        """
        重写参数检查
        """
        if self.option('fasta').is_set:
            pass
        else:
            raise OptionError('没有提供fasta')
        return True

    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 2
        self._memory = '3G'

    def end(self):
        super(CodonBiasAgent, self).end()


class CodonBiasTool(Tool):
    def __init__(self, config):
        super(CodonBiasTool, self).__init__(config)
        self._version = '1.0'

        self.codonw = self.config.SOFTWARE_DIR + '/bioinfo/tool_lab/CodonW/codonw'
        self.perl_src = self.config.PACKAGE_DIR + '/tool_lab/code.pl'
        #self.r_bin = self.config.SOFTWARE_DIR + '/program/R-3.3.1/bin'
        self.perl_bin = self.config.SOFTWARE_DIR + '/miniconda2/bin/perl'
        self.set_environ(LD_LIBRARY_PATH="{}/gcc/5.1.0/lib64".format(self.config.SOFTWARE_DIR))


    def run(self):
        """
        运行
        """
        super(CodonBiasTool, self).run()
        self.run_codon()

    def run_codon(self):
        if os.path.exists('fasta.result'):
            os.remove('fasta.result')

        if os.path.exists('fasta.bulkoutfile'):
            os.remove('fasta.bulkoutfile')

        in_fasta = self.option("fasta").path
        in_name = os.path.basename(in_fasta)
        if os.path.exists(in_name):
            os.remove(in_name)
        os.link(in_fasta, in_name)

        cmd = self.codonw + " %s  fasta.result fasta.bulkoutfile -all_indices -nomenu" %(in_name)
        try:
            self.logger.info(cmd)
            subprocess.check_output(cmd, shell=True)
        except  Exception as e:
            self.set_error('运行codon失败')

        cmd2 = self.perl_bin +' '+ self.perl_src + ' -out %s' % self.work_dir
        try:
            self.logger.info(cmd2)
            subprocess.check_output(cmd2, shell=True)
        except  Exception as e:
            self.set_error('生成R脚本失败')


        cmd3 = self.config.SOFTWARE_DIR + '/program/R-3.3.1_gcc5.1/bin/Rscript %s/codon.r' % self.work_dir
        self.logger.info(cmd3)
        try:
            subprocess.check_output(cmd3, shell=True)
        except Exception as e:
            self.set_error("运行R脚本失败")


        # r_cmd = self.add_command('run_r',cmd3).run()
        # if r_cmd.return_code == 0:
        #     self.logger.info('%s 运行成功'%cmd3)
        # else:
        #     self.set_error('%s 运行失败'%cmd3)

        self.set_output()
        self.end()

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        for f in ['result.txt', 'fasta.result','Sample_codon.png']:
            out = self.output_dir + '/' + f
            if os.path.exists(out):
                os.remove(out)
            os.link(self.work_dir+'/'+ f, out)




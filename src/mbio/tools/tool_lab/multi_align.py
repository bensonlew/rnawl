# -*- coding: utf-8 -*-

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
# import copy
import collections
import pandas as pd
import subprocess
import math
import shutil


class MultiAlignAgent(Agent):

    def __init__(self, parent):
        super(MultiAlignAgent, self).__init__(parent)
        options = [
            {"name": "fasta", "type": "infile", "format": "sequence.fasta"},
            {"name": "align_method", "type": "string", "default": "muscle"},  #muscle ,mafft, clustalo
            {"name": "out_format", "type": "string", "default": "fasta"}
        ]
        self.add_option(options)

    def check_options(self):
        if not self.option('fasta').is_set:
            raise OptionError('必须提供输入序列文件')
        self.option('fasta').get_info()
        if int(self.option('fasta').prop['seq_number']) < 2:
            raise OptionError('序列条数必须大于等于2')
        #if int(self.option('fasta').prop['seq_number']) > 300:
        #    raise OptionError('%s序列条数必须小于等于300', variables=(self.option('fasta').prop['seq_number'],), code="32300603")
        # if self.option('sequence_type') not in ["coding", "no_coding", "amino_acid"]:
        #     raise OptionError('序列类型只能传入coding、no_coding、amino_acid', code="32300606")
        #seq_name = self.option('fasta').get_all_seq_name()
        if self.option('align_method') == 'mafft':
            #self.out_format = 'fasta'
            if self.option('out_format')  in ['fasta','clustal']:
                #self.out_format = self.option('out_format')
                pass
            else:
                raise OptionError('mafft out format: fasta,clustal')
        elif self.option('align_method') == 'muscle':
            if self.option('out_format')  in ['fasta','msf','clustal']:
                #self.out_format = self.option('out_format')
                pass
            else:
                raise OptionError('muscle out format: fasta,msf,clustal')
        elif self.option('align_method') == 'clustalo':
            clustal_out_f = ['fasta','msf','clustal','phylip','selex','stockholm','vienna']
            if self.option('out_format')  in clustal_out_f:
                #a2m=fa[sta],clu[stal],msf,phy[lip],selex,st[ockholm],vie[nna]}
                #self.out_format = self.option('out_format')
                pass
            else:
                raise OptionError('clustalo out format: %s'%(','.join(clustal_out_f)))



    def set_resource(self):
        """
        设置所需资源
        """
        self._cpu = 3
        total = os.path.getsize(self.option("fasta").prop["path"])
        total = int(math.ceil(total / (1024 * 1024 * 1024)))
        total = 10 + int(total * 45)
        #total = '40G'
        self._memory = "{}G".format(total)


    def end(self):

        super(MultiAlignAgent, self).end()


class MultiAlignTool(Tool):
    def __init__(self, config):
        super(MultiAlignTool, self).__init__(config)
        self.muscle = '/bioinfo/phylogenetic/muscle-3.8.31-release/muscle'
        self.mafft_path = self.config.SOFTWARE_DIR+'/bioinfo/align/mafft-7.299-with-extensions/bin/mafft'
        self.clustalo = self.config.SOFTWARE_DIR + '/bioinfo/tool_lab/clustal-omega-1.2.4/bin/clustalo'

        self.trimal ='/bioinfo/phylogenetic/trimal-trimAl/source/trimal'


    def run(self):
        """
        运行
        """
        super(MultiAlignTool, self).run()
        self.run_tool()

        self.end()



    def run_tool(self):
        self.out_format = self.option('out_format')

        fasta_file = self.option('fasta').prop['path']
        cmd= "sed -i 's/\s.*//g' "+fasta_file
        try :
            subprocess.check_output(cmd,shell=True)
            self.logger.info("fasta文件修改成功")
        except subprocess.CalledProcessError:
            self.set_error("fasta文件生成失败", code="32300603")

        align = self.output_dir + '/align'
        if self.option("align_method") in ["mafft"]:  #add 20191015
            if self.out_format == 'fasta':
                align = align + '.fasta'
                cmd = "{} {} > {}".format(self.mafft_path, fasta_file, align)
            else:
                align = align + '.clustal'
                cmd = "{} {} --clustalout > {}".format(self.mafft_path, fasta_file, align)
            self.logger.info("开始运行{}，进行比对".format(cmd))
            command = subprocess.Popen(cmd, shell=True)
            command.communicate()
            if command.returncode == 0:
                self.logger.info("完成比对！")
            else:
                self.set_error("mafft运行出错！", code="32300610")

        elif self.option("align_method") in ["muscle"]:
            if self.out_format == 'msf':
                align = align+'.msf'
                cmd1 = '%s -in %s -out %s -msf' % (self.muscle, fasta_file, align)
            elif self.out_format == 'clustal':
                align = align+'.clustal'
                cmd1 = '%s -in %s -out %s -clw' % (self.muscle, fasta_file, align)
            else:
                align = align+'.fasta'
                cmd1 = '%s -in %s -out %s ' % (self.muscle, fasta_file, align)

            command1 = self.add_command('cmd_1', cmd1)
            command1.run()
            self.wait(command1)
            if command1.return_code == 0:
                self.logger.info("align succeed")
            else:
                self.set_error("align failed")

        else:
            align = align +'.'+self.out_format
            cmd = '%s -i %s -o %s --outfmt=%s --force'%(self.clustalo, fasta_file, align, self.out_format)
            self.logger.info(cmd)
            #command = self.add_command('cmd_2', cmd)
            #command.run()
            #if command.return_code == 0:
            command = subprocess.Popen(cmd, shell=True)
            command.communicate()
            if command.returncode == 0:
                self.logger.info("align succeed")
            else:
                self.set_error("align failed")

        if self.out_format in ['fasta', 'clustal']:  #'phylip'
            trimal_out = self.output_dir + '/trimal.'+self.out_format
            cmd_trimal = '%s -in %s -out %s'%(self.trimal, align, trimal_out)
            command1 = self.add_command('cmd_trimal', cmd_trimal)
            command1.run()
            self.wait(command1)
            if command1.return_code == 0:
                self.logger.info("trimal succeed")
            else:
                self.set_error("trimal failed")

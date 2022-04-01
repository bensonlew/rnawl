# -*- coding: utf-8 -*-
# __author__ = 'zhouxuan'
import os
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool
import subprocess
import re


class SoapfuseAgent(Agent):
    """
    基于RNA-data的基因融合分析SOAPfuse
    version v1
    author: zhouxuan 
    last_modify: 2016.11.20  2017.1.16
    """
    def __init__(self, parent):
        super(SoapfuseAgent, self).__init__(parent)
        options = [
            {"name": "sample_data", "type": "infile", "format": "ref_rna.gene_fusion.sample_data_dir"},  # 存放样本数据的文件夹
            {"name": "sample_list", "type": "infile", "format": "ref_rna.gene_fusion.sample_list"},
            # 样本信息文件(其中只包含一个样本的信息)
            {"name": "sample_name", "type": "string"},  # 所进行分析的样本名称
            {"name": "reads_number", "type": "int", "default": 3},  # 支持融合的reads数目的最小值 3 5 10
            {"name": "fusions_otu", "type": "outfile", "format": "ref_rna.gene_fusion.fusionsout_dir"}
            #  基因融合的单个样本的结果文件夹
        ]
        self.add_option(options)
        self.step.add_steps("soapfuse")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.soapfuse.start()
        self.step.update()

    def stepfinish(self):
        self.step.soapfuse.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        self.logger.info(self.option('sample_name'))
        if not self.option('sample_data'):
            raise OptionError('必须输入样本数据文件夹，文件夹里的文件为fq.gz格式')
        if not self.option('sample_list') :
            raise OptionError('必须输入样本信息文件，文件里的信息必须完整无误')
        if not self.option('sample_name'):
            raise OptionError('必须选择要分析的样本')
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 10
        self._memory = "100G"

    # def end(self):  #tool 中暂时不设置该处上传目录的问题
    #     result_dir = self.add_upload_dir(self.output_dir)
    #     result_dir.add_relpath_rules([
    #         [".", "", "结果输出目录"],
    #     ])
    #     result_dir.add_regexp_rules([
    #        # ["_out.gtf", "gtf", "样本拼接之后的gtf文件"]
    #     ])
    #
    #     super(SoapfuseAgent, self).end()


class SoapfuseTool(Tool):
    def __init__(self, config):
        super(SoapfuseTool, self).__init__(config)
        self._version = "v1.26"
        self.soapfuse_path = '/mnt/ilustre/users/sanger-dev/app/bioinfo/rna/SOAPfuse-v1.26/'
        self.perl_path = 'program/perl-5.24.0/bin/perl '

    def run(self):
        """
        运行
        :return:
        """
        super(SoapfuseTool, self).run()
        self.run_soapfuse()
        self.set_output()
        # self.end() # 和def end() 函数一起使用

    def run_soapfuse(self):
        """
        运行soapfuse软件，为了使得分析结果能尽快得出，此处命令只能进行单个样本分析，多样本在module中设置
        """
        # cmd = self.perl_path+self.soapfuse_path + ('SOAPfuse-RUN.pl -c /mnt/ilustre/users/sanger-dev/app/bioinfo'
        #                                            '/rna/SOAPfuse-v1.26/config/config-%s.txt '
        #                                            '-fd %s -l %s -o %s' % (self.option('reads_number'),
        #                                                                    self.option('sample_data').prop['path'],
        #                                                                    self.option('sample_list').prop['path'],
        #                                                                    self.work_dir + '/' +
        #                                                                    str(self.option('sample_name'))))
        # cmd = self.perl_path + self.soapfuse_path + ('SOAPfuse-RUN.pl -c /mnt/ilustre/users/sanger-dev/app/bioinfo'
        #                                              '/rna/SOAPfuse-v1.26/config/config-%s.txt '
        #                                              '-fd %s -l %s -o %s' % (self.option('reads_number'),
        #                                                                      self.option('sample_data'),
        #                                                                      self.option('sample_list'),
        #                                                                      self.work_dir + '/' +
        #                                                                      str(self.option('sample_name'))))
        cmd = self.perl_path + self.soapfuse_path + ('SOAPfuse-RUN.pl -c /mnt/ilustre/users/sanger-dev/app/bioinfo'
                                                     '/rna/SOAPfuse-v1.26/config/config-%s.txt '
                                                     '-fd %s -l %s -o %s' % (self.option('reads_number'),
                                                                             self.option('sample_data'),
                                                                             self.option('sample_list'),
                                                                             self.work_dir + '/' +
                                                                             str(self.option('sample_name'))))
        self.logger.info('运行soapfuse软件,进行基因融合的分析')
        command = self.add_command("soapfuse_cmd", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("soapfuse运行完成!")
            result_dir = self.work_dir + '/' + self.option('sample_name')
            if os.path.exists(result_dir + '/TEMP'):
                the_config_list = os.listdir(result_dir + '/TEMP')
                number_list = []
                for i in the_config_list:
                   if os.path.exists(result_dir + '/TEMP' + '/' + i):
                       number_list.append(len(os.listdir(result_dir + '/TEMP' + '/' + i)))
                number_list.sort()
                if number_list[-1] != 19:
                    self.logger.info("此样本在此分析参数下不存在融合位点！")
        else:
            self.set_error("soapfuse运行出错!")

    def set_output(self):
            """
            将结果文件复制到output文件夹下面
            :return:
            """
            self.logger.info("设置结果目录")
            try:
                shutil.copytree(self.work_dir + '/' + self.option('sample_name') + '/final_fusion_genes/' + self.option(
                    'sample_name') + '/analysis/figures/fusions_expression', self.output_dir + "/pdf_result")
                shutil.copy(self.work_dir + "/" + self.option('sample_name') + "/final_fusion_genes/" + self.option(
                    'sample_name') + "/" + self.option('sample_name') + ".final.Fusion.specific.for.trans",
                            self.output_dir + "/final.Fusion.specific.for.trans")
                shutil.copy(self.work_dir + "/" + self.option('sample_name') + "/final_fusion_genes/" + self.option(
                    'sample_name') + "/" + self.option('sample_name') + ".final.Fusion.specific.for.genes",
                            self.output_dir + "/final.Fusion.specific.for.genes")
                self.option('fusions_otu').set_path(self.work_dir + "/OUT/final_fusion_genes/")
                self.logger.info("设置基因融合分析结果目录成功")
            except Exception as e:
                self.logger.info("设置基因融合分析结果目录失败{}".format(e))
                self.set_error("设置基因融合分析结果目录失败{}".format(e))
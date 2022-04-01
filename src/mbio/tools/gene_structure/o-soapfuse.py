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
    last_modify: 2016.10.13
    """
    def __init__(self, parent):
        super(SoapfuseAgent, self).__init__(parent)
        options = [
            {"name": "sample_data", "type": "infile", "format": "ref_rna.sample_data_dir"},  # 存放样本数据的文件夹
            {"name": "sample_list", "type": "infile", "format": "ref_rna.sample_list_dir"},  # 样本信息文件夹
			
            {"name": "sample_name", "type": "string"},  # 所进行分析的样本名称
			{"name": "reads_number", "type": "int", "default" : 3}, # 支持融合的reads数目的最小值

            {"name": "fusionsotu", "type": "outfile","format":"ref_rna.fusionsout_dir"}# 基因水平上的融合位点详细信息
            # {"name": "finalfusion_for_trans", "type": "outfile", "format": "ref_rna.finalfusion_for_trans"},  # 转录本水平上的融合位点详细信息
            # {"name": "finalfusion_expression", "type": "outfile", "format": "ref_rna.finalfusion_expression_dir"},  # gene_fusion_expression
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
        if not self.option('sample_data'):
            raise OptionError('必须输入样本数据文件夹，文件夹里的文件为fq.gz格式')
        if not self.option('sample_list') :
            raise OptionError('必须输入样本信息文件夹，文件夹里的信息必须完整无误')
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

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
           # ["_out.gtf", "gtf", "样本拼接之后的gtf文件"]
		  
        ])
        super(SoapfuseAgent, self).end()


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
        self.end()

    def run_soapfuse(self):
        """
        运行soapfuse软件，基因融合的分析
        """
        if self.option('reads_number') == 3:
            cmd = self.perl_path+self.soapfuse_path + ('SOAPfuse-RUN.pl -c /mnt/ilustre/users/sanger-dev/app/bioinfo/rna/SOAPfuse-v1.26/config/config-1.txt -fd %s -l %s -o %s' %(self.option('sample_data').prop['path'],self.option('sample_list').prop['path']+'/'+self.option('sample_name'),self.work_dir + '/'+self.option('sample_name')))
        elif self.option('reads_number') == 5:
            cmd = self.perl_path+self.soapfuse_path + ('SOAPfuse-RUN.pl -c /mnt/ilustre/users/sanger-dev/app/bioinfo/rna/SOAPfuse-v1.26/config/config-2.txt -fd %s -l %s -o %s' %(self.option('sample_data').prop['path'],self.option('sample_list').prop['path']+'/'+self.option('sample_name'),self.work_dir + '/'+self.option('sample_name')))
        else:
            cmd = self.perl_path+self.soapfuse_path + ('SOAPfuse-RUN.pl -c /mnt/ilustre/users/sanger-dev/app/bioinfo/rna/SOAPfuse-v1.26/config/config-3.txt -fd %s -l %s -o %s' %(self.option('sample_data').prop['path'],self.option('sample_list').prop['path']+'/'+self.option('sample_name'),self.work_dir + '/'+self.option('sample_name')))
        self.logger.info('运行soapfuse软件,进行基因融合的分析')
        command = self.add_command("soapfuse_cmd", cmd).run()
        self.wait(command)
        if os.path.exists(self.work_dir + "/" + self.option('sample_name') + "/final_fusion_genes/" + self.option('sample_name') + "/" + self.option('sample_name') + ".final.Fusion.specific.for.genes"):
            if os.path.exists(self.work_dir + "/" + self.option('sample_name') + "/final_fusion_genes/" + self.option('sample_name') + "/" + self.option('sample_name') + ".final.Fusion.specific.for.trans"):
                if os.path.exists(self.work_dir + "/" + self.option('sample_name') + "/final_fusion_genes/" + self.option('sample_name') + "/analysis/figures/fusions_expression"):
                    self.logger.info("soapfuse运行成功!")
        else:
            self.set_error("soapfuse运行出错!")
        #if command.return_code == 0:
            #self.logger.info("soapfuse运行完成!")
        #else:
            #self.set_error("soapfuse运行出错!")

    def set_output(self):
            """
            将结果文件复制到output文件夹下面
            :return:
            """
            self.logger.info("设置结果目录")
            try:
                shutil.copy2(self.work_dir + "/" + self.option('sample_name') + "/final_fusion_genes/" + self.option('sample_name') + "/" + self.option('sample_name') + ".final.Fusion.specific.for.genes" ,self.output_dir + "/final.Fusion.specific.for.genes")
                shutil.copy2(self.work_dir + "/" + self.option('sample_name') + "/final_fusion_genes/" + self.option('sample_name') + "/" + self.option('sample_name') + ".final.Fusion.specific.for.trans",self.output_dir + "/final.Fusion.specific.for.trans")
                shutil.copytree(self.work_dir + "/" + self.option('sample_name') + "/final_fusion_genes/" + self.option('sample_name') + "/analysis/figures/fusions_expression" , self.output_dir + "/fusions_expression")
                self.logger.info("设置基因融合分析结果目录成功")

            except Exception as e:
                self.logger.info("设置基因融合分析结果目录失败{}".format(e))
                self.set_error("设置基因融合分析结果目录失败{}".format(e))



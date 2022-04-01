# -*- coding: utf-8 -*-
# __author__ = "qiuping"
# last_modify:20161031

from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import re


class NetworkAgent(Agent):
    """
    调用wgcna脚本，进行网络共表达计算分析
    version v1.0
    author: qiuping
    last_modify: 2016.07.12
    """
    def __init__(self, parent):
        super(NetworkAgent, self).__init__(parent)
        options = [
            {"name": "diff_fpkm", "type": "infile", "format": "rna.express_matrix"},  # 输入文件，差异基因表达量矩阵
            {"name": "gene_file", "type": "infile", "format": "rna.gene_list"},  # 差异基因名称文件
            {"name": "softpower", "type": "int", "default": 9},
            {"name": "dissimilarity", "type": "float", "default": 0.25},
            {"name": "module", "type": "float", "default": 0.1},
            {"name": "network", "type": "float", "default": 0.2}
        ]
        self.add_option(options)
        self.step.add_steps("network")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.network.start()
        self.step.update()

    def stepfinish(self):
        self.step.network.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option("diff_fpkm").is_set:
            raise OptionError("必须设置输入文件:差异基因fpkm表")
        if self.option("softpower") > 20 or self.option("softpower") < 1:
            raise OptionError("softpower值超出范围")
        if self.option('dissimilarity') > 1 or self.option("dissimilarity") < 0:
            raise OptionError("模块dissimilarity相异值超出范围")
        if not self.option("gene_file").is_set:
            raise OptionError("必须设置输入文件:基因名字列表")
        if self.option('module') > 1 or self.option("module") < 0:
            raise OptionError("模块module相异值超出范围")
        if self.option('network') > 1 or self.option("network") < 0:
            raise OptionError("模块network相异值超出范围")
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 10
        self._memory = '5G'

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
            ["all_edges.txt", "txt", "edges结果信息"],
            ["all_nodes.txt ", "txt", "nodes结果信息"],
            ["removeGene.xls ", "xls", "移除的基因信息"],
            ["removeSample.xls ", "xls", "移除的样本信息"],
            ["softPower.pdf", "pdf", "softpower相关信息"],
            ["ModuleTree.pdf", "pdf", "ModuleTree图"],
            ["eigengeneClustering.pdf", "pdf", "eigengeneClustering图"],
            ["eigenGeneHeatmap.pdf", "pdf", "eigenGeneHeatmap图"],
            ["networkHeatmap.pdf", "pdf", "networkHeatmap图"],
            ["sampleClustering.pdf", "pdf", "sampleClustering图"]
        ])
        result_dir.add_regexp_rules([
            [r"^CytoscapeInput.*", "txt", "Cytoscape作图数据"]
        ])
        super(NetworkAgent, self).end()


class NetworkTool(Tool):
    """
    表达量差异检测tool
    """
    def __init__(self, config):
        super(NetworkTool, self).__init__(config)
        self._version = '1.0.1'
        self.r_path = '/program/R-3.3.1/bin/Rscript'
        self.script_path = self.config.SOFTWARE_DIR + '/bioinfo/rna/scripts/'
        self.gcc = self.config.SOFTWARE_DIR + '/gcc/5.1.0/bin'
        self.gcc_lib = self.config.SOFTWARE_DIR + '/gcc/5.1.0/lib64'
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)

    def run_wgcna_one(self):
        one_cmd = self.r_path + " %sInModuleWGCNA-step01.r --args %s %s %s %s" % (self.script_path, self.option('diff_fpkm').prop['path'], 'wgcna_result', self.option('softpower'), self.option('dissimilarity'))
        self.logger.info(one_cmd)
        self.logger.info("开始运行one_cmd")
        cmd = self.add_command("one_cmd", one_cmd).run()
        self.wait(cmd)
        if cmd.return_code == 0:
            self.logger.info("运行one_cmd成功")
        else:
            self.logger.info("运行one_cmd出错")

    def run_wgcna_two(self):
        gene_file_path = self.work_dir + '/gene_file'
        self.option('gene_file').get_network_gene_file(gene_file_path)
        two_cmd = self.r_path + " %sInModuleWGCNA-step02.r --args %s %s %s %s" % (self.script_path, 'wgcna_result', gene_file_path, self.option('module'), self.option('network'))
        self.logger.info("开始运行two_cmd")
        cmd = self.add_command("two_cmd", two_cmd).run()
        self.wait(cmd)
        if cmd.return_code == 0:
            self.logger.info("运行two_cmd成功")
        else:
            self.logger.info("运行two_cmd出错")

    def convert_pdf_to_png(self, olds, news):
        self.image_magick = '/program/ImageMagick/bin/convert'
        convert_commands = []
        for index, i in enumerate(olds):
             cmd = self.image_magick + ' -flatten -quality 100 -density 130 -background white ' + i + ' ' + news[index]
             command = self.add_command('convert_{}'.format(index), cmd)
             command.run()
             convert_commands.append(command)
             self.wait()
             for i in convert_commands:
                  if i.return_code == 0:
                      pass
                  else:
                      return False
        return True
 
    def set_output(self):
        """
        将结果文件link到output文件夹下面
        :return:
        """
        for root, dirs, files in os.walk(self.output_dir):
            for names in files:
                os.remove(os.path.join(root, names))
        self.logger.info("设置结果目录")
        results = os.listdir(self.work_dir + '/wgcna_result/')
        try:
            for f in results:
                if re.search(r'.*Rdata$', f):
                    pass
                else:
                    os.link(self.work_dir + '/wgcna_result/' + f, self.output_dir + '/' + f)
            self.logger.info('设置文件夹路径成功')
        except Exception as e:
            self.logger.info("设置network分析结果目录失败{}".format(e))
        pdf_path = [self.output_dir+"/softPower.pdf", self.output_dir+"/ModuleTree.pdf"]
        png_path = [self.output_dir+"/softPower.png", self.output_dir+"/ModuleTree.png"]
        info = self.convert_pdf_to_png(pdf_path, png_path)
        if info:
            self.logger.info("pdf转png文件成功！")
        else:
            raise Exception("pdf转png文件失败！")
        
        

    def run(self):
        super(NetworkTool, self).run()
        self.run_wgcna_one()
        self.run_wgcna_two()
        self.set_output()
        self.end()

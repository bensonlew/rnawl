#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division
import math
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
import subprocess
from biocluster.core.exceptions import OptionError


class PaupPhylipAgent(Agent):
    """
    version 1.0
    """

    def __init__(self, parent):
        super(PaupPhylipAgent, self).__init__(parent)
        options = [
            {"name": "fasta_file", "type": "infile", "format": "sequence.fasta"},  # 输入文件
            {"name": "phylo_tre", "type": "outfile", "format": "meta.beta_diversity.newick_tree"},  # 输出结果
            {"name": "method", "type": "string", "default": "mafft"},  # 比对方法
            {"name": "bootstrap", "type": "int", "default": 500},
            {"name": "tree_software", "type": "string", "default": "paup"}  #只做mp方法建树
        ]
        self.add_option(options)
        self.step.add_steps('phylo_tree')

    def phylo_tree_start_callback(self):
        self.step.phylo_tree.start()
        self.step.update()

    def phylo_tree_end_callback(self):
        self.step.phylo_tree.finish()
        self.step.update()

    def check_options(self):
        """
        检查参数是否正确
        """
        if not self.option("fasta_file").is_set:
            raise OptionError("请传入fasta序列文件", code="33200101")
        if self.option("method") not in ["mafft", "clustalw2"]:
            raise OptionError("请选择正确的比对方法", code="33200101")

    def set_resource(self):
        """
        所需资源
        """
        self._cpu = 10
        total = os.path.getsize(self.option("fasta_file").prop["path"])
        total = int(math.ceil(total / (1024 * 1024 * 1024)))
        total = int(total * 45)
        self._memory = "{}G".format(total)

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
            ["./phylo.tre", "tre", "进化树树文件"]
        ])
        #print self.get_upload_files()
        super(PaupPhylipAgent, self).end()


class PaupPhylipTool(Tool):
    """
    version 1.0
    """

    def __init__(self, config):
        super(PaupPhylipTool, self).__init__(config)
        self.clustalw2_path = self.config.SOFTWARE_DIR+'/bioinfo/align/clustalw-2.1/src/'
        self.python_path = '/miniconda2/bin/'
        self.mafft_path = self.config.SOFTWARE_DIR+'/bioinfo/align/mafft-7.299-with-extensions/bin/'
        #self.FastTree_path = os.path.join(self.config.SOFTWARE_DIR, "bioinfo/phylogenetic/fasttree2.1.9/FastTreeMP")
        self.paup = self.config.SOFTWARE_DIR + "/bioinfo/phylogenetic/paup4.0a166/paup4a166_centos64"
        self.phylip = self.config.SOFTWARE_DIR + "/bioinfo/compare_genome/software/phylip-3.697/exe"
        self.phylip_sh = self.config.PACKAGE_DIR + "/meta/phylip.sh"
        self.jmodetest = self.config.SOFTWARE_DIR + "/bioinfo/phylogenetic/jmodeltest-2.1.10/jModelTest.jar"

    def align(self):
        """
        比对，根据method参数，选择不同的比对软件进行比对，结果文件为phylo.align
        """
        if self.option("method") in ["mafft"]:
            cmd = "{}mafft {} > phylo.align".\
                format(self.mafft_path, self.option('fasta_file').prop['path'])
        else:
            cmd = self.clustalw2_path + "clustalw2 -ALIGN -INFILE=%s -OUTFILE=phylo.align  -OUTPUT=FASTA" % \
                                    self.option('fasta_file').prop['path']
        print cmd
        self.add_state('phylo_tree_start', data='开始运行程序生成树文件')
        # os.system(cmd)
        self.logger.info(cmd)
        self.logger.info("开始运行{}软件，进行比对".format(self.option("method")))
        command = subprocess.Popen(cmd, shell=True)
        command.communicate()
        if command.returncode == 0:
            self.logger.info("完成比对！")
        else:
            self.set_error("运行出错！", code="33200101")
            raise Exception("比对出错")
        # self.add_state('clustalw_end', data='done')

    def change_paup_format(self):
        params = "ToNEXUS format=FASTA toFile=align.nex fromFile=phylo.align;\n"
        with open("to_nex.params", 'w') as fw:
            fw.write(params)
        cmd = '{} -n < to_nex.params'.format(self.paup)
        self.logger.info("change_paup_format:" + cmd)
        command = subprocess.Popen(cmd, shell=True)
        command.communicate()
        self.logger.info('change_paup_format成功！')
        # if command.returncode == 0:
        #     self.logger.info("change_paup_format成功！")
        # else:
        #     self.set_error("change_paup_format运行出错！", code="33200101")

    def run_paup(self):
        params = 'HSearch;\nSaveTrees trees=firstOnly file=align.tre format=Newick brLens=yes;\nBootstrap nreps=%s;\nSaveTrees file=align_boot.tre format=Newick;\n'%self.option("bootstrap")
        with open('params.txt','w') as fw:
            fw.write(params)
        cmd = "{} -n align.nex  <params.txt".format(self.paup)
        self.logger.info("run_paup: " + cmd )
        command = subprocess.Popen(cmd, shell=True)
        command.communicate()
        if command.returncode == 0:
            self.logger.info("run_paup成功！")
        else:
            self.set_error("run_paup运行出错！", code="33200101")


    def get_phylip_format(self):
        cmd = "java -jar {} -d phylo.align -getPhylip".format(self.jmodetest)
        command = subprocess.Popen(cmd, shell=True)
        command.communicate()
        if command.returncode == 0:
            self.logger.info("生成Phylip文件phylo.align.phy")
        else:
            self.set_error("get_phylip_format运行出错！", code="33200101")

    def run_phylip(self):
        cmd = "{} {} {} {}".format(self.phylip_sh,'phylo.align.phy' ,self.phylip, self.option('bootstrap'))
        self.logger.info(cmd)
        command = subprocess.Popen(cmd, shell=True)
        command.communicate()
        if command.returncode == 0:
            self.logger.info("run_phylip成功，生成进化树！")
        else:
            self.set_error("run_phylip运行出错！", code="33200101")


    def set_output(self):
        """
        设置输出文件
        """
        self.logger.info("set out put")
        for f in os.listdir(self.output_dir):
            os.remove(os.path.join(self.output_dir, f))
        if self.option('tree_software') == 'phylip':
            os.link(self.work_dir+'/mp.bootstrap', self.output_dir+'/phylo.tre')
            self.option('phylo_tre').set_path(self.output_dir+'/phylo.tre')
        elif  self.option('tree_software') == 'paup':
            os.link(self.work_dir+"/align_boot.tre", self.output_dir+"/phylo.tre")
            self.option('phylo_tre').set_path(self.output_dir+'/phylo.tre')
        self.logger.info("done")

    def run(self):
        """
        运行
        """
        super(PaupPhylipTool, self).run()
        self.align()
        if self.option('tree_software') == 'paup':
            self.change_paup_format()
            self.run_paup()
        else:
            self.get_phylip_format()
            self.run_phylip()

        self.set_output()
        self.end()


if __name__ == '__main__':
    from mbio.workflows.single import SingleWorkflow
    from biocluster.wsheet import Sheet
    data={
        "name" : "phylo.paup_phylip",
        "id" : "test_paup_3",
        "type" : "tool",
        "options" : {
            "fasta_file" : "/mnt/ilustre/users/sanger-dev/sg-users/zouguanqing/tree_test/otu_reps_num4.fasta",
            "bootstrap" : 100,
            "tree_software" : 'paup'
        }
    }

    wsheet = Sheet(data=data)
    wf = SingleWorkflow(wsheet)
    wf.run()
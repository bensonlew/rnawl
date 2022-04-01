# -*- coding: utf-8 -*-
from biocluster.agent import Agent
from biocluster.tool import Tool
import os
import re
import shutil
from biocluster.config import Config
from biocluster.core.exceptions import OptionError
from operator import itemgetter

class PseudogeneAgent(Agent):
    """
    author: guanqing.zou
    last_modify: 20180528
    """
    def __init__(self, parent):
        super(PseudogeneAgent, self).__init__(parent)
        options = [
            # {"name": "query", "type": "infile", },
            {"name": "query", "type": "string", "default": ""},
            # {"name": "database", "type": "infile", }
            {"name": "ref", "type": "string", "default": ""},
            {"name": "sample", "type": "string","default":"sample"}
            ]
        self.add_option(options)
        self.step.add_steps('pseudogene')
        self.on('start', self.step_start)
        self.on('end', self.step_end)
        #self.queue = 'BLAST'  # 投递到指定的队列BLAST

    def step_start(self):
        self.step.pseudogene.start()
        self.step.update()

    def step_end(self):
        self.step.pseudogene.finish()
        self.step.update()

    def check_options(self):
        #if not self.option("bsnout").is_set:
        if not self.option("query"):
            raise OptionError("必须设置参数query", code="32201601")
        if not self.option("ref"):
            raise OptionError("必须设置参数ref", code="32201602")
        return True

    def set_resource(self):
        self._cpu = 32
        self._memory = '60G'

    def end(self):
        super(PseudogeneAgent, self).end()


class PseudogeneTool(Tool):
    def __init__(self, config):
        super(PseudogeneTool, self).__init__(config)

        self.perl_path = "program/perl-5.24.0/bin/perl"
        genewise_path = os.path.join(self.config.SOFTWARE_DIR,"bioinfo/Genomic/Sofware/wise-2.4.1/src/bin")
        self.set_environ(WISECONFIGDIR=self.config.SOFTWARE_DIR+'/bioinfo/Genomic/Sofware/wise-2.4.1/wisecfg')
        self.set_environ(PATH=genewise_path)
        #os.system("export WISECONFIGDIR={}".format(self.config.SOFTWARE_DIR+'/bioinfo/Genomic/Sofware/wise-2.4.1/wisecfg/'))
        #self.logger.info(os.system("echo $WISECONFIGDIR"))
        #self.logger.info("wiseconfigwiseconfig")
        
    def run_pseudogene(self):
        cmd_path = os.path.join(self.config.PACKAGE_DIR, "gene_structure/Runall_for_pseudogenes_identification.pl")
        genblast_path = os.path.join(self.config.SOFTWARE_DIR,"bioinfo/align/genBlast_v138_linux_x86_64/")
        genewise_path = os.path.join(self.config.SOFTWARE_DIR,"bioinfo/Genomic/Sofware/wise-2.4.1/src/bin/")
        self.result = os.path.join(self.work_dir,"pseudogene_information.xls")
        cmd_str = "{} {} -scaffold_file {} -proteins_file {} -proteins_dir Ref -genblast_dir Genblast -genewise_dir Genewise -genblast_bin {} -genewise_bin {}".format(self.perl_path,cmd_path, self.option('query'), self.option("ref"), genblast_path, genewise_path)

        self.logger.info("开始运行{}".format(cmd_str))
        command = self.add_command("pseudogene",cmd_str)
        command.run()
        self.wait(command)
        if command.return_code == 0:
                self.logger.info("{} 运行完成".format(cmd_str))
        else:
                self.set_error("%s 运行失败", variables=(cmd_str), code="32201601")

    def add_pse_id(self):
        new_file = os.path.join(self.output_dir,self.option('sample')+'_pseudogene.xls')
        fw = open(new_file,'w')
        fr = open(self.result)
        lines = fr.readlines()
        fw.write("Gene ID\tref\tScaffold ID\tBegin\tEnd\tFrameshifts\tPremature stop codons\n")

        plist = []
        pat_len = len('scaffold')
        for line in lines[1:]:
            spl = line.split('\t')
            spl.append(int(spl[1][pat_len:]))
            spl.append(int(spl[2]))
            plist.append(spl)
        s_plist = sorted(plist, key=itemgetter(-2,-1))
        id = 1
        for i in s_plist:
             pse_id = "Pse"+'0'*(4-len(str(id)))+str(id)
             fw.write(pse_id+'\t'+'\t'.join(i[:-2]))
             id += 1

    def run(self):
        """
        运行
        :return:
        """
        super(PseudogeneTool, self).run()
        self.run_pseudogene()
        self.add_pse_id()
        if os.path.exists(self.work_dir + '/Ref'):   #guanqing 20180718 删除过大的中间文件夹
            shutil.rmtree(self.work_dir + '/Ref')
        if os.path.exists(self.work_dir + '/Genblast'):
            shutil.rmtree(self.work_dir + '/Genblast')
        if os.path.exists(self.work_dir + '/Genewise'):
            shutil.rmtree(self.work_dir + '/Genewise')
        self.end()




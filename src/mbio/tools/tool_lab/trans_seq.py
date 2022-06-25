# -*- coding: utf-8 -*-
# __author__ = 'fwy'
import os
import shutil
from biocluster.core.exceptions import OptionError
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.config import Config
from mbio.packages.ref_rna.filter_gtf import FilterGtf
import re
import pandas as pd
import unittest


class TransSeqAgent(Agent):
    """
    Used for cds to protein code

    detail:
        table:
        OrderedDict([('Standard', '0'),
             ('Standard (with alternative initiation codons)', '1'),
             ('Vertebrate Mitochondrial', '2'),
             ('Yeast Mitochondrial', '3'),
             ('Mold, Protozoan,   Coelenterate Mitochondrial and Mycoplasma/Spiroplasma','4'),
             ('Invertebrate   Mitochondrial', '5'),
             ('Ciliate Macronuclear and  Dasycladacean', '6'),
             ('Echinoderm    Mitochondrial', '9'),
             ('Euplotid Nuclear', '10'),
             ('Bacterial', '11'),
             ('Alternative Yeast Nuclear', '12'),
             ('Ascidian Mitochondrial', '13'),
             ('Flatworm Mitochondrial', '14'),
             ('Blepharisma Macronuclear', '15'),
             ('Chlorophycean Mitochondrial', '16'),
             ('Trematode Mitochondrial', '21'),
             ('Scenedesmus obliquus', '22'),
             ('Thraustochytrium Mitochondrial', '23')])

    """
    def __init__(self, parent):
        super(TransSeqAgent, self).__init__(parent)
        options = [
            # 参考基因文件
            {"name": "fasta_file", "type": "infile", "format": "ref_rna_v2.fasta"},
            #'Frame(s) to translate,value for 1,2,3,F(1+2+3),-1,-2,-3,R(-1and-2and-3),6(F+R)'
            {"name": "frame", "type": "int","default":"1" },
            #Code to use with species of your choose
            {"name": "table", "type": "int", "default": "0"},
            # 翻译方向 ["5","3","both"]
            {"name": "direction", "type": "string", "default": "5"},
            # 翻译起始位点 ["1","2","3"]
            {"name": "initial_position", "type": "string", "default": "1"},
            # {"name": "regions", "type": "string"},
            # {"name": "trim", "type": "int", "default": "1"},
            # {"name": "clean", "type": "int", "default": "1"},

        ]
        self.add_option(options)
        self.step.add_steps("seq_trans")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.seq_trans.start()
        self.step.update()

    def stepfinish(self):
        self.step.seq_trans.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option('fasta_file').is_set:
            raise OptionError('必须输入参考基因组序列文件')
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 2
        self._memory = "5G"

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        super(TransSeqAgent, self).end()


class TransSeqTool(Tool):
    def __init__(self, config):
        super(TransSeqTool, self).__init__(config)
        self._version = "v1.0.1"
        self.python_path = '/miniconda2/bin/python'
        self.perl =  'miniconda2/bin/'
        self.tool_path=self.config.PACKAGE_DIR+"/tool_lab/trans_seq.py"
        self.transeq_path = self.config.SOFTWARE_DIR + "/bioinfo/tool_lab/transeq/"
        self.libnucleus_path=self.config.SOFTWARE_DIR+"/bioinfo/seq/EMBOSS-6.6.0/lib/"
        self.emboss_path=self.libnucleus_path=self.config.SOFTWARE_DIR+"/bioinfo/seq/EMBOSS-6.6.0/bin/"
        self.set_environ(PATH=self.emboss_path,LD_LIBRARY_PATH=self.libnucleus_path)


    def run(self):
        """
        运行
        :return:
        """
        super(TransSeqTool, self).run()
        transtype = {("5","1"):[1],("5","2"):[2],("5","3"):3,("5","both"):[1,2,3],("3","1"):[-1],("3","2"):[-2],("3","3"):[-3],("3","both"):[-1,-2,-3],
                     ("both", "1"):[-1,-2],("both", "2"):[2,-2],("both", "3"):[3,-3],("both", "both"):[1,2,3,-1,-2,-3]}
        type=transtype[(self.option("direction"),self.option('initial_position'))]
        self.run_seq_trans(type)
        self.set_output()
        self.end()

    def run_seq_trans(self,type):
        self.logger.info("开始对fasta文件进行序列翻译")
        table = self.option("table")
        name = {1:"_5_1",2:"_5_2",3:"_5_3",-1:"_3_1",-2:"_3_2",-3:"_3_3"}
        for i in type:
            outfa_dir=os.path.join(self.work_dir,"pep{}.fasta".format(name[i]))
            seq_trans_cmd = '{} {} --infa {} --oufa {} --table {} --trim {} '\
                           .format(self.python_path, self.tool_path, self.option("fasta_file").prop["path"],outfa_dir,table,i)
            commond_name="seq_trans{}".format(name[i])
            command = self.add_command(commond_name, seq_trans_cmd, ignore_error=True)
            command.run()
            self.wait(command)
            if command.return_code == 0:
                self.logger.info("运行seq_trans完成")
            else:
                self.set_error("运行seq_trans运行出错!")
                return False



    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        try:
                results = [i for i in os.listdir(self.work_dir) if i.endswith(".fasta") ]
                result_files = [os.path.join(self.work_dir, i) for i in results]
                link_files = [os.path.join(self.output_dir, i) for i in results]
                for i in link_files:
                    if os.path.exists(i):
                        os.remove(i)
                for n, i in enumerate(result_files):
                    os.link(result_files[n], link_files[n])
        except Exception as e:
            self.logger.info("设置结果目录失败{}".format(e))


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        test_dir = "/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/fungi/Saccharomyces_cerevisiae/Ensemble_release_39/"
        data = {
            "id": "seq_trans" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "tool_lab.trans_seq",
            "instant": True,
            "options": dict(
                # fasta_file=test_dir + "/" + "dna/Saccharomyces_cerevisiae.dna.toplevel.fa",
                fasta_file ="/mnt/ilustre/users/sanger-dev/sg-users/fuwenyao/tool_lab/transeq/input/Homo_sapiens.GRCh38.cds.all.fa" ,
                table="Standard",
                direction="both",
                initial_position="1",
                # min_len=500,
                # max_len=10000
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()

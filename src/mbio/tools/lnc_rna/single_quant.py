# -*- coding: utf-8 -*-
import os
import glob
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import unittest
__author__ = 'gdq'


class SingleQuantAgent(Agent):
    """
    Transcript/Gene Expression Quantification based on Salmon, RSEM, Kallisto
    """
    def __init__(self, parent):
        super(SingleQuantAgent, self).__init__(parent)
        options = [
            dict(name="transcriptome", type="infile", format="lnc_rna.fasta"),
            dict(name="fastq", type="string"),
            dict(name="method", type="string", default="rsem"),
            dict(name="libtype", type="string", default=None),
            dict(name="t2g", type="string"),
            dict(name="pool", type="int", default=7),
            dict(name="thread", type="int", default=6),
            dict(name="output", type="string", default=None),
            dict(name="read_len", type="int", default=149),
            dict(name="read_len_sd", type="int", default=30),
            dict(name="map_tool", type="string", default="bowtie2"),
        ]
        self.add_option(options)
        self._memory_increase_step = 50

    def check_options(self):
        if self.option("method").lower() not in ["rsem", "salmon", "kallisto"]:
            raise OptionError("Method is incorrect", code = "33707901")
        if not self.option("transcriptome").is_set:
            raise OptionError("transcriptome not exist", code = "33707902")
        if self.option("map_tool").lower() not in ["bowtie", "bowtie2", "star"]:
            raise OptionError("map_tool/aligner is incorrect", code = "33707903")
        if self.option("libtype") is not None:
            if self.option("libtype").lower() not in ["fr", "rf"]:
                raise OptionError("libtype argument is not in [None, 'fr', 'rf']", code = "33707904")

    def set_resource(self):
        self._cpu = self.option("thread")
        self._memory = "36G"
        self.logger.info("set_cpu: " + str(self._cpu)+'; '+"set_memory: "+str(self._memory))

    def end(self):
        super(SingleQuantAgent, self).end()


class SingleQuantTool(Tool):
    """
    RNAseq expression quantification tool, Rsem, Kallisto, Salmon are supported
    """
    def __init__(self, config):
        super(SingleQuantTool, self).__init__(config)
        software_dir = self.config.SOFTWARE_DIR
        self.python_path = 'program/Python/bin/python'
        # self.quant_toolbox = software_dir + '/bioinfo/rna/scripts/single_quant_toolbox.py'
        self.quant_toolbox = self.config.PACKAGE_DIR + '/denovo_rna_v2/single_quant_toolbox.py'
        self.rsem_path = software_dir + '/bioinfo/rna/RSEM-1.3.1/'
        self.kallisto_path = software_dir + '/bioinfo/rna/kallisto_linux-v0.43.1/'
        self.salmon_path = software_dir + '/bioinfo/rna/Salmon-0.8.2_linux_x86_64/bin/'
        self.bowtie2_path = software_dir + '/bioinfo/align/bowtie2-2.3.4.3-linux-x86_64/'
        self.bowtie_path = software_dir + '/bioinfo/align/bowtie-1.1.2/'
        self.star_path = software_dir + "/bioinfo/rna/star-2.5/bin/Linux_x86_64/"
        self.samtools_path = software_dir + "/program/Python/bin/samtools"
        self.gcc = software_dir + '/gcc/5.1.0/bin'
        self.gcc_lib = software_dir + '/gcc/5.1.0/lib64'
        self.perl = software_dir + '/program/perl/perls/perl-5.24.0/bin/'
        python_path = self.config.SOFTWARE_DIR + '/program/Python/bin/'
        self.set_environ(PATH=python_path)
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)
        self.set_environ(PATH=self.rsem_path)
        self.set_environ(PATH=self.perl)

    def quant(self):
        self.logger.info(self.perl)
        self.logger.info(os.environ["PATH"])
        rsem_mapper = self.option("map_tool")
        if rsem_mapper == "bowtie":
            mapper_path = self.bowtie_path
        elif rsem_mapper == "bowtie2":
            mapper_path = self.bowtie2_path
        else:
            mapper_path = self.star_path
        cmd = '{} {} '.format(self.python_path, self.quant_toolbox)
        cmd += '-t {} '.format(self.option("transcriptome").prop['path'])
        cmd += '-fq "{}" '.format(self.option("fastq"))
        cmd += '-m {} '.format(self.option("method"))
        cmd += '-t2g {} '.format(self.option("t2g"))
        if self.option("output") is None:
            self.option("output", self.work_dir)
        else:
            if not os.path.exists(self.option("output")):
                os.mkdir(self.option("output"))
        cmd += '-o {} '.format(self.option("output"))
        if self.option("libtype"):
            cmd += '-strand {} '.format(self.option("libtype"))
        cmd += '-pool {} '.format(self.option("pool"))
        cmd += '-thread {} '.format(self.option('thread'))
        cmd += '-salmon {} '.format(self.salmon_path)
        cmd += '-kallisto {} '.format(self.kallisto_path)
        cmd += '-rsem {} '.format(self.rsem_path)
        cmd += '-rl {} '.format(self.option("read_len"))
        cmd += '-sd {} '.format(self.option("read_len_sd"))
        cmd += '-mapper {} '.format(self.option("map_tool"))
        cmd += '-mapper_path {} '.format(mapper_path)
        command = self.add_command('quant', cmd, ignore_error=True)
        command.run()
        self.wait()
        if self.option("method")=="rsem":
            try:
                glob.glob(self.work_dir + "/*_quant/*.transcript.sorted.bam")[0]
            except:
                self.add_state("memory_limit", "memory is low!")
        if command.return_code == 0:
            self.logger.info("quant 运行成功")
        #elif command.return_code in [1, -9]:  # add memory limit by shicaiping at 20180724
            #self.add_state("memory_limit", "memory is low!")
        elif command.return_code is None:
            self.logger.warn("运行Quant出错，返回值为None，尝试重新运行")
            command.rerun()
            self.wait()
            if command.return_code is 0:
                self.logger.info("quant 运行成功")
            else:
                self.set_error("运行>>>%s出错", variables = (cmd), code = "33707905")
        else:
            self.set_error("运行>>>%s出错", variables = (cmd), code = "33707906")

    def set_output(self):
        pass

    def run(self):
        super(SingleQuantTool, self).run()
        self.quant()
        self.set_output()
        self.end()


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "Quant" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "lnc_rna.single_quant",
            "instant": False,
            "options": dict(
                transcriptome="/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/test_denovo/test_data2/Trinity.fasta",
                fastq="S1;/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/test_denovo/test_data2/rawdata/S1.1.fq;/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/test_denovo/test_data2/rawdata/S1.2.fq",
                method="rsem",
                t2g="/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/test_denovo/test_data2/t2g",
                pool=1,
                thread=6,
                output=None,
                read_len=149,
                read_len_sd=30,
                map_tool="bowtie2",
            )
        }
        data['options']['method'] = 'rsem'
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()
        #
        data['id'] += '1'
        data['options']['method'] = 'salmon'
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()
        #
        data['id'] += '2'
        data['options']['method'] = 'kallisto'
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()



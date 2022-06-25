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
    mRNA, sRNA, mRNAsRNA Expression Quantification based on Salmon, RSEM, Kallisto
    """
    def __init__(self, parent):
        super(SingleQuantAgent, self).__init__(parent)
        options = [
            dict(name="transcriptome", type="infile", format="prok_rna.fasta"),
            dict(name="fastq", type="string"),
            dict(name="method", type="string", default="rsem"),
            dict(name="libtype", type="string", default=None),
            dict(name="pool", type="int", default=7),
            dict(name="thread", type="int", default=6),
            dict(name="output", type="string", default=None),
            dict(name="read_len", type="int", default=149),
            dict(name="read_len_sd", type="int", default=30),
            dict(name="map_tool", type="string", default="bowtie2"),
        ]
        self.add_option(options)

    def check_options(self):
        if self.option("method").lower() not in ["rsem", "salmon", "kallisto"]:
            raise OptionError("method is incorrect", code = "35004501")
        if not self.option("transcriptome").is_set:
            raise OptionError("transcriptome not exist", code = "35004502")
        if self.option("map_tool").lower() not in ["bowtie", "bowtie2", "star"]:
            raise OptionError("map_tool/aligner is incorrect", code = "35004503")
        if self.option("libtype") is not None:
            if self.option("libtype").lower() not in ["fr", "rf"]:
                raise OptionError("libtype argument is not in [None, 'fr', 'rf']", code = "35004504")

    def set_resource(self):
        self._cpu = self.option("thread") + 1
        if self.option("method").lower() == "rsem":
            self._memory = "{}G".format(30)
        else:
            self._memory = "{}G".format(20)
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
        self.python_path = 'miniconda2/bin/python'
        self.quant_toolbox = self.config.PACKAGE_DIR + '/prok_rna/single_quant_toolbox.py'
        self.rsem_path = software_dir + '/bioinfo/ref_rna_v3/RSEM/miniconda3/bin/'
        #self.rsem_path = software_dir + '/bioinfo/rna/RSEM-1.3.1/'
        # self.rsem_path = software_dir + '/bioinfo/rna/RSEM-1.2.31/bin/'
        # self.kallisto_path = software_dir + '/bioinfo/rna/kallisto_linux-v0.43.1/'
        self.kallisto_path = software_dir + '/bioinfo/ref_rna_v2/kallisto/'
        # self.salmon_path = software_dir + '/bioinfo/rna/Salmon-0.8.2_linux_x86_64/bin/'
        self.salmon_path = software_dir + '/bioinfo/ref_rna_v2/salmon-latest_linux_x86_64/bin/'
        # self.bowtie2_path = software_dir + '/bioinfo/align/bowtie2-2.3.4.3-linux-x86_64/'
        self.bowtie2_path = software_dir + '/miniconda2/bin/'
        # self.bowtie_path = software_dir + '/bioinfo/align/bowtie-1.1.2/'
        self.bowtie_path = software_dir + '/miniconda2/bin/'
        # self.star_path = software_dir + "/bioinfo/rna/star-2.5/bin/Linux_x86_64/"
        self.star_path = software_dir + "/miniconda2/bin/"
        self.samtools_path = software_dir + "/miniconda2/bin/samtools"
        self.gcc = software_dir + '/gcc/5.1.0/bin'
        self.gcc_lib = software_dir + '/gcc/5.1.0/lib64'
        python_path = self.config.SOFTWARE_DIR + '/miniconda2/bin/'
        self.set_environ(PATH=python_path)
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)
        self.set_environ(PATH=self.rsem_path)
        self.perl = software_dir + '/miniconda2/bin/'
        self.set_environ(PATH=self.perl)

    def quant(self):
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
        command = self.add_command('quant', cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            if self.option("method").lower() == "rsem":
                try:
                    glob.glob(self.work_dir + "/*_quant/*.transcript.sorted.bam")[0]
                except:
                    self.set_error("运行Quant出错")
            self.logger.info("quant 运行成功")
            self.logger.info("quant 运行成功")
        elif command.return_code is None:
            self.logger.warn("运行Quant出错，返回值为None，尝试重新运行")
            command.rerun()
            self.wait()
            if command.return_code is 0:
                self.logger.info("quant 运行成功")
            else:
                self.set_error("Failed. >>> %s", variables = (cmd), code = "35004505")
        else:
            self.set_error("Failed. >>> %s", variables = (cmd), code = "35004506")

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
            "name": "prok_rna.single_quant",
            "instant": False,
            "options": dict(
                transcriptome="/mnt/ilustre/users/sanger-dev/i-sanger_workspace2/20191209/Prokrna_majorbio_226063/ExtractMrna2/gene.fa",
                fastq="MT_1;/mnt/ilustre/users/sanger-dev/i-sanger_workspace2/20191209/Prokrna_majorbio_226063/HiseqQc/output/sickle_dir/MT1_sickle_l.fastq;/mnt/ilustre/users/sanger-dev/i-sanger_workspace2/20191209/Prokrna_majorbio_226063/HiseqQc/output/sickle_dir/MT1_sickle_r.fastq",
                method="rsem",
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



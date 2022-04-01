# coding=utf-8
import os
import glob
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import unittest
# import pandas as pd
__author__ = 'gdq'


class TfPredictAgent(Agent):
    """
    tf_predict description
    """
    def __init__(self, parent):
        super(TfPredictAgent, self).__init__(parent)
        options = [
            {'name': 's', 'type': 'string', 'default': 'plant', 'format': 'None'},
            {'name': 'organism', 'type': 'string', 'default': 'unknown', 'format': 'None'},
            {'name': 'blast_all', 'type': 'string', 'default': 'yes', 'format': 'None'},
            {'name': 'hmmscan', 'type': 'string', 'default': '/mnt/ilustre/users/sanger-dev/sg-users/litangjian/hmm/hmmer-3.1b2-linux-intel-x86_64/binaries/hmmscan', 'format': 'None'},
            {'name': 'hmmdb', 'type': 'string', 'default': '/mnt/ilustre/users/sanger-dev/app/database/pfam_31/Pfam-A.hmm', 'format': 'None'},
            {'name': 'seqfile', 'type': 'infile', 'default': '/mnt/ilustre/users/sanger-dev/sg-users/deqing/TestFiles/TF/example.fa', 'format': 'sequence.fasta'},
            {'name': 'E', 'type': 'float', 'default': '0.001', 'format': 'None'},
            {'name': 'domE', 'type': 'float', 'default': '0.0001', 'format': 'None'},
            {'name': 'cpu', 'type': 'int', 'default': '12', 'format': 'None'},
            {'name': 'tfdb', 'type': 'string', 'default': '/mnt/ilustre/users/sanger-dev/app/database/TFDB/', 'format': 'None'},
            {'name': 'diamond', 'type': 'string', 'default': '/mnt/ilustre/users/sanger-dev/app/bioinfo/align/diamond-0.8.35/diamond', 'format': 'None'},
            {'name': 'evalue', 'type': 'float', 'default': '0.0001', 'format': 'None'},
        ]
        self.add_option(options)

    def check_options(self):
        if float(self.option('E')) == 0:
            raise OptionError('The E value is too low')

    def set_resource(self):
        self._cpu = int(self.option('cpu'))
        self._memory = "{}G".format('20G')

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "TfPredict"]
            ])
        """
        # more detail
        result_dir.add_regexp_rules([
            [r"*.xls", "xls", "xxx"],
            [r"*.list", "", "xxx"],
        ])
        """
        super(TfPredictAgent, self).end()


class TfPredictTool(Tool):
    """
    tf_predict description
    """
    def __init__(self, config):
        super(TfPredictTool, self).__init__(config)
        software_dir = self.config.SOFTWARE_DIR
        self.python_path = 'program/Python/bin/python'
        self.tf_predict = self.config.PACKAGE_DIR + '/transcription_factor/tf_predict.py'
        self.hmmscan = software_dir + '/bioinfo/align/hmmer-3.1b2-linux-intel-x86_64/binaries/hmmscan'
        self.diamond = software_dir + '/bioinfo/align/diamond-0.8.35/diamond'
        self.pfam_db = software_dir + "/database/Pfam/Pfam-A.hmm"  # hmm参考库
        self.tfdb = software_dir + "/database/TFDB/"
        self.gcc = software_dir + '/gcc/5.1.0/bin'
        self.gcc_lib = software_dir + '/gcc/5.1.0/lib64'
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)
        self.r_path = software_dir + "/program/R-3.3.1/bin:$PATH"
        self._r_home = software_dir + "/program/R-3.3.1/lib64/R/"
        self._LD_LIBRARY_PATH = software_dir + "/program/R-3.3.1/lib64/R/lib:$LD_LIBRARY_PATH"
        self.set_environ(PATH=self.r_path, R_HOME=self._r_home, LD_LIBRARY_PATH=self._LD_LIBRARY_PATH)

    def run_tf_predict(self):
        cmd = '{} {} '.format(self.python_path, self.tf_predict)
        cmd += '-{} {} '.format("s", self.option("s"))
        cmd += '-{} {} '.format("organism", self.option("organism"))
        cmd += '-{} {} '.format("blast_all", self.option("blast_all"))
        cmd += '-{} {} '.format("hmmscan", self.hmmscan)
        cmd += '-{} {} '.format("hmmdb", self.pfam_db)
        cmd += '-{} {} '.format("seqfile", self.option("seqfile").prop['path'])
        cmd += '-{} {} '.format("E", self.option("E"))
        cmd += '-{} {} '.format("domE", self.option("domE"))
        cmd += '-{} {} '.format("cpu", self.option("cpu"))
        cmd += '-{} {} '.format("tfdb", self.tfdb)
        cmd += '-{} {} '.format("diamond", self.diamond)
        cmd += '-{} {} '.format("evalue", self.option("evalue"))
        cmd_name = 'tf_predict'
        command = self.add_command(cmd_name, cmd)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("{} Finished successfully".format(cmd_name))
        elif command.return_code is None:
            self.logger.warn("{} Failed and returned None, we will try it again.".format(cmd_name))
            command.rerun()
            self.wait()
            if command.return_code is 0:
                self.logger.info("{} Finished successfully".format(cmd_name))
            else:
                self.set_error("{} Failed. >>>{}".format(cmd_name, cmd))
        else:
            self.set_error("{} Failed. >>>{}".format(cmd_name, cmd))

    def set_output(self):
        target_files = glob.glob(self.work_dir + '/domain_predict.txt')
        target_files += glob.glob(self.work_dir + '/final_tf_predict.xls')
        target_files += glob.glob(self.work_dir + '/predicted_TFs.fa')
        target_files += glob.glob(self.work_dir + '/diamond.out.txt')
        for each in target_files:
            name = os.path.basename(each)
            link = os.path.join(self.output_dir, name)
            if os.path.exists(link):
                os.remove(link)
            os.link(each, link)

    def run(self):
        super(TfPredictTool, self).run()
        self.run_tf_predict()
        self.set_output()
        self.end()


class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "TfPredict" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "rna.tf_predict",
            "instant": False,
            "options": dict(
                s="animal",
                organism="unknown",
                blast_all="yes",
                hmmscan="/mnt/ilustre/users/sanger-dev/sg-users/litangjian/hmm/hmmer-3.1b2-linux-intel-x86_64/binaries/hmmscan",
                hmmdb="/mnt/ilustre/users/sanger-dev/app/database/pfam_31/Pfam-A.hmm",
                seqfile="/mnt/ilustre/users/sanger-dev/sg-users/deqing/TestFiles/TF/example.fa",
                E="0.001",
                domE="0.0001",
                cpu="12",
                tfdb="/mnt/ilustre/users/sanger-dev/app/database/TFDB/",
                diamond="/mnt/ilustre/users/sanger-dev/app/bioinfo/align/diamond-0.8.35/diamond",
                evalue="0.0001",
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()
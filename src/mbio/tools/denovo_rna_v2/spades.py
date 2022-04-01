# coding=utf-8
import os
import glob
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import unittest
from Bio import SeqIO
# import pandas as pd
__author__ = '刘彬旭'


class SpadesAgent(Agent):
    """
    spades assemble
    """
    def __init__(self, parent):
        super(SpadesAgent, self).__init__(parent)
        options = [
            {'type': 'string', 'name': 'type', 'default': 'rna'},
            {'type': 'infile', 'name': 'l', 'format': 'denovo_rna_v2.common'},
            {'type': 'infile', 'name': 'r', 'format': 'denovo_rna_v2.common'},
            {'type': 'infile', 'name': 's', 'format': 'denovo_rna_v2.common'},
            {'type': 'string', 'name': 'o', 'default': 'spades_out'},
            {'type': 'int', 'name': 't', 'default': 20},
            {'type': 'int', 'name': 'm', 'default': 200},
            {'type': 'int', 'name': 'min_len', 'default': 200},
            {'type': 'string', 'name': 'k', 'default': "auto"},
            {'type': 'string', 'name': 'c', 'default': "auto"},
            {'type': 'int', 'name': 'Q', 'default': 'auto'},
            {'type': 'string', 'name': 'n', 'default': 'NODE'},
        ]
        self.add_option(options)

    def check_options(self):
        self._memory_increase_step = 200
        pass

    def set_resource(self):
        file_size = float(os.path.getsize(self.option('l').prop['path'])) / 1024 / 1024 /1024
        if file_size <= 10:
            mem = int(file_size) * 10 + 40 + 2
        elif file_size > 10 and file_size < 100:
            mem = int(file_size - 10) * 4 + 140
        else:
            mem = 500

        self._memory = "{}G".format(mem)
        self._cpu = self.option("t") + 1

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", ""]
            ])
        """
        # more detail
        result_dir.add_regexp_rules([
            [r"*.xls", "xls", "xxx"],
            [r"*.list", "", "xxx"],
            ])
        """
        super(SpadesAgent, self).end()


class SpadesTool(Tool):
    """
    spades assemble
    """
    def __init__(self, config):
        super(SpadesTool, self).__init__(config)
        software_dir = self.config.SOFTWARE_DIR
        self.python_path = 'program/Python/bin/python'
        self.spades = software_dir + '/bioinfo/denovo_rna_v2/SPAdes-3.13.1-Linux/bin/spades.py'
        self.gcc = software_dir + '/gcc/5.1.0/bin'
        self.gcc_lib = software_dir + '/gcc/5.1.0/lib64'
        self.set_environ(PATH=self.gcc, LD_LIBRARY_PATH=self.gcc_lib)
        self.r_path = software_dir + "/program/R-3.3.1/bin:$PATH"
        self._r_home = software_dir + "/program/R-3.3.1/lib64/R/"
        self._LD_LIBRARY_PATH = software_dir + "/program/R-3.3.1/lib64/R/lib:$LD_LIBRARY_PATH"
        self.set_environ(PATH=self.r_path, R_HOME=self._r_home, LD_LIBRARY_PATH=self._LD_LIBRARY_PATH)

    def get_memory(self):
        # get memory from sbatch
        with open(glob.glob("Spades.sbatch")[0], "r") as f:
            lines = f.readlines()
            for line in lines:
                if line.startswith("#SBATCH --mem="):
                    mem = line.strip().split("=")[1][:-1]
                    return mem


    def run_spades(self):
        try:
            memory_limit = self.get_memory()
        except:
            file_size = float(os.path.getsize(self.option('l').prop['path'])) / 1024 / 1024 /1024

            if file_size <= 10:
                mem = int(file_size) * 10 + 40 + 2
            elif file_size > 10 and file_size < 100:
                mem = int(file_size - 10) * 4 + 140
            else:
                mem = 500
            memory_limit = mem

        cmd = '{} {} '.format(self.python_path, self.spades)
        cmd += '--{} '.format(self.option("type"))
        if self.option("s").is_set:
            cmd += '-{} {} '.format("s", self.option("s").prop['path'])
        else:
            cmd += '-{} {} '.format("1", self.option("l").prop['path'])
            cmd += '-{} {} '.format("2", self.option("r").prop['path'])
        cmd += '-{} {} '.format("o", self.option("o"))
        cmd += '-{} {} '.format("t", self.option("t"))
        cmd += '-{} {} '.format("m", memory_limit)
        if self.option("k") != "auto":
            cmd += '-{} {} '.format("k", self.option("k"))
        if self.option("c") != "auto":
            cmd += '-{} {} '.format("-cov-cutoff", self.option("c"))
        if self.option("Q") != "auto":
            cmd += '-{} {} '.format("-phred-offset", self.option("Q"))

        if os.path.exists(self.work_dir + "/" + self.option("o")):
            cmd = '{} {} '.format(self.python_path, self.spades)
            cmd += "--restart-from last "
            cmd += '-{} {} '.format("o", self.option("o"))
            cmd += '-{} {} '.format("t", self.option("t"))
            cmd += '-{} {} '.format("m", memory_limit)

        cmd_name = 'spades'
        command = self.add_command(cmd_name, cmd, ignore_error=True)
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("{} Finished successfully".format(cmd_name))
        elif command.return_code == 6:
            # self.logger.info("内存不足，重新运行")
            self.set_error("内存不足，重新运行", code="32007001")
        else:
            self.logger.warn("{} Failed and returned None, we will try it again.".format(cmd_name))
            a = glob.glob("spades_out/K*")
            ks = [int(x[-2:]) for x in a]
            # ks_str = ",".join(map(str, sorted(ks)[:-1]))
            ks_str = str(sorted(ks)[0])

            if self.option("k") == "auto":
                cmd += ' -{} {}'.format("k", ks_str)
            else:
                cmd = cmd.replace("-k {}".format(self.option("k")), "-k 25")
            cmd_name = 'spades_rerun'
            command2 = self.add_command(cmd_name, cmd)
            command2.run()
            self.wait()
            if command2.return_code is 0:
                self.logger.info("{} Finished successfully".format(cmd_name))
            else:
                self.set_error("%s Failed. >>>%s", variables=(cmd_name, cmd), code="32007002")


    def rename(self):
        ass_fa = self.option("o") + '/hard_filtered_transcripts.fasta'
        with open(self.work_dir + '/assemble.fa', 'w') as fa_fo, open(self.work_dir + '/assemble.g2t', 'w') as g2t_fo:
            for seq in SeqIO.parse(ass_fa, "fasta"):
                new_name = self.option("n") + seq.name.split("NODE")[-1]
                if len(seq.seq) > self.option("min_len"):
                    fa_fo.write(">{}\n{}\n".format(new_name, seq.seq))
                    g2t_fo.write("{}\t{}\n".format(new_name, new_name))


    def set_output(self):
        pass
        '''Example:
        diff_files = glob.glob(self.option("output") + '/*_vs_*.xls')
        diff_list = glob.glob(self.option("output") + '/*.DE.list')
        diff_summary = glob.glob(self.option("output") + '/*summary.xls')
        all_files = diff_files + diff_list + diff_summary
        '''
        self.rename()
        all_files = [self.work_dir + '/assemble.fa', self.work_dir + '/assemble.g2t']
        for each in all_files:
            fname = os.path.basename(each)
            link = os.path.join(self.output_dir, fname)
            if os.path.exists(link):
                os.remove(link)
            os.link(each, link)


    def run(self):
        super(SpadesTool, self).run()
        self.run_spades()
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
        test_dir='/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/ref_denovo/'
        data = {
            "id": "spades" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "denovo_rna_v2.spades",
            "instant": False,
            "options": dict(
                l=test_dir + "C_0_5R_1_S211_L004_R1_001.fastq",
                r=test_dir + "C_0_5R_1_S211_L004_R2_001.fastq",
                n='sample1'
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()

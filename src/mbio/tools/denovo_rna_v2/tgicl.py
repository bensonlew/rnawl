# coding=utf-8
import os
import glob
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import unittest
from Bio import SeqIO
# import pandas as pd
__author__ = 'liubinxu'


class TgiclAgent(Agent):
    """
    merge assemble result
    """
    def __init__(self, parent):
        super(TgiclAgent, self).__init__(parent)
        options = [
            {'type': 'string', 'name': 'transcript' ,'default': ''},
            {'type': 'float', 'name': 'p', 'default': 88},
            {'type': 'int', 'name': 'c', 'default': 40}, # 该数量在1和16之间
            {'type': 'bool', 'name': 's', 'default': True},
        ]
        self.add_option(options)

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = 17
        self._memory = "{}G".format(self.option('c'))

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
        super(TgiclAgent, self).end()


class TgiclTool(Tool):
    """
    merge assemble result
    """
    def __init__(self, config):
        super(TgiclTool, self).__init__(config)
        software_dir = self.config.SOFTWARE_DIR
        self.perl = 'program/perl/perls/perl-5.24.0/bin/perl'
        self.tgicl = software_dir + '/bioinfo/denovo_rna_v2/TGICL-2.1/TGICL-2.1/blib/script/tgicl'
        self.tgicl_lib = software_dir + '/bioinfo/denovo_rna_v2/TGICL-2.1/TGICL-2.1/lib'
        self.tgicl_cfg = software_dir + '/bioinfo/denovo_rna_v2/TGICL-2.1/TGICL-2.1/conf/tgicl.cfg'
        self.set_environ(PERL5LIB=self.tgicl_lib)

    def merge_fasta(self):
        os.system("cat {} > {}".format(" ".join(self.option("transcript").split(",")), self.work_dir + "/all_fa"))

    def merge_result(self):
        seq_records = SeqIO.parse(self.work_dir + '/all_fa', 'fasta')
        diff_files = glob.glob(self.work_dir + '/asm*/contigs')
        ace_files = glob.glob(self.work_dir + '/asm*/ACE')
        with open(self.work_dir + "/merge.fa", 'w') as fo:
            for afile in diff_files:
                with open(afile, 'r') as fi:
                    fo.write(fi.read())

            '''
            # 不根据all_fa.singletons 判断没有合并的序列原因是有时候没有相关的文件出来
            if self.option('s') and os.path.exists(self.work_dir + '/all_fa.singletons'):
                with open(self.work_dir + '/all_fa.singletons', 'r') as fs:
                    sings = list()
                    for line in fs:
                        sings.append(line.strip())
            '''
            contigs_list = list()
            for afile in ace_files:
                with open(afile, 'r') as fi:
                    for line in fi:
                        if line.startswith("RD"):
                            cols = line.split()
                            gene = cols[1].split("_i")[0]
                            contigs_list.append(gene)

            for seq_record in seq_records:
                seq_seq = seq_record.seq
                seq_name = seq_record.name
                gene_id = seq_name.split("_i")[0]
                if gene_id in contigs_list:
                    pass
                else:
                    if "_cov_" in gene_id:
                        try:
                            print gene_id, seq_name
                            cov = float(gene_id.split("_cov_")[1].split("_")[0])
                            if cov < 20:
                                continue
                        except:
                            pass
                    line = '>{}\n{}\n'.format(seq_name, str(seq_seq))
                    fo.write(line)
                    contigs_list.append(gene_id)



    def run_tgicl(self):
        os.system("cp {} {}".format(self.tgicl_cfg, "tgicl.cfg"))
        '''
        if self.option("c") > 16:
            psx_cpu = 16
        else:
            psx_cpu = self.option("c")
        '''
        cmd = '{} {} '.format(self.perl, self.tgicl)
        cmd += '-{} {} '.format("F", self.work_dir + "/all_fa")
        cmd += '-{} {} '.format("p", self.option("p"))
        cmd += '-{} {} '.format("l", self.option("c"))
        cmd += '-{} {} '.format("c", 16)
        cmd_name = 'tgicl'
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
                self.set_error("%s Failed. >>>%s", variables=(cmd_name, cmd), code="32007101")
        else:
            self.set_error("%s Failed. >>>%s", variables=(cmd_name, cmd), code="32007102")

    def get_g2t(self):
        with open(self.work_dir + '/merge.g2t', 'w') as g2t_fo:
            for seq in SeqIO.parse('merge.fa', "fasta"):
                g2t_fo.write("{}\t{}\n".format(seq.name, seq.name))

    def set_output(self):
        pass
        '''Example:
        diff_files = glob.glob(self.option("output") + '/*_vs_*.xls')
        diff_list = glob.glob(self.option("output") + '/*.DE.list')
        diff_summary = glob.glob(self.option("output") + '/*summary.xls')
        all_files = diff_files + diff_list + diff_summary
        '''
        self.get_g2t()
        for each in ['merge.fa', 'merge.g2t']:
            fname = os.path.basename(each)
            link = os.path.join(self.output_dir, fname)
            if os.path.exists(link):
                os.remove(link)
            os.link(each, link)


    def run(self):
        super(TgiclTool, self).run()
        self.merge_fasta()
        self.run_tgicl()
        self.merge_result()
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
        test_dir='/mnt/ilustre/users/sanger-dev/sg-users/liubinxu/ref_denovo/denovo_result/sub_merge/'
        data = {
            "id": "Tgicl" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "denovo_rna_v2.tgicl",
            "instant": False,
            "options": dict(
                transcript=test_dir + "Trinity.sub1.fasta," + test_dir + "Trinity.sub2.fasta," + test_dir + "Trinity.sub3.fasta",
                p=88,
                c=10,
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()

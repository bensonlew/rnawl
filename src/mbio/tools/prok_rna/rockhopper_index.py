# coding=utf-8
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import unittest
import shutil
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
__author__ = 'fengyitong'


class RockhopperIndexAgent(Agent):
    """
    fasta files remove duplication
    """
    def __init__(self, parent):
        super(RockhopperIndexAgent, self).__init__(parent)
        options = [
            {'name': 'fna', 'type': 'string', 'default': ''},
            {'name': 'input_file', 'type': 'string'},
            {'name': 'type', 'type': 'string', 'default': 'gff'},
            {"name": "ptt", "type": "outfile", "format": "prok_rna.common"},
            {"name": "query", "type": "outfile", "format": "prok_rna.fasta"},
            {"name": "querypep", "type": "outfile", "format": "prok_rna.fasta"},
            {"name": "gtf", "type": "outfile", "format": "gene_structure.gtf"},
            {"name": "gtf_exon", "type": "outfile", "format": "gene_structure.gtf"},
        ]
        self.add_option(options)

    def check_options(self):
        pass

    def set_resource(self):
        self._cpu = 1
        self._memory = "{}G".format('10')

    def end(self):
        # result_dir = self.add_upload_dir(self.output_dir)
        # result_dir.add_relpath_rules([
        #     [".", "", ""]
        #     ])
        """
        # more detail
        result_dir.add_regexp_rules([
            [r"*.xls", "xls", "xxx"],
            [r"*.list", "", "xxx"],
            ])
        """
        super(RockhopperIndexAgent, self).end()

class RockhopperIndexTool(Tool):
    """
    fasta files remove duplication
    """
    def __init__(self, config):
        super(RockhopperIndexTool, self).__init__(config)
        self.python_path = self.config.SOFTWARE_DIR + '/program/Python/bin/python'
        self.python_path_faa = self.config.SOFTWARE_DIR + '/bioinfo/prok_rna/miniconda2/bin/python'
        self.rock_index = self.config.PACKAGE_DIR + "/prok_rna/rockhopper_index.py"
        self.bedtool_path = self.config.SOFTWARE_DIR + '/bioinfo/seq/bedtools-2.25.0/bin/bedtools'
        self.gff_read = self.config.SOFTWARE_DIR + '/bioinfo/rna/cufflinks-2.2.1/gffread'
        self.faa = self.config.PACKAGE_DIR + '/prok_rna/faa.py'

    def rockhopper(self):
        cmd = '{} {} '.format(self.python_path, self.rock_index)
        cmd += '-{} {} '.format("fna", self.option("fna"))
        cmd += '-{} {} '.format("input", self.option("input_file"))
        cmd += '-{} {} \n'.format("type", self.option("type"))
        cmd_name = 'rockhopper'
        command = self.add_command(cmd_name, cmd)
        command.software_dir = ""
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
                self.set_error("%s Failed. >>> %s", variables = (cmd_name, cmd), code = "35004001")
        else:
            self.set_error("%s Failed. >>> %s", variables = (cmd_name, cmd), code = "35004002")

    def get_fasta_from_ptt(self):
        chr_list = os.listdir(self.work_dir + '/rock_index')
        if os.path.exists(self.work_dir + '/Rockhopper_Results/'):
            shutil.rmtree(self.work_dir + '/Rockhopper_Results/')
            os.mkdir(self.work_dir + '/Rockhopper_Results/')
        else:
            os.mkdir(self.work_dir + '/Rockhopper_Results/')
        for i in chr_list:
            with open(self.work_dir + '/rock_index/' + i + '/' + i + '.ptt', 'r') as ptt_r, \
                open(self.work_dir + '/Rockhopper_Results/' + 'cds.bed', 'a') as bed_w, \
                open(self.work_dir + '/Rockhopper_Results/' + 'reshape.gtf', 'a') as gtf_w, \
                open(self.work_dir + '/Rockhopper_Results/' + 'reshape_exon.gtf', 'a') as gtf_exon, \
                open(self.work_dir + '/Rockhopper_Results/' + 'ptt.bed', 'a') as ptt_w:
                _ = ptt_r.readline()
                _ = ptt_r.readline()
                _ = ptt_r.readline()
                for line in ptt_r.readlines():
                    line = line.strip('\n').split('\t')
                    bed_w.write(
                        i + '\t' + str(int(line[0].split('..')[0]) - 1) + '\t' + line[0].split('..')[1] + '\t' + line[
                            5] + '\t0' + '\t' + line[1] + '\n')
                    ptt_w.write(i + '\t' + '\t'.join(line) + '\n')
                    gtf_w.write(i + '\t' + 'RefSeq\tCDS' + '\t' + line[0].split('..')[0] + '\t' + line[0].split('..')[1] + '\t' + '.' + '\t' + line[1] + '\t' + '0' + '\t' + 'transcript_id "' + line[5] + '"; gene_id "' + line[5] + '"; gene_name "' + line[4] + '";' + '\n')
                    gtf_exon.write(i + '\t' + 'RefSeq\texon' + '\t' + line[0].split('..')[0] + '\t' + line[0].split('..')[1] + '\t' + '.' + '\t' + line[1] + '\t' + '0' + '\t' + 'transcript_id "' + line[5] + '"; gene_id "' + line[5] + '"; gene_name "' + line[4] + '";' + '\n')
        # 使用samtools建索引， bedtools 强制建的索引在id后有空格可能出错
        if os.path.exists(self.option("fna") + ".fai"):
            pass
        else:
            cmd = "samtools faidx {}".format(self.option("fna"))
            os.system(cmd)
        cmd = "{bedtool_path} getfasta -fi {fna} -bed {work_dir}/Rockhopper_Results/cds.bed -s -name -fo {work_dir}/Rockhopper_Results/cds.fa".format(
            bedtool_path=self.bedtool_path, fna=self.option("fna"), work_dir = self.work_dir)
        self.logger.info("cmd {}".format(cmd))
        os.system(cmd)

    def run_faa(self):
        cmd = '{} {}'.format(self.python_path_faa, self.faa)
        cmd += ' -c {} -o {}'.format(self.work_dir +'/Rockhopper_Results/cds.fa', self.work_dir + '/Rockhopper_Results/cds.faa')
        cmd_name = 'run_faa'
        command = self.add_command(cmd_name, cmd, shell=True)
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
                self.set_error("%s Failed. >>> %s", variables = (cmd_name, cmd), code = "35004001")
        else:
            self.set_error("%s Failed. >>> %s", variables = (cmd_name, cmd), code = "35004002")

    def get_rna_seq(self):
        with open(self.option("input_file"), "r") as f, open("rrna.bed", "w") as w:
            rna_num = 0
            for line in f:
                cols = line.split("\t")
                if len(cols) == 9 and cols[2].lower() == "rrna":
                    if "16S " in cols[-1] or "23S " in cols[-1]:
                        rna_num += 1
                        if cols[-1].startswith("ID="):
                            gene_id = cols[-1].split(";")[0].strip().split("=")[-1]
                        else:
                            gene_id = cols[0] + "_" + cols[3] + "_" + cols[4]
                        w.write("\t".join([
                            cols[0],
                            str(int(cols[3]) - 1),
                            cols[4],
                            gene_id + "_" + str(rna_num),
                            "0",
                            cols[6]
                        ]) + "\n")

        self.getfasta_path  = "/bioinfo/rna/bedtools2-master/bin/bedtools"
        out_fa = self.output_dir + "/rrna.fa"
        if os.path.getsize("rrna.bed") > 0:
            cmd = "%s getfasta -fi %s -bed %s -fo %s -name -s" % (self.getfasta_path, self.option("fna"), "rrna.bed", out_fa)
            self.logger.info("开始打印cmd命令!")
            self.logger.info(cmd)
            fa_cmd = self.add_command("rrna_fa",cmd).run()
            self.wait(fa_cmd)
            if fa_cmd.return_code == 0:
                self.logger.info("%s运行完成" % fa_cmd)
            else:
                self.set_error("%s运行出错", variables = (fa_cmd), code = "33705602")
            print "生成rrna的fa序列"


    def get_gtf(self):
        if self.option('type') == 'gff':
            to_gtf_cmd = '%s %s -T -o %s  ' % (self.gff_read, self.option('input_file'), self.work_dir + '/Rockhopper_Results/ref.gtf')
        else:
            to_gtf_cmd = 'cp %s %s'%(self.option('input_file'), self.work_dir + '/Rockhopper_Results/ref.gtf')
        cmd_name = 'gffread'
        command = self.add_command(cmd_name, to_gtf_cmd)
        command.software_dir = ""
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
                self.set_error("{} Failed. >>>{}".format(cmd_name, to_gtf_cmd))
        else:
            self.set_error("{} Failed. >>>{}".format(cmd_name, to_gtf_cmd))

    def set_output(self):
        if os.path.exists(self.work_dir + '/Rockhopper_Results/reshape.gtf'):
            self.option("gtf").set_path(self.work_dir + '/Rockhopper_Results/reshape.gtf')
            self.option("gtf_exon").set_path(self.work_dir + '/Rockhopper_Results/reshape_exon.gtf')
            if os.path.exists(self.output_dir + '/reshape.gtf'):
                os.remove(self.output_dir + '/reshape.gtf')
            os.link(self.work_dir + '/Rockhopper_Results/reshape.gtf', self.output_dir + '/reshape.gtf')
        if os.path.exists(self.work_dir + '/Rockhopper_Results/cds.fa'):
            self.option("query").set_path(self.work_dir + '/Rockhopper_Results/cds.fa')
            if os.path.exists(self.output_dir + '/cds.fa'):
                os.remove(self.output_dir + '/cds.fa')
            os.link(self.work_dir + '/Rockhopper_Results/cds.fa', self.output_dir + '/cds.fa')
        if os.path.exists(self.work_dir + '/Rockhopper_Results/cds.faa'):
            self.option("querypep").set_path(self.work_dir + '/Rockhopper_Results/cds.faa')
            if os.path.exists(self.output_dir + '/cds.faa'):
                os.remove(self.output_dir + '/cds.faa')
            os.link(self.work_dir + '/Rockhopper_Results/cds.faa', self.output_dir + '/cds.faa')
        if os.path.exists(self.work_dir + '/Rockhopper_Results/ptt.bed'):
            self.option("ptt").set_path(self.work_dir + '/Rockhopper_Results/ptt.bed')
            if os.path.exists(self.output_dir + '/ptt.bed'):
                os.remove(self.output_dir + '/ptt.bed')
            os.link(self.work_dir + '/Rockhopper_Results/ptt.bed', self.output_dir + '/ptt.bed')

    def run(self):
        super(RockhopperIndexTool, self).run()
        self.rockhopper()
        self.get_fasta_from_ptt()
        self.run_faa()
        self.get_rna_seq()
        # self.get_gtf()
        self.set_output()
        self.end()

class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run this script to do test.
    """
    def test(self):
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        import datetime
        test_dir='/mnt/ilustre/users/sanger-dev/sg-users/fengyitong/prok_rna/pipline/ref'
        data = {
            "id": "Rockhopper_gff" + datetime.datetime.now().strftime('%H-%M-%S'),
            "type": "tool",
            "name": "prok_rna.rockhopper_index",
            "instant": False,
            "options": dict(
                fna = "/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/sanger_test/ref.fasta",
                input_file = "/mnt/ilustre/users/sanger-dev/sg-users/zoujiaxun/sanger_test/PROKKA_03222021.gff",
                type = "gff",
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()

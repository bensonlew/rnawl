# coding=utf-8
import os
import glob
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
import os
import unittest
from mbio.packages.rna.annot_config import AnnotConfig
__author__ = 'fengyitong'


class RockhopperAgent(Agent):
    """
    fasta files remove duplication
    """
    def __init__(self, parent):
        super(RockhopperAgent, self).__init__(parent)
        options = [
            {'name': 'fna', 'type': 'string', 'default': ''},
            {'name': 'input_file', 'type': 'string'},
            {'name': 'type', 'type': 'string', 'default': 'feature'},
            {'name': 'group_list', 'type': 'string', 'default': 'group_list'},
            {'name': 'trimPairFq', 'type': 'string', 'default': 'trimPairFq.list'},
            {'name': 'special', 'type': 'string', 'default': 'true'},
            {"name": "predict_fa", "type": "outfile", "format": "prok_rna.common"},
            {"name": "genome_bed", "type": "outfile", "format": "prok_rna.common"},
            {"name": "genome_fa", "type": "outfile", "format": "prok_rna.common"},
            {"name": "feature_fa", "type": "outfile", "format": "prok_rna.common"},
            {"name": "evalue", "type": "float", "default": 1e-5},  # evalue值,blast参数
        ]
        self.add_option(options)

    def check_options(self):
        pass

    def set_resource(self):
        with open(self.option('trimPairFq'), 'r') as fq:
            memory = len(fq.readlines()) * 20
        if memory < 80:
            memory = 80
        if memory > 180:
            memory = 300
        self._cpu = 20
        self._memory = "{}G".format(str(memory))

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
        super(RockhopperAgent, self).end()

class RockhopperTool(Tool):
    """
    fasta files remove duplication
    """
    def __init__(self, config):
        super(RockhopperTool, self).__init__(config)
        python_path = self.config.SOFTWARE_DIR + '/miniconda2/bin/'
        self.set_environ(PATH=python_path)
        self.python_rel = '/miniconda2/bin/python'
        self.python_path = self.config.SOFTWARE_DIR + '/miniconda2/bin/python'
        self.rock_index = self.config.PACKAGE_DIR + "/prok_rna/rockhopper_index.py"
        self.rock_run = self.config.PACKAGE_DIR + "/prok_rna/rockhopper_run.py"
        self.rock2bed = self.config.PACKAGE_DIR + "/prok_rna/rockhopper_to_bed.py"
        self.python_path_new = self.config.SOFTWARE_DIR + '/bioinfo/prok_rna/miniconda2/bin/python'
        self.script_file = self.config.SOFTWARE_DIR
        self.choose_predictrna = self.config.PACKAGE_DIR + "/prok_rna/rockhopper_choose_predict.py"

    def rockhopper(self):
        with open(self.option('trimPairFq'), 'r') as fq:
            # memory = len(fq.readlines()) * 20 - 2
            memory = len(fq.readlines()) * 20 - 10      # modified by zhangyitong on 20211220
        if memory < 80:
            memory = 80 - 2
        if memory > 180:
            memory = 300 - 10
        # tsg框架限制内存,临时修改
        if memory > 40:
            memory = 37
        cmd = '{} {} '.format(self.python_path, self.rock_index)
        cmd += '-{} {} '.format("fna", self.option("fna"))
        cmd += '-{} {} '.format("input", self.option("input_file"))
        cmd += '-{} {} \n'.format("type", self.option("type"))
        cmd += ' {} {} '.format(self.python_path, self.rock_run)
        cmd += '{} {} {} {} {}\n'.format(self.option("group_list"), self.option("trimPairFq"), self.option("special"), str(memory), self.work_dir)
        
        with open(self.work_dir + '/rockhopper.sh', 'w') as rock_sh:
            rock_sh.write(cmd)
        cmd = 'bash %s/rockhopper.sh' % self.work_dir
        cmd_name = 'rockhopper'
        command = self.add_command(cmd_name, cmd)
        command.software_dir = "/bin"
        command.run()
        self.wait()
        if command.return_code == 0:
            self.logger.info("{} Finished successfully".format(cmd_name))
        elif command.return_code in [1]:
            self.logger.info("return code: %s" % command.return_code)
            self.add_state('memory_limit', 'memory is low!')

        elif command.return_code is None:
            self.logger.warn("{} Failed and returned None, we will try it again.".format(cmd_name))
            command.rerun()
            self.wait()
            if command.return_code is 0:
                self.logger.info("{} Finished successfully".format(cmd_name))
            else:
                self.set_error("%s Failed. >>> %s", variables = (cmd_name, cmd), code = "35003901")
        else:
            self.set_error("%s Failed. >>> %s", variables = (cmd_name, cmd), code = "35003902")
        os.system("touch rockhopper.finished")

    def choose_predictRNA(self):
        cmd = '{} {} {} {}'.format(self.python_rel, self.choose_predictrna, self.option("fna"), self.script_file)
        choose_rna_command = self.add_command("extract_rna", cmd)
        choose_rna_command.run()
        self.wait()

    def run_diamond(self, db_name="bacteria"):
        self.cmd_path = "bioinfo/align/diamond-2.0.13"
        self.db_path = AnnotConfig().get_dmnd_path(db_name="bacteria",
                                                   version="2019",
                                                   nr_version="202110")
        hit_num = 1
        outfmt = "6 qseqid qlen sseqid slen evalue sallseqid stitle" # 表格
        query_name = os.path.join(self.work_dir, "Rockhopper_Results/genome.predicted_RNA.fa")
        cmd = os.path.join(self.cmd_path, "diamond")

        # outputfile = os.path.join(self.work_dir, "genome.predicted_cds.fa")
        
        outputfile = 'genome.predicted_RNA_vs_NR.txt'  # outfmt默认为5
        cmd += " {} -q {} -d {} -o {} -f {} -p {} -k {}".format(
            "blastx", "genome.predicted_RNA.fa", self.db_path,
            outputfile, outfmt, 20, hit_num)
        if self.option("evalue") != None:
            cmd += " -e {}".format(self.option("evalue"))

        cmd += " --more-sensitive"
        self.logger.info("开始运行blast")
        blast_command = self.add_command("diamond", cmd)
        blast_command.run()
        self.wait()
        if blast_command.return_code == 0:
            self.logger.info("运行diamond完成")
            self.logger.info(outputfile)
        elif blast_command.return_code == None:
            self.logger.info("重新运行diamond")
            blast_command.rerun()
            self.wait(blast_command)
            if blast_command.return_code == 0:
                self.logger.info("重新运行diamond成功")
            if blast_command.return_code in [1]:
                self.add_state("memory_limit")
            else:
                self.set_error("diamond运行出错!")
        else:
            self.set_error("diamond运行出错", code = "35000406")
        os.system("touch diamond.finished")


    def rockhopper_to_bed(self):
        cmd = '{} {} {} {} {}'.format(self.python_rel, self.rock2bed, self.option("fna"), self.script_file, 'genome.predicted_RNA_vs_NR.txt')
        blast_command = self.add_command("extract_rockhopper_result", cmd)
        blast_command.run()
        self.wait()


    def set_output(self):
        all_files = os.listdir(self.work_dir + '/Rockhopper_Results')
        all_files = [self.work_dir + '/Rockhopper_Results/' + each for each in all_files ]
        for each in all_files:
            if each.endswith('.fa') or each.endswith('.txt') or each.endswith('.bed') or each.endswith('.xls') or each.endswith('.faa') or each.endswith('.cds'):
                fname = os.path.basename(each)
                link = os.path.join(self.output_dir, fname)
                if os.path.exists(link):
                    os.remove(link)
                os.link(each, link)

        if os.path.exists(self.output_dir + '/genome.predicted_RNA.fa'):
            self.option("predict_fa").set_path(self.output_dir + '/genome.predicted_RNA.fa')
        if os.path.exists(self.output_dir + '/genome.gene.bed'):
            self.option("genome_bed").set_path(self.output_dir + '/genome.gene.bed')
        if os.path.exists(self.output_dir + '/genome.gene.fa'):
            self.option("genome_fa").set_path(self.output_dir + '/genome.gene.fa')
        if os.path.exists(self.output_dir + '/genome.feature.fa'):
            self.option("feature_fa").set_path(self.output_dir + '/genome.feature.fa')

    def run(self):
        super(RockhopperTool, self).run()
        if os.path.exists(self.work_dir + '/rockhopper.finished'):
            pass
        else:
            self.rockhopper()
        self.choose_predictRNA()
        if os.path.exists(self.work_dir + '/diamond.finished'):
            pass
        else:
            self.run_diamond()
        self.rockhopper_to_bed()
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
            "name": "prok_rna.rockhopper",
            "instant": False,
            "options": dict(
                fna = "/mnt/lustre/users/sanger/workspace/20181127/Prokrna_sanger_142069/remote_input/genome_db/P_87_scaf.fna",
                input_file = "/mnt/lustre/users/sanger/workspace/20181127/Prokrna_sanger_142069/remote_input/gff_or_gtf_file/ref_genome.gtf",
                type = "gtf",
                group_list = "/mnt/lustre/users/sanger/workspace/20181127/Prokrna_sanger_142069/remote_input/group_table/group.txt",
                trimPairFq = "/mnt/lustre/users/sanger/workspace/20181127/Prokrna_sanger_142069/HiseqQc/output/sickle_dir/fq_list.txt",
                special = 'true'
            )
           }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()


if __name__ == '__main__':
    unittest.main()

# -*- coding: utf-8 -*-
# __author__ = 'wangzhaoyue'
# last modified by shicaiping at 20180522
from biocluster.agent import Agent
from biocluster.tool import Tool
from biocluster.core.exceptions import OptionError
from mbio.packages.ref_rna.trans_step import merged_add_code
from mbio.packages.ref_rna.express.single_sample import *
import os
import unittest
import shutil

class NewTranscriptsAgent(Agent):
    """
    新转录本预测
    version v1.0.1
    """
    def __init__(self, parent):
        super(NewTranscriptsAgent, self).__init__(parent)
        options = [
            {"name": "tmap", "type": "infile", "format": "assembly.tmap"},  # compare后的tmap文件
            {"name": "ref_fa", "type": "infile", "format": "ref_rna_v2.fasta"},  # 参考基因文件
            {"name": "ref_gtf", "type": "infile", "format": "gene_structure.gtf"},  # 参考序列的注释文件
            {"name": "merged_gtf", "type": "string"},  # 输入文件，拼接后的注释文件
            {"name": "new_trans_gtf", "type": "outfile", "format": "ref_rna_v2.gtf"},  # 新转录本注释文件
            {"name": "new_genes_gtf", "type": "outfile", "format": "ref_rna_v2.gtf"},  # 新基因gtf文件
            {"name": "new_trans_fa", "type": "outfile", "format": "ref_rna_v2.fasta"},  # 新转录本注释文件
            {"name": "new_genes_fa", "type": "outfile", "format": "ref_rna_v2.fasta"},  # 新基因注释文件
            {"name": "add_code_gtf", "type": "outfile", "format": "ref_rna_v2.gtf"},  # 增加class_code后的注释文件
            {"name": "change_id_gtf", "type": "outfile", "format": "ref_rna_v2.gtf"},  # 转换ID后的注释文件
            {"name": "change_id_fa", "type": "outfile", "format": "ref_rna_v2.fasta"},  # 转换ID后的序列文件文件
            {"name": "ref_and_new_gtf", "type": "outfile", "format": "ref_rna_v2.gtf"},
            {"name": "trans2gene", "type": "outfile", "format": "ref_rna_v2.common"},
        ]
        self.add_option(options)
        self.step.add_steps("newtranscripts")
        self.on('start', self.stepstart)
        self.on('end', self.stepfinish)

    def stepstart(self):
        self.step.newtranscripts.start()
        self.step.update()

    def stepfinish(self):
        self.step.newtranscripts.finish()
        self.step.update()

    def check_options(self):
        """
        重写参数检测函数
        :return:
        """
        if not self.option('tmap'):
            raise OptionError('必须输入compare后tmap文件', code = "33703701")
        if not self.option('ref_fa'):
            raise OptionError('必须输入参考序列ref.fa', code = "33703702")
        if not self.option('ref_gtf'):
            raise OptionError('必须输入参考序列的注释文件', code = "33703703")
        if not self.option('merged_gtf'):
            raise OptionError('必须输入参考序列merged_gtf', code = "33703704")
        return True

    def set_resource(self):
        """
        设置所需资源，需在之类中重写此方法 self._cpu ,self._memory
        :return:
        """
        self._cpu = 1
        self._memory = "10G"

    def end(self):
        result_dir = self.add_upload_dir(self.output_dir)
        result_dir.add_relpath_rules([
            [".", "", "结果输出目录"],
        ])
        result_dir.add_regexp_rules([
            ["new_transcripts.gtf", "gtf", "新转录本注释文件"],
            ["new_genes.gtf", "gtf", "新基因注释文件"],
            ["new_transcripts.fa", "fa", "新转录本序列文件"],
            ["new_genes.fa", "fa", "新基因序列文件"],
        ])
        super(NewTranscriptsAgent, self).end()


class NewTranscriptsTool(Tool):
    def __init__(self, config):
        super(NewTranscriptsTool, self).__init__(config)
        self._version = "v1.0.1"
        self.perl_path = '/program/perl/perls/perl-5.24.0/bin/perl '
        self.Python_path = 'program/Python/bin/python '
        self.newtranscripts_gtf_path = self.config.SOFTWARE_DIR + '/bioinfo/rna/scripts/assembly_stat.py'
        self.change_id_path = self.config.SOFTWARE_DIR + '/bioinfo/rna/scripts/gtfmerge.pl'
        self.gffread_path = "bioinfo/rna/cufflinks-2.2.1/"

    def run(self):
        """
        运行
        :return:
        """
        super(NewTranscriptsTool, self).run()
        self.run_change_id()
        self.run_newtranscripts_gtf()
        self.run_change_newtranscripts_gtf()
        self.cat()
        self.run_newtranscripts_fa()
        self.trans2gene()
        self.run_alltranscripts_fa()
        self.set_output()
        self.end()

    def run_change_id(self):
        """
        运行perl脚本，将class_code为“=”全部替换掉
        """
        merged_gtf = self.work_dir + '/add_code_merged.gtf'
        merged_add_code(self.option('merged_gtf'), self.option('tmap').prop['path'], merged_gtf)
        cmd = self.perl_path + self.change_id_path + " -i %s -compare %s -ref %s -o %schange_id_merged.gtf" % (
            merged_gtf, self.option('tmap').prop['path'], self.option('ref_gtf').prop['path'], self.work_dir + "/")
        self.logger.info('运行perl，转换id')
        command = self.add_command("change_id_cmd", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("perl运行完成")
        else:
            self.set_error("perl运行出错!", code = "33703705")

    def run_newtranscripts_gtf(self):
        """
        运行python，挑出新转录本gtf文件(merged.gtf)
        """
        olddir = os.path.join(self.work_dir, 'merged')
        if not os.path.exists(olddir):
            os.mkdir(olddir)
        cmd = self.Python_path + self.newtranscripts_gtf_path \
            + " -tmapfile %s -transcript_file %s -out_new_trans %snew_transcripts.gtf -out_new_genes %snew_genes.gtf -out_old_trans %sold_trans.gtf -out_old_genes %sold_genes.gtf" % (
                self.option('tmap').prop['path'], self.work_dir + "/add_code_merged.gtf",
                self.work_dir+"/merged/", self.work_dir+"/merged/", self.work_dir+"/merged/", self.work_dir+"/merged/")
        self.logger.info('运行python，挑出新转录本gtf文件')
        command = self.add_command("merged_newtranscripts_gtf_cmd", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("python运行完成")
        else:
            self.set_error("python运行出错!", code = "33703706")

    def run_change_newtranscripts_gtf(self):
        """
        运行python，挑出新转录本gtf文件(change_id_merged.gtf)
        """
        newdir = os.path.join(self.work_dir, 'changed')
        if not os.path.exists(newdir):
            os.mkdir(newdir)
        cmd = self.Python_path + self.newtranscripts_gtf_path \
            + " -tmapfile %s -transcript_file %s -out_new_trans %snew_transcripts.gtf -out_new_genes %snew_genes.gtf -out_old_trans %sold_trans.gtf -out_old_genes %sold_genes.gtf" % (
                self.option('tmap').prop['path'], self.work_dir + "/change_id_merged.gtf",self.work_dir+"/changed/",
                self.work_dir+"/changed/", self.work_dir+"/changed/", self.work_dir+"/changed/")
        self.logger.info('运行python，挑出新转录本gtf文件')
        command = self.add_command("change_newtranscripts_gtf_cmd", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("python运行完成")
        else:
            self.set_error("python运行出错!", code = "33703707")

    def cat(self):
        cmd = "cat {} {} > ref_and_new.gtf".format(self.work_dir + "/changed/new_transcripts.gtf", self.option('ref_gtf').prop['path'])
        os.system(cmd)
        self.option("ref_and_new_gtf", self.work_dir + "/ref_and_new.gtf")

    def run_newtranscripts_fa(self):
        """
        提取新转录本fa序列，用于注释
        """
        cmd = self.gffread_path + "gffread %s -g %s -w new_transcripts.fa" % (
            self.work_dir + "/changed/" + "new_transcripts.gtf", self.option('ref_fa').prop['path'])
        self.logger.info('运行gtf_to_fasta，提取新转录本fa序列')
        command = self.add_command("newtranscripts_fa_cmd", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("gtf_to_fasta运行完成")
        else:
            self.set_error("gtf_to_fasta运行出错!", code = "33703708")

    def trans2gene(self):
        trans2gene_tmp = os.path.join(self.work_dir, 'trans2gene_tmp')
        trans2gene = os.path.join(self.work_dir, 'trans2gene')
        gtf(self.work_dir + "/ref_and_new.gtf", trans2gene_tmp)
        with open(trans2gene_tmp, "r") as r, open(trans2gene, "w") as w:
            for line in r:
                line = line.strip().split()
                w.write(line[1] + "\t" + line[0] + "\n")
        self.option('trans2gene', self.work_dir + "/trans2gene")
        self.logger.info("提取gene_transcript成功！")

    def run_alltranscripts_fa(self):
        """
        提取所有新转录本fa序列，用于定量
        """
        cmd = self.gffread_path + "gffread %s -g %s -w all_transcripts.fa" % (
            self.work_dir + "/ref_and_new.gtf", self.option('ref_fa').prop['path'])
        self.logger.info('运行gtf_to_fasta，提取所有新转录本fa序列')
        command = self.add_command("alltranscripts_fa_cmd", cmd).run()
        self.wait(command)
        if command.return_code == 0:
            self.logger.info("gtf_to_fasta运行完成")
        else:
            self.set_error("gtf_to_fasta运行出错!", code = "33703709")

    def set_output(self):
        """
        将结果文件复制到output文件夹下面
        :return:
        """
        self.logger.info("设置结果目录")
        try:
            if os.path.exists(self.output_dir + "/new_transcripts.gtf"):
                os.remove(self.output_dir + "/new_transcripts.gtf")
            if os.path.exists(self.output_dir + "/new_genes.gtf"):
                os.remove(self.output_dir + "/new_genes.gtf")
            if os.path.exists(self.output_dir + "/old_trans.gtf"):
                os.remove(self.output_dir + "/old_trans.gtf")
            if os.path.exists(self.output_dir + "/old_genes.gtf"):
                os.remove(self.output_dir + "/old_genes.gtf")
            if os.path.exists(self.output_dir + "/change_id_merged.gtf"):
                os.remove(self.output_dir + "/change_id_merged.gtf")
            if os.path.exists(self.output_dir + "/ref_and_new.gtf"):
                os.remove(self.output_dir + "/ref_and_new.gtf")
            if os.path.exists(self.output_dir + "/new_transcripts.fa"):
                os.remove(self.output_dir + "/new_transcripts.fa")
            if os.path.exists(self.output_dir + "/all_transcripts.fa"):
                os.remove(self.output_dir + "/all_transcripts.fa")
            if os.path.exists(self.output_dir + "/add_code_merged.gtf"):
                os.remove(self.output_dir + "/add_code_merged.gtf")
            if os.path.exists(self.output_dir + "/trans2gene"):
                os.remove(self.output_dir + "/trans2gene")
            os.link(self.work_dir + "/changed/new_transcripts.gtf", self.output_dir + "/new_transcripts.gtf")
            os.link(self.work_dir + "/changed/new_genes.gtf", self.output_dir + "/new_genes.gtf")
            os.link(self.work_dir + "/merged/old_trans.gtf", self.output_dir + "/old_trans.gtf")
            os.link(self.work_dir + "/merged/old_genes.gtf", self.output_dir + "/old_genes.gtf")
            os.link(self.work_dir + "/change_id_merged.gtf", self.output_dir + "/change_id_merged.gtf")
            os.link(self.work_dir + "/ref_and_new.gtf", self.output_dir + "/ref_and_new.gtf")
            os.link(self.work_dir + "/new_transcripts.fa", self.output_dir + "/new_transcripts.fa")
            os.link(self.work_dir + "/all_transcripts.fa", self.output_dir + "/all_transcripts.fa")
            os.link(self.work_dir + "/add_code_merged.gtf", self.output_dir + "/add_code_merged.gtf")
            os.link(self.work_dir + "/trans2gene", self.output_dir + "/trans2gene")
            self.option('new_trans_gtf').set_path(self.output_dir + "/new_transcripts.gtf")
            self.option('new_genes_gtf').set_path(self.output_dir + "/new_genes.gtf")
            self.option('change_id_gtf').set_path(self.output_dir + "/change_id_merged.gtf")
            self.option('ref_and_new_gtf').set_path(self.output_dir + "/ref_and_new.gtf")
            self.option('new_trans_fa').set_path(self.output_dir + "/new_transcripts.fa")
            self.option('trans2gene').set_path(self.output_dir + "/trans2gene")
            self.logger.info("设置拼接比较结果目录成功")

        except Exception as e:
            self.logger.info("设置拼接比较分析结果目录失败{}".format(e))
            self.set_error("设置拼接比较分析结果目录失败%s", variables = (e), code = "33703710")

class TestFunction(unittest.TestCase):
    """
    This is test for the tool. Just run script to do test.
    """
    def test(self):
        import random
        from mbio.workflows.single import SingleWorkflow
        from biocluster.wsheet import Sheet
        data = {
            "id": "NewTranscripts_" + str(random.randint(1, 10000)),
            "type": "tool",
            "name": "ref_rna_v2.assemble.new_transcripts",
            "instant": False,
            "options": dict(
                tmap="/mnt/ilustre/users/sanger-dev/workspace/20180514/Refrna_RefrnaV2_7320/RefrnaAssemble/StringtieMerge/output/cuffcmp.merged.gtf.tmap",
                merged_gtf="/mnt/ilustre/users/sanger-dev/workspace/20180514/Refrna_RefrnaV2_7320/RefrnaAssemble/StringtieMerge/output/merged.gtf",
                ref_gtf="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Mus_musculus/Ensemble_release_89/gtf/Mus_musculus.GRCm38.89.gtf",
                ref_fa="/mnt/ilustre/users/sanger-dev/app/database/Genome_DB_finish/vertebrates/Mus_musculus/Ensemble_release_89/dna/Mus_musculus.GRCm38.dna_rm.toplevel.clean.fa",
            )
        }
        wsheet = Sheet(data=data)
        wf = SingleWorkflow(wsheet)
        wf.run()

if __name__ == '__main__':
    unittest.main()